# Visium HD（细胞分割）多样本：读取 → 合并 QC → CellTypist 注释 → 细胞类型叠 HE 图；人类基因符号。
required_package = {
    "scanpy": "1.12.5",
    "squidpy": "1.8.4",
    "anndata": "0.13.4",
    "celltypist": "1.7.1",
    "scikit-learn": "1.6.1",
    "geopandas": "1.2.0",
    "pandas": "3.1.1",
    "numpy": "2.5.4",
    "matplotlib": "3.11.2",
}

import os

# 线程数按资源确认结果调整；须在导入 numpy 前设置。
os.environ["OMP_NUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["MKL_NUM_THREADS"] = "10"

import gc
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import celltypist

sc.settings.n_jobs = 10
sc.settings.verbosity = 1

os.chdir(os.path.expanduser("~/VisiumHD_project"))
out_dir = "Results/QC_celltypist"
os.makedirs(out_dir, exist_ok=True)
db = "Database/GSE123456"

# 1. 读取：每个样本目录下有 segmented_outputs/（矩阵、geojson、scalefactors）和 spatial/（HE 图）
samples = sorted(d for d in os.listdir(db) if os.path.isdir(os.path.join(db, d, "segmented_outputs")))
print(samples)

adatas, spatial_uns = [], {}
for sid in samples:
    seg = os.path.join(db, sid, "segmented_outputs")
    ad = sc.read_10x_h5(os.path.join(seg, "filtered_feature_cell_matrix.h5"))
    ad.var_names_make_unique()

    # geojson 细胞多边形质心作为空间坐标（全分辨率像素）
    gj = gpd.read_file(os.path.join(seg, "cell_segmentations.geojson"))
    cen = gj.geometry.centroid
    coord = pd.DataFrame({"cx": cen.x.values, "cy": cen.y.values}, index=gj["cell_id"].values)
    ad.obs["cell_int"] = ad.obs_names.str.extract(r"cellid_0*(\d+)-")[0].astype(int).values
    mask = ad.obs["cell_int"].isin(coord.index).values
    print(f"{sid}: 矩阵 {ad.n_obs} 细胞, geojson {len(coord)} 细胞, 匹配 {mask.sum()}")
    ad = ad[mask].copy()
    ad.obsm["spatial"] = coord.loc[ad.obs["cell_int"].values, ["cx", "cy"]].values
    ad.obs["sample"] = sid
    ad.obs_names = sid + "_" + ad.obs_names

    # HE 图与 scalefactors；PNG 读入为 0–1 浮点，转 uint8 省内存
    with open(os.path.join(seg, "spatial", "scalefactors_json.json")) as f:
        sf = json.load(f)
    sf.setdefault("spot_diameter_fullres", 10)   # squidpy 叠图需要该字段，缺失时补为 10
    img = (plt.imread(os.path.join(db, sid, "spatial", "tissue_hires_image.png")) * 255).astype(np.uint8)
    spatial_uns[sid] = {"images": {"hires": img}, "scalefactors": sf}

    adatas.append(ad)
    del gj, cen, coord, ad, img
    gc.collect()

adata = sc.concat(adatas, join="inner", merge="same")
del adatas
gc.collect()
adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=samples)
adata.uns["spatial"] = spatial_uns
print(adata)

# 2. QC 指标；小鼠改为 mt- 和 ^Rp[sl]
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, inplace=True)
# 总 counts 为 0 的细胞百分比为 NaN，记为 0
adata.obs[["pct_counts_mt", "pct_counts_ribo"]] = adata.obs[["pct_counts_mt", "pct_counts_ribo"]].fillna(0)

feats = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
fig, axes = plt.subplots(2, 2, figsize=(4 + 1.2 * len(samples), 9))
for k, ax in zip(feats, axes.flatten()):
    sc.pl.violin(adata, k, groupby="sample", jitter=False, stripplot=False, rotation=45, ax=ax, show=False)
    ax.set_title(k); ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "QC_violin_before.png"), dpi=200, bbox_inches="tight")
plt.close()

# 3. 过滤：阈值看筛选前 QC 图后确定
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_cells(adata, max_genes=1000)
sc.pp.filter_cells(adata, max_counts=1000)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()
sc.pp.filter_genes(adata, min_cells=3)

fig, axes = plt.subplots(2, 2, figsize=(4 + 1.2 * len(samples), 9))
for k, ax in zip(feats, axes.flatten()):
    sc.pl.violin(adata, k, groupby="sample", jitter=False, stripplot=False, rotation=45, ax=ax, show=False)
    ax.set_title(k); ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "QC_violin_after.png"), dpi=200, bbox_inches="tight")
plt.close()

# 4. 保存原始 counts；与训练一致按 1e4 + log1p 标准化，存为 adata.raw
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# 5. CellTypist 注释：模型来自 Celltypist.py；逐样本预测以降低单次内存占用
model = celltypist.models.Model.load("Results/Celltypist/scRNA_celltypist.pkl")
print("模型类别:", list(model.cell_types))
print("模型基因与 HD 交集:", len(set(model.features) & set(adata.var_names)), "/", len(model.features))

adata.obs["celltypist_label"] = ""
adata.obs["celltypist_conf"] = np.nan
for sid in samples:
    m = (adata.obs["sample"] == sid).values
    ad = sc.AnnData(X=adata.X[m], obs=pd.DataFrame(index=adata.obs_names[m]),
                    var=pd.DataFrame(index=adata.var_names))
    pred = celltypist.annotate(ad, model=model, majority_voting=False)
    adata.obs.loc[ad.obs_names, "celltypist_label"] = pred.predicted_labels["predicted_labels"].astype(str).values
    adata.obs.loc[ad.obs_names, "celltypist_conf"] = pred.probability_matrix.max(axis=1).values
    print(f"{sid}: 注释完成 {ad.n_obs} 细胞")
    del ad, pred
    gc.collect()

# 置信度 < 0.5 记为 Unassigned
adata.obs["celltype"] = np.where(adata.obs["celltypist_conf"] >= 0.5, adata.obs["celltypist_label"], "Unassigned")
adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"], categories=list(model.cell_types) + ["Unassigned"])
adata.obs["celltype"] = adata.obs["celltype"].cat.remove_unused_categories()
print(adata.obs["celltype"].value_counts())

# 固定配色，所有样本一致；Unassigned 为浅灰
cats = list(adata.obs["celltype"].cat.categories)
pal = sc.pl.palettes.default_20 if len(cats) <= 20 else sc.pl.palettes.default_102
ct_colors = {c: pal[i] for i, c in enumerate([c for c in cats if c != "Unassigned"])}
ct_colors["Unassigned"] = "#d3d3d3"
adata.uns["celltype_colors"] = [ct_colors[c] for c in cats]

# 6. 每个样本：细胞类型叠 HE 图
for sid in samples:
    sq.pl.spatial_scatter(adata, library_key="sample", library_id=sid, img_res_key="hires",
                          color="celltype", shape=None, size=0.5, figsize=(10, 10), title=sid)
    plt.savefig(os.path.join(out_dir, f"{sid}_celltype_HE.png"), dpi=300, bbox_inches="tight")
    plt.close()

# 7. 保存最终对象：X 与 raw 为 log1p，layers["counts"] 为原始 counts
adata.write_h5ad(os.path.join(out_dir, "HD_annotated.h5ad"), compression="lzf")
