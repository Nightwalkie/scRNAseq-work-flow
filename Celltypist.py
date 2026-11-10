# 用已注释的单细胞参考训练 CellTypist 模型，基因与 Visium HD 数据取交集；人类基因符号。
required_package = {
    "numpy": "待核实",
    "pandas": "待核实",
    "scipy": "待核实",
    "scanpy": "待核实",
    "celltypist": "待核实",
}

import os

# 线程数按资源确认结果调整；须在导入 numpy 前设置。
os.environ["OMP_NUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["MKL_NUM_THREADS"] = "10"

import numpy as np
import pandas as pd
import scipy.io
import scanpy as sc
import celltypist

sc.settings.verbosity = 1

os.chdir(os.path.expanduser("~/VisiumHD_project"))
out_dir = "Results/Celltypist"
os.makedirs(out_dir, exist_ok=True)
ref_dir = os.path.expanduser("~/scRNAseq_project/Results/Annotation")   # Annotation.R 导出的训练数据

# 1. 读取单细胞参考：counts.mtx 为基因 × 细胞，meta.csv 含 barcode 和 cell_type 列
X = scipy.io.mmread(os.path.join(ref_dir, "counts.mtx")).T.tocsr().astype(np.float32)
genes = pd.read_csv(os.path.join(ref_dir, "genes.tsv"), header=None)[0].astype(str).values
barcodes = pd.read_csv(os.path.join(ref_dir, "barcodes.tsv"), header=None)[0].astype(str).values
meta = pd.read_csv(os.path.join(ref_dir, "meta.csv"), dtype={"barcode": str}).set_index("barcode")

ref = sc.AnnData(X=X, obs=meta.loc[barcodes].copy(), var=pd.DataFrame(index=genes))
ref.obs["cell_type"] = ref.obs["cell_type"].astype("category")
print(ref)

# 2. 与 HD 基因取交集；HD 为 .h5 时改用 sc.read_10x_h5()
hd_genes = sc.read_h5ad("Results/QC/adata_QC.h5ad", backed="r").var_names
shared = ref.var_names.intersection(hd_genes)
# 交集明显偏少时检查基因名格式：Seurat 会把 "_" 改成 "-"，重名基因后缀 ".1" 与 "-1" 不同
print(f"参考基因 {ref.n_vars} | HD 基因 {len(hd_genes)} | 交集 {len(shared)}")
ref = ref[:, shared].copy()

# 3. 查看各类型细胞数；少于 50 的类型由用户决定合并或删除，再继续
counts_per_type = ref.obs["cell_type"].value_counts()
print(counts_per_type)
print(counts_per_type[counts_per_type < 50])

# 每类最多随机抽 3000 个细胞
rng = np.random.default_rng(44)
idx = []
for pos in ref.obs.groupby("cell_type", observed=True).indices.values():
    idx.extend(pos if len(pos) <= 3000 else rng.choice(pos, 3000, replace=False))
ref = ref[np.sort(idx)].copy()
print(ref.obs["cell_type"].value_counts())

# 4. 标准化：CellTypist 要求每细胞总数 1e4 + log1p
sc.pp.normalize_total(ref, target_sum=1e4)
sc.pp.log1p(ref)

# 5. 训练：先训练一轮挑每类 top 基因，再用这些基因重训
model = celltypist.train(ref, labels="cell_type", n_jobs=10, feature_selection=True,
                         top_genes=300, max_iter=1000, check_expression=True)
# 不抽样、使用全部细胞时可改用 SGD mini-batch：
# model = celltypist.train(ref, labels="cell_type", n_jobs=10, use_SGD=True, mini_batch=True,
#                          balance_cell_type=True, feature_selection=True)
model.write(os.path.join(out_dir, "scRNA_celltypist.pkl"))
print(list(model.cell_types))
