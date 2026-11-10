# 单个捕获文库的 10X RNA counts，使用人类基因符号。
required_package <- c(
  Seurat = "4.3.0.1",
  SeuratObject = "4.1.3",
  future = "1.31.0",
  Matrix = "1.6-1",
  DoubletFinder = "2.0.3",
  harmony = "1.2.0",
  ggplot2 = "3.4.2",
  sctransform = "待核实",
  presto = "待核实"
)

library(Seurat)
library(DoubletFinder)
library(harmony)
library(ggplot2)
library(presto)  # 加速 FindAllMarkers 的 wilcox 检验
future::plan("sequential")
set.seed(44)

setwd("~/scRNAseq_project")
out_dir <- "Results/QC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 兼容当前 irlba：fastpath = FALSE
# 以下局部适配也覆盖 DoubletFinder 内部的 PCA，不修改已安装的包。
df <- new.env(parent = asNamespace("DoubletFinder"))
df$RunPCA <- function(object, ...) Seurat::RunPCA(object, ..., fastpath = FALSE)
for (name in c("paramSweep_v3", "parallel_paramSweep_v3", "doubletFinder_v3")) {
  df[[name]] <- get(name, envir = asNamespace("DoubletFinder"))
  environment(df[[name]]) <- df
}

# 1. 读取、QC、去双胞
obj <- CreateSeuratObject(
  counts = Read10X("Database/GSE123456/sample1"), project = "sample1",
  min.cells = 3, min.features = 200
)
obj$sample_id <- "sample1"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
obj$all_cells <- "All cells"

p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
             group.by = "all_cells", pt.size = 0, ncol = 2) &
  NoLegend() &
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path(out_dir, "QC_before.png"), p, width = 10, height = 8, dpi = 200)

# 阈值按组织和测序深度调整；小鼠基因前缀改为 ^mt- 和 ^Rp[sl]。
obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                nCount_RNA >= 500 & nCount_RNA <= 20000 &
                percent.mt <= 20 & percent.ribo <= 35)
p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
             group.by = "all_cells", pt.size = 0, ncol = 2) &
  NoLegend() &
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path(out_dir, "QC_after.png"), p, width = 10, height = 8, dpi = 200)

# 临时聚类仅用于估计同型双胞比例。
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 20, fastpath = FALSE, seed.use = 44, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 44, verbose = FALSE)

# num.cores 可根据实际内存和线程调整；每个进程会复制对象，Windows 下只能为 1。
sweep <- df$paramSweep_v3(obj, PCs = 1:20, sct = FALSE, num.cores = 6)
sweep_stats <- summarizeSweep(sweep, GT = FALSE)
bcmvn <- find.pK(sweep_stats)
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 7.5% 是示例双胞率，按实际平台和上样量修改。
nExp <- round(round(0.075 * ncol(obj)) * (1 - modelHomotypic(obj$seurat_clusters)))
obj <- df$doubletFinder_v3(
  obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp,
  reuse.pANN = FALSE, sct = FALSE
)
obj$DFclass <- obj@meta.data[[grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)]]
obj <- subset(obj, subset = DFclass == "Singlet")

# 2. SCT
rm(sweep, sweep_stats, bcmvn)
obj <- SCTransform(obj, assay = "RNA", vst.flavor = "v2",
                   vars.to.regress = "percent.mt", seed.use = 44, verbose = FALSE)
DefaultAssay(obj) <- "SCT"

# 3. PCA；单批次直接使用 PCA 嵌入。
obj <- RunPCA(obj, assay = "SCT", npcs = 50, fastpath = FALSE, seed.use = 44, verbose = FALSE)
# 多批次时填入真实 sample_id，启用 Harmony，并把下方两处 reduction 改为 "harmony"。
# obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca",
#                   dims.use = 1:30, project.dim = FALSE)

# 4. 聚类和 UMAP
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 44, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 44, verbose = FALSE)
ggsave(file.path(out_dir, "UMAP_clusters.png"), DimPlot(obj, label = TRUE),
       width = 7, height = 5, dpi = 200)

# 5. 使用 RNA 表达量找各群 marker；默认 wilcox 检验，presto 加速。
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, assay = "RNA", slot = "data",
                          only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(markers, file.path(out_dir, "FindAllMarkers.csv"), row.names = FALSE)
write.csv(obj[[]], file.path(out_dir, "cell_metadata.csv"))
obj@misc$markers <- markers
obj@misc$session_info <- sessionInfo()
saveRDS(obj, file.path(out_dir, "scRNA_final.rds"))
