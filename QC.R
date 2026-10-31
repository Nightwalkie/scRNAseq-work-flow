# scRNAseq_v4: Seurat 4.3.0.1 / DoubletFinder 2.0.3 / harmony 1.2.0
# 单个捕获文库的 10X RNA counts，使用人类基因符号。
library(Seurat)
library(DoubletFinder)
library(harmony)
library(ggplot2)
future::plan("sequential")
set.seed(44)

setwd("E:/Bioinformatics/Learning/R WorkFlow")
out_dir <- "results/QC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 兼容当前 irlba：fastpath = FALSE，PCA 用 20/30 维（50 维实测崩溃）。
# 以下局部适配也覆盖 DoubletFinder 内部的 PCA，不修改已安装的包。
df <- new.env(parent = asNamespace("DoubletFinder"))
df$RunPCA <- function(object, ...) Seurat::RunPCA(object, ..., fastpath = FALSE)
for (name in c("paramSweep_v3", "parallel_paramSweep_v3", "doubletFinder_v3")) {
  df[[name]] <- get(name, envir = asNamespace("DoubletFinder"))
  environment(df[[name]]) <- df
}

# 1. 读取、QC、去双胞
obj <- CreateSeuratObject(
  counts = Read10X("data/sample1"), project = "sample1",
  min.cells = 3, min.features = 200
)
obj$sample_id <- "sample1"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

ggsave(file.path(out_dir, "QC_before.png"),
       VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
               pt.size = 0, ncol = 4), width = 12, height = 4, dpi = 200)

# 阈值按组织和测序深度调整；小鼠基因前缀改为 ^mt- 和 ^Rp[sl]。
obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                nCount_RNA >= 500 & nCount_RNA <= 20000 &
                percent.mt <= 20 & percent.ribo <= 35)
ggsave(file.path(out_dir, "QC_after.png"),
       VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
               pt.size = 0, ncol = 4), width = 12, height = 4, dpi = 200)

# 临时聚类仅用于估计同型双胞比例。
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 20, fastpath = FALSE, seed.use = 44, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 44, verbose = FALSE)

sweep <- df$paramSweep_v3(obj, PCs = 1:20, sct = FALSE, num.cores = 1)
sweep_stats <- summarizeSweep(sweep, GT = FALSE)
pdf(file.path(out_dir, "pK.pdf"))
bcmvn <- find.pK(sweep_stats)
dev.off()
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
write.csv(bcmvn, file.path(out_dir, "pK.csv"), row.names = FALSE)

# 7.5% 是示例双胞率，按实际平台和上样量修改。
nExp <- round(round(0.075 * ncol(obj)) * (1 - modelHomotypic(obj$seurat_clusters)))
obj <- df$doubletFinder_v3(
  obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp,
  reuse.pANN = FALSE, sct = FALSE
)
obj$DFclass <- obj@meta.data[[grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)]]
write.csv(obj[[]], file.path(out_dir, "doublets.csv"))
obj <- subset(obj, subset = DFclass == "Singlet")

# 2. SCT
rm(sweep, sweep_stats)
obj <- SCTransform(obj, assay = "RNA", vst.flavor = "v2",
                   vars.to.regress = "percent.mt", seed.use = 44, verbose = FALSE)
DefaultAssay(obj) <- "SCT"

# 3. PCA；单批次直接使用 PCA 嵌入。
obj <- RunPCA(obj, assay = "SCT", npcs = 50, fastpath = FALSE, seed.use = 44, verbose = FALSE)
ggsave(file.path(out_dir, "PCA_elbow.png"), ElbowPlot(obj, ndims = 30),
       width = 6, height = 4, dpi = 200)
# 多批次时填入真实 sample_id，启用 Harmony，并把下方两处 reduction 改为 "harmony"。
# obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca",
#                   dims.use = 1:30, project.dim = FALSE)

# 4. 聚类和 UMAP
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 44, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 44, verbose = FALSE)
ggsave(file.path(out_dir, "UMAP_clusters.png"), DimPlot(obj, label = TRUE),
       width = 7, height = 5, dpi = 200)

# 5. 使用 RNA 表达量找各群 marker；沿用 sc_function.R 的 LR 检验。
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, assay = "RNA", slot = "data", test.use = "LR",
                          only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(markers, file.path(out_dir, "FindAllMarkers.csv"), row.names = FALSE)
write.csv(obj[[]], file.path(out_dir, "cell_metadata.csv"))
obj@misc$markers <- markers
obj@misc$session_info <- sessionInfo()
saveRDS(obj, file.path(out_dir, "scRNA_final.rds"))
