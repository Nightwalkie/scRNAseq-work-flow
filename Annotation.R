# 读取 QC.R 的 final RDS，按 marker 手动注释 cluster；人类基因符号。
required_package <- c(
  Seurat = "4.3.0.1",
  SeuratObject = "4.1.3",
  ggplot2 = "3.4.2",
  dplyr = "待核实"
)

library(Seurat)
library(ggplot2)
library(dplyr)
set.seed(44)

setwd("~/scRNAseq_project")
out_dir <- "Results/Annotation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. 读取 QC 结果，查看各 cluster marker
obj <- readRDS("Results/QC/scRNA_final.rds")
Idents(obj) <- "seurat_clusters"
table(obj$seurat_clusters)

markers <- read.csv("Results/QC/FindAllMarkers.csv")
top10 <- markers[markers$p_val_adj < 0.05, ] %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
top10

# 大类 marker：来自 Markers.R 的 PDAC_Maintype，其他组织从 Markers.R 选择后替换。
mainType_marker <- list(
  T_cell        = c("CD3D", "CD3E"),
  CD8T          = c("CD8A", "CD8B"),
  CD4T          = c("CD4"),
  Treg          = c("FOXP3", "IL2RA", "IL7R"),
  NKcell        = c("KLRD1", "NKG7", "NCAM1"),
  B_cell        = c("CD19", "CD79A", "MS4A1"),
  Plasma        = c("IGKC", "IGLC2", "TNFRSF17"),
  Macrophage    = c("APOC1", "SPP1", "C1QC"),
  Monocyte      = c("FCN1", "S100A8", "S100A9"),
  cDC           = c("FCER1A", "CD207", "XCR1"),
  pDC           = c("IL3RA", "LILRA4"),
  mregDC        = c("LAMP3", "CCR7"),
  Mast          = c("KIT", "MS4A2"),
  Neutrophil    = c("FCGR3B", "CEACAM8", "CSF3R"),
  Fibroblast    = c("COL1A1", "COL1A2"),
  Endothelial   = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  Acinar        = c("PRSS1", "CPA1", "CTRB", "AMY2A", "PNLIP"),
  Proliferating = c("MKI67", "TOP2A"),
  Epithelial    = c("EPCAM", "KRT8", "KRT18", "KRT19", "SOX9"),
  Endocrine     = c("CHGA", "INS", "GCG", "SST", "PPY"),
  SMC           = c("MYH11", "MYOCD", "DES"),
  Pericyte      = c("RGS5", "NDUFA4L2", "HIGD1B")
)
setdiff(unlist(mainType_marker), rownames(obj[["RNA"]]))   # 缺失基因，下方自动去掉
gene_list <- lapply(mainType_marker, intersect, rownames(obj[["RNA"]]))
gene_list <- gene_list[lengths(gene_list) > 0]

p <- DotPlot(obj, assay = "RNA", group.by = "seurat_clusters", features = gene_list,
             scale = TRUE, col.min = -1, col.max = 1) +
  scale_color_gradientn(colors = colorRampPalette(c("#0070b2", "#009bc7", "#b8e3ea", "#f3f3f1",
                                                    "#fccdb9", "#f15e4c", "#da1735"))(100)) +
  theme_bw() +
  theme(legend.position = "right", legend.box = "vertical",
        legend.margin = margin(t = 0, unit = "cm"),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black")) +
  labs(x = "", y = "", title = "")
ggsave(file.path(out_dir, "marker_dotplot_seurat_clusters.png"), p,
       width = 4 + 0.35 * length(unlist(gene_list)), height = 2 + 0.4 * nlevels(obj$seurat_clusters),
       dpi = 300, limitsize = FALSE)

# 2. 手动注释：查看上方 DotPlot 和 top10 后填写。
# 每个 cluster 都要写，多个 cluster 可写同一类型；书写顺序即图中细胞类型顺序。以下为示例。
anno <- c(
  "0" = "Epithelial",
  "1" = "Endothelial",
  "2" = "Epithelial",
  "3" = "Epithelial",
  "4" = "Fibroblast",
  "5" = "Pericyte",
  "6" = "Fibroblast",
  "7" = "Macrophage",
  "8" = "Macrophage",
  "9" = "Epithelial",
  "10" = "SMC",
  "11" = "Acinar",
  "12" = "Endothelial",
  "13" = "T/NK",
  "14" = "Endothelial",
  "15" = "Epithelial",
  "16" = "B_cell",
  "17" = "Epithelial",
  "18" = "Plasma",
  "19" = "Monocyte",
  "20" = "Endocrine",
  "21" = "Pericyte",
  "22" = "Epithelial",
  "23" = "Low quality",
  "24" = "Epithelial",
  "25" = "Endothelial"
)
setdiff(levels(obj$seurat_clusters), names(anno))   # 应为 character(0)，否则有 cluster 漏写
obj$celltype <- factor(unname(anno[as.character(obj$seurat_clusters)]), levels = unique(anno))
table(obj$celltype)

# 3. 按注释画 UMAP；坐标轴与边框样式沿用 plot.R 的 featureplot_new。
p <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE,
             label.size = 4, pt.size = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text(colour = "black", size = 15),
        plot.title = element_text(size = 17),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        legend.text = element_text(size = 12, color = "black")) +
  labs(x = "UMAP1", y = "UMAP2", title = "")
ggsave(file.path(out_dir, "UMAP_celltype.png"), p, width = 8, height = 6, dpi = 300)

# 4. 按注释画自选基因 DotPlot；dot_genes 按需修改，同一基因不要出现在两个组里。
dot_genes <- list(
  Epithelial  = c("KRT19", "KRT7", "KRT18", "KRT8"),
  Endocrine   = c("CHGA", "INS", "GCG"),
  Endothelial = c("CLDN5", "VWF", "CDH5"),
  Fibroblast  = c("COL1A1", "DCN", "COL1A2"),
  Myeloid     = c("CD68", "LYZ", "ITGAM", "ITGAX"),
  "T/NK"      = c("CD3D", "CD3E", "NKG7", "GNLY"),
  Bcell       = c("MS4A1", "CD79A", "CD79B"),
  Plasma      = c("MZB1", "IGJ"),
  Mural       = c("NOTCH3", "CALD1", "MYL9"),
  "Mesenchymal-like" = c("VIM", "LGALS1"),
  Acinar      = c("PRSS1", "CPA1", "AMY2A", "PNLIP")
)
setdiff(unlist(dot_genes), rownames(obj[["RNA"]]))   # 缺失基因，下方自动去掉
gene_list <- lapply(dot_genes, intersect, rownames(obj[["RNA"]]))
gene_list <- gene_list[lengths(gene_list) > 0]

p <- DotPlot(obj, assay = "RNA", group.by = "celltype", features = gene_list,
             scale = TRUE, col.min = -1, col.max = 1) +
  scale_color_gradientn(colors = colorRampPalette(c("#0070b2", "#009bc7", "#b8e3ea", "#f3f3f1",
                                                    "#fccdb9", "#f15e4c", "#da1735"))(100)) +
  theme_bw() +
  theme(legend.position = "right", legend.box = "vertical",
        legend.margin = margin(t = 0, unit = "cm"),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black")) +
  labs(x = "", y = "", title = "")
ggsave(file.path(out_dir, "marker_dotplot_celltype.png"), p,
       width = 4 + 0.35 * length(unlist(gene_list)), height = 2 + 0.4 * nlevels(obj$celltype),
       dpi = 300, limitsize = FALSE)

saveRDS(obj, file.path(out_dir, "scRNA_annotated.rds"))
