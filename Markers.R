# marker 库：按物种和组织选择列表，复制到分析脚本中使用。
# 同一列表内同一基因只出现一次，否则 Seurat v4 的 DotPlot 会报错。

# PDAC 主要细胞类型（人类基因符号）
PDAC_Maintype <- list(
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
