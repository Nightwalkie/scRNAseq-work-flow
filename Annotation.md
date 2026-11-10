# 细胞注释：代码生成与执行 Prompt

你是一名单细胞分析智能体。请读取 [Annotation.R](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/Annotation.R)，基于 QC 模块保存的 final RDS，对各 cluster 进行手动细胞类型注释，按注释画 UMAP 和 marker DotPlot，并导出 CellTypist 训练数据。

- `mainType_marker` 可使用 [Markers.R](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/Markers.R) 中与物种、组织匹配的列表更换；模板中为 `PDAC_Maintype`。所选列表直接写入脚本，不用 `source()` 引入。
- 先画按 `seurat_clusters` 分组的 DotPlot 并查看 FindAllMarkers 结果（`p_val_adj < 0.05`），暂停；`anno` 由用户看图后给出，不擅自给 cluster 命名。`anno` 需覆盖全部 cluster。
- 导出 CellTypist 训练数据时使用 RNA assay 的原始 counts，不用 SCT；不应进入模型的类型（模板中为 Low quality）由用户确认。
- 版本基线、参数写法、注释风格、分阶段运行、依赖核对及输出说明等其余要求参考 [QC.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/QC.md)；文件位置遵循 [files.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/files.md)。
