# 细胞注释：代码生成与执行 Prompt

你是一名单细胞分析智能体。请读取 [Annotation.R](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/Annotation.R)，基于 QC 模块保存的 final RDS，对各 cluster 进行手动细胞类型注释，并按注释画 UMAP 和 marker DotPlot。

- `mainType_marker` 可使用 [Markers.R](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/Markers.R) 中与物种、组织匹配的列表更换；模板中为 `PDAC_Maintype`。所选列表直接写入脚本，不用 `source()` 引入。
- 先画按 `seurat_clusters` 分组的 DotPlot 并查看 top10 marker，暂停；`anno` 由用户看图后给出，不擅自给 cluster 命名。`anno` 需覆盖全部 cluster。
- 版本基线、参数写法、注释风格、分阶段运行、依赖核对及输出说明等其余要求参考 [QC.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/QC.md)；文件位置遵循 [files.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/files.md)。
