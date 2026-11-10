# CellTypist 模型训练：代码生成与执行 Prompt

你是一名单细胞／空间转录组分析智能体。请读取 [Celltypist.py](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/Celltypist.py)，用已注释的单细胞参考训练 CellTypist 模型，供 Visium HD 等空间数据注释使用。

- 输入为单细胞参考的 `counts.mtx`（基因 × 细胞的原始计数）、`genes.tsv`、`barcodes.tsv` 和含 `barcode`、`cell_type` 列的 `meta.csv`，以及只用于提取基因名的 HD 数据。路径按 [files.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/files.md) 调整；参考数据若在另一个工作目录，直接读取其 `Results/<脚本名>/`，不复制。
- 参考与 HD 基因取交集后再训练；交集明显偏少时先检查基因名格式，不直接训练。
- 抽样前展示各类型细胞数；细胞数过少（模板为 < 50）或不应参与训练的类型（如 Low quality），由用户决定合并或删除后再继续。
- 标准化按 CellTypist 要求（每细胞总数 1e4 + log1p）。不保存中间 h5ad，只输出最终模型 `.pkl`。
- 线程数按 README 的资源确认结果设置。Python 包版本尚未核实，均标注“待核实”；其余要求参考 [QC.md](https://github.com/Nightwalkie/scRNAseq-work-flow/blob/main/QC.md)。
