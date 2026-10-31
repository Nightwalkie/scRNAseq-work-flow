# 分析工作目录管理规范

生成和运行分析代码时，必须遵循以下目录结构。模板中的示例路径应按本规范和用户实际工作目录调整。

## 1. 工作目录

工作目录下只允许存在以下三个一级文件夹，名称与大小写固定：

- `Database`：存放输入数据。
- `Results`：存放各脚本生成的所有文件。
- `Script`：存放 R 或 Python 脚本。

分析使用的数据、脚本及产生的文件必须分别归入上述目录，不散放在工作目录根部，也不创建其他一级文件夹。

## 2. Database

- `Database` 下只允许存在以“数据集序列号／数据库登录号 + 可选备注”命名的文件夹。
- 例如：`GSE123456`、`GSE123456_肺组织`、`CRA001160`。
- 数据放在相应数据集文件夹内部，不直接放在 `Database` 根目录。
- 数据集文件夹内部可按原始数据组织方式保留样本目录及数据文件。

## 3. Results

- `Results` 下只允许存在以脚本文件名去掉扩展名后命名的文件夹。
- 例如，`QC.R` 的全部输出写入 `Results/QC/`；`Annotation.py` 的全部输出写入 `Results/Annotation/`。
- 图片、CSV、最终 RDS、日志及其他生成文件，都必须直接放在对应脚本的结果文件夹中。
- **脚本结果文件夹内不允许再创建任何子文件夹。** 例如，`Results/QC/` 内不得再按图片、表格或 RDS 分层。
- 多样本、多批次或不同类型的结果通过文件名区分，不通过新增子文件夹区分。
- `Results` 根目录不直接存放结果文件。

## 4. Script

- `Script` 下不允许创建任何子文件夹。
- 只允许直接存放 R 或 Python 脚本文件，即 `.R`、`.r` 或 `.py` 文件。
- 运行脚本时，将输出路径明确指向对应的 `Results/<脚本名>/`，避免把结果、日志或缓存写入 `Script`。

## 5. 示例

```text
工作目录/
├── Database/
│   └── GSE123456_肺组织/
│       ├── sample1/
│       │   ├── matrix.mtx.gz
│       │   ├── features.tsv.gz
│       │   └── barcodes.tsv.gz
│       └── sample2/
│           ├── matrix.mtx.gz
│           ├── features.tsv.gz
│           └── barcodes.tsv.gz
├── Results/
│   └── QC/
│       ├── sample1_QC_before.png
│       ├── sample2_QC_before.png
│       ├── PCA_elbow.png
│       ├── UMAP_clusters.png
│       ├── FindAllMarkers.csv
│       ├── cell_metadata.csv
│       └── scRNA_final.rds
└── Script/
    └── QC.R
```

本规范只规定文件位置，不要求额外生成文件。单细胞 QC 流程仍遵循 QC.md：不保存中间 checkpoint，只在最后保存一个完整的 final RDS。
