# 单细胞绘图：代码生成与执行 Prompt

你是一名单细胞分析智能体。请结合用户的绘图目标，读取 [plot.R](plot.R) 中对应的范例，并遵循 [files.md](files.md) 生成绘图代码。本模块包含 28 个绘图／导出函数和 1 个趋势聚类辅助函数，不包含 `dimplot_new`、`export_dimplot`，也不负责重新进行 QC、差异分析或拟时序推断。

## 1. 函数是模板的存储形式，应用时必须展开

**plot.R 用函数组织不同画法，方便检索。实际生成分析代码时，应提取所需函数体，放在分析流程中数据准备完成的位置，直接调用 Seurat、ggplot2、ComplexHeatmap 等包的绘图函数，并按实际数据修改。**

- 不把整份 plot.R 复制进分析脚本，不用 `source("plot.R")` 加模板函数调用作为最终实现，也不再定义 `my_dotplot()` 等业务封装。
- 同时展开所选函数内部的辅助逻辑，以及它调用的其他模板函数。例如，借鉴 `exp_distribution_ecdf` 时也需要读取 `ecdf_plot`；轨迹热图需要变量分组时，按需借鉴 `varhclust_for_hmap`。
- 展开后移除函数外壳、无关分支及 `return()`；函数内的 `on.exit()` 等收尾操作改为对应的顺序代码。按实际物种、数据类型和图形样式保留需要的部分，不复制一整套参数分支。
- 批量绘图可以使用简单的 `for` 循环。绘图接口要求的回调，例如 ComplexHeatmap 的 `cell_fun = function(...)`，可以保留；这不等于重新封装整段绘图流程。
- 固定参数直接写在使用处，不用 `pt_size`、`fig_width`、`top_n` 等变量集中保存。对象、基因集输入、结果表和重复使用的路径可以使用变量。
- 先参考 plot.R 中相应画法的默认参数，再按用户要求及实际数据调整；不要将无关的全部图形一次性生成。

例如，用户要求在 QC 后绘制指定 marker 点图，应在已有 `obj` 和分组信息后插入类似以下的代码，基因和分组字段按实际情况修改：

```r
p <- Seurat::DotPlot(
  obj, assay = "RNA", group.by = "seurat_clusters",
  features = list(Tcell = c("CD3D", "CD3E"), Bcell = c("CD79A", "MS4A1")),
  scale = TRUE, col.min = -1, col.max = 1
) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(file.path(out_dir, "marker_dotplot.png"), p,
               width = 10, height = 5, dpi = 300)
```

这里的 `out_dir` 由当前分析脚本确定：代码位于 `Script/QC.R` 时写入 `Results/QC/`；独立写成 `Script/plot.R` 时写入 `Results/plot/`。

## 2. 按用户目标选择画法

| 需求 | 参考函数 | 主要输入及说明 |
|---|---|---|
| QC 指标分布 | `qc_plot` | 含 RNA counts 的 Seurat 对象；四个 QC 指标直方图 |
| 基因或评分在嵌入上的分布 | `featureplot_new` | 已有 UMAP 等嵌入、基因名或 metadata 评分字段 |
| 多图拼接、导出 | `export_featureplot` | 图对象列表；拼接前需先完成各张图 |
| marker 点图 | `dotplot_marker` | Seurat 对象、分组字段、基因向量或命名基因列表 |
| 内置 marker 面板的表达图 | `featureplot_marker` | 已有嵌入的对象；内置人／鼠基因仅作为可修改的例子 |
| 各群 top marker 排序图 | `diff_loadingplot` | 包含 gene、cluster、avg_log2FC 或 logfc 的差异基因表 |
| 多组富集结果热图 | `enrichr_cluster_plot` | 包含 id、type、pvalue／adjustp 的已有富集结果 |
| 按组展示嵌入空间密度 | `prop_density`、`prop_density2` | 已有嵌入和分组字段；两种配色方案 |
| 两组组成对比 | `prop_back2back`、`prop_back2back_lollipop` | 分组和细胞类型；分别展示比例及比例的 log2 比值 |
| 基因表达分布比较 | `exp_distribution`、`exp_distribution_boxplot`、`exp_distribution_ecdf` | Seurat 对象、基因和分组；分别为均值／散点图、箱线图、累积分布图 |
| 通用累积分布图 | `ecdf_plot` | 包含 value、type 的长表 |
| top marker 平均表达热图 | `heatmap_diff_seurat` | Seurat 对象、差异基因表和分组字段 |
| 堆积细胞比例图 | `prop_plot_hca` | Var1 为样本或组，Var2 为细胞类型，Freq 为细胞数的长表 |
| 按样本比较细胞比例 | `prop_sample_compare` | 样本、实验组及细胞类型字段；先算各样本比例，再比较两个实验组 |
| 数值矩阵热图 | `heatmap_text` | 行列有名称的数值矩阵；标注绝对值超过阈值的数值 |
| 细胞组成环形图 | `seurat_circ_prop` | 分组字段和细胞类型 |
| 评分散点与误差线 | `sc_pointplot_stat` | value、type 长表；展示均值 ± SD 或 SEM，并可进行组间检验 |
| 亚群中的组间变化 | `diff_signal_acrossclust` | 一个基因或仅含一个基因集的命名列表、两个组及细胞类型；展示两组均值之差 |
| 趋势热图前的变量分组 | `varhclust_for_hmap` | 变量 × 细胞矩阵及已有参考值；返回 id、type 分组表，是辅助分析步骤 |
| 拟时序／连续状态热图 | `heatmap_correlation` | 变量 × 细胞矩阵、reference、varcluster；不生成拟时序 |
| 变量模块趋势曲线 | `hmap_cluster_vis` | 同上；展示各变量模块沿参考值变化的平滑曲线 |
| 单变量趋势曲线 | `hmap_variable_vis` | 表达／评分矩阵、reference 和所选变量；可另给 id、group 分组表 |
| 表达热图与比例圆圈 | `heatmap_dotplot_marker` | 对象、分组、命名 marker 列表或含 gene、term 的 CSV |
| 通用散点与均值线 | `common_boxplot3` | value、type 长表；保留原函数名，但图形实际是散点和均值线 |
| 双指标四象限密度图 | `seurat_cytometry` | metadata 中两个数值指标及可选分组字段 |

`reference` 应包含 `id`、`value` 两列，id 与矩阵列名对应；`varcluster` 应包含 `id`、`type`，id 与矩阵行名对应。普通 marker 展示不需要先运行轨迹或评分模块。

## 3. 与上游模块衔接

- 用户同时要求 QC 和画图时，先按 QC.md 完成相应阶段，再加入绘图代码。已有满足条件的对象或结果表时直接使用，不重复运行 SCT、Harmony、聚类、FindAllMarkers 或基因集评分。
- QC 图必须在 QC 的 `subset()` 前展示并按 QC.md 等待用户选择阈值；多样本的筛选前出图与筛选使用不同循环。
- PCA 维度选择仍遵循 QC.md 的 ElbowPlot 顺序，不因增加绘图模块而跳过暂停点。
- 绘制 UMAP 表达图需要已有 UMAP；绘制富集热图需要已有富集结果；绘制趋势图需要已有参考值。缺少上游结果时说明具体缺项，由总 prompt 调度对应分析模块，不能伪造分析结果或把细胞顺序当成推断出的拟时序。
- marker 图可以依据 cluster 绘制，尚未完成细胞类型注释时使用真实的 cluster 字段，不擅自补上细胞类型标签。

## 4. 借鉴时的注意事项

- 默认从 RNA assay 的标准化表达取值。QC 结束时默认 assay 可能是 SCT，因此要明确 assay；绘制已有评分时直接使用正确的 metadata 字段。分组字段必须使用对象中的真实列名。
- 内置 marker 包含特定组织的例子，不能当作所有组织通用的注释规则。按物种、组织和用户需求选择基因；报告实际缺失的目标基因。`dotplot_marker` 中自定义 marker 对应 `species = NULL` 分支，展开时直接保留自定义基因即可。
- 细胞、样本、基因和分组注释按 ID 对齐，不能只按行列位置拼接。热图的平均表达、表达比例和分组注释必须使用同一套行列顺序。
- 展示基因表达时区分原始 counts、标准化表达和跨组 Z-score。`heatmap_diff_seurat` 展示组内平均 RNA 标准化表达的行 Z-score；`heatmap_dotplot_marker` 使用 DotPlot 的缩放平均表达和表达比例，二者数值定义不同。
- FeaturePlot 的分位数截断仅用于色阶展示，不删除细胞，也不改写原始表达或评分。不同组需要比较颜色时使用一致的色阶。
- `diff_loadingplot` 的纵轴是 log2FC；`diff_signal_acrossclust` 计算的是两组均值之差，不能标成 log2FC。`prop_back2back_lollipop` 才计算两组比例的 log2 比值；有零比例时应明确处理方式，不能悄悄加伪计数后隐藏说明。
- `prop_sample_compare` 中每个点代表一个样本。其他表达／评分比较函数中每个点通常代表一个细胞；若目标是实验组之间的推断，应先确定统计单位、重复样本及配对关系，不能把细胞级检验当成样本级检验。
- 检验方法、分组顺序及多重比较处理按研究设计明确设置；图上的 SD、SEM、均值、P 值含义与实际计算保持一致。富集热图只展示输入结果，不重新计算或更改显著性结论。
- 趋势函数的 `nbin = 50` 沿用原代码，含义是每箱约 50 个细胞，实际箱数约为细胞数除以 nbin。两处趋势处理使用一致的分箱和平滑设置，保留足够的箱数。热图顶部的 0–100 标度表示排序后的相对位置，不是重新估计的时间距离。
- `heatmap_dotplot_marker` 中圆面积随表达细胞比例变化，`dot.scale = 1` 是本模板的默认尺寸倍率。四象限图中等于阈值的点归入低值侧；`median.thres = TRUE` 会逐组取中位数，跨组比较固定阈值时使用共同的明确边界。
- 每种图按实际分组数和基因数调整尺寸、标签和图例。不要为了改变图形结果而修改原始数据或人为改变细胞组成。

## 5. 依赖与兼容性

遵循 QC.md 的 Seurat v4 版本基线。plot.R 采用显式的 `包名::函数名`，定义函数时不加载全部依赖，不执行分析，也不自动安装包。只为实际选中的画法核对所需依赖。

2026-10-31 本机核对：R 4.1.3、Seurat 4.3.0.1、SeuratObject 4.1.3、future 1.31.0、Matrix 1.6-1.1、ggplot2 3.4.2。Matrix 的实际版本与 QC.md 的目标 1.6-1 存在差异，不能将二者写成同一版本。

当前缺失的直接或关键间接依赖如下。建议版本依据相应版本的包声明选择，**尚未在本机完成这些缺失包的安装及联合运行验证**；其他传递依赖也需核对，不能为安装绘图包而自动升级 R、Seurat、Matrix 或 ggplot2。

| 包 | 建议版本 | 用途及版本依据 |
|---|---|---|
| ComplexHeatmap | 2.10.0 | 热图、趋势变量聚类；[Bioconductor 3.14 版本手册](https://mirrors.dotsrc.org/bioconductor-releases/3.14/bioc/manuals/ComplexHeatmap/man/ComplexHeatmap.pdf) |
| circlize | 0.4.15 | 热图色阶；[版本声明](https://raw.githubusercontent.com/cran/circlize/0.4.15/DESCRIPTION) |
| ggpubr | 0.6.0 | 统计检验与图上标注；要求 ggplot2 ≥ 3.4.0；[版本声明](https://raw.githubusercontent.com/cran/ggpubr/0.6.0/DESCRIPTION) |
| rstatix | 0.7.2 | ggpubr 的统计依赖；[版本声明](https://raw.githubusercontent.com/cran/rstatix/0.7.2/DESCRIPTION) |
| ggbeeswarm | 0.7.2 | 分类散点；[版本声明](https://raw.githubusercontent.com/cran/ggbeeswarm/0.7.2/DESCRIPTION) |
| ggrastr | 1.0.2 | 表达分布图或启用 raster 的散点图；还依赖 Cairo、ragg 等；[版本声明](https://raw.githubusercontent.com/cran/ggrastr/1.0.2/DESCRIPTION) |
| signal | 0.7-7 | 仅选择 Savitzky–Golay 平滑时需要；默认 loess 不需要；[版本声明](https://raw.githubusercontent.com/cran/signal/0.7-7/DESCRIPTION) |

## 6. 保存与汇报

- 所有图和导出的表写入执行脚本对应的 `Results/<脚本名>/`，目录下不再创建子文件夹。图中多样本、多基因和不同画法用文件名区分。
- 纯绘图模块不保存新的 RDS 或 checkpoint。与 QC 联合执行时，只保留整体流程最终的完整 RDS。
- ggplot 图用 `ggsave()` 保存；基础 R 图先打开图形设备、绘制后 `dev.off()`；ComplexHeatmap 图用图形设备配合 `draw()`。已经执行 `draw()` 的模板在展开后只保留一次实际绘制，避免重复页面。
- 完成时说明使用了哪些画法、输入和分组字段、必要的模板调整及输出位置。对依赖缺失、未实际运行的分支和运行警告如实说明，不将语法检查表述为全部绘图验证通过。
