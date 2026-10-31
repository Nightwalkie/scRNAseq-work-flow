# 绘图模板索引：见 plot.md。28 个绘图/导出函数 + 1 个趋势聚类辅助函数。
# 不包含 DimPlot；不自动执行分析、安装依赖或保存 RDS。
# 函数仅用于组织范例；实际应用请展开所需函数体，固定参数写在调用处。
# 基线：Seurat 4.3.0.1 / SeuratObject 4.1.3 / ggplot2 3.4.2。
# 可选依赖仅在相应函数使用时加载；输出路径遵循 files.md。

# qc_plot：QC 指标直方图；Seurat 对象 ----
qc_plot <- function(data, species = "human") {
    counts <- Seurat::GetAssayData(object = data, slot = "counts", assay = "RNA")
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    if (species == "human") {
        mt.genes <- grep(pattern = "^MT-", x = rownames(counts), value = TRUE)
        rb.genes <- grep(pattern = "^RP[SL]", x = rownames(counts), value = TRUE)
    }
    else if (species == "mouse") {
        mt.genes <- grep(pattern = "^mt-", x = rownames(counts), value = TRUE)
        rb.genes <- grep(pattern = "^Rp[sl]", x = rownames(counts), value = TRUE)
    }
    data$percent_mito <- Matrix::colSums(counts[mt.genes, , drop = FALSE])/Matrix::colSums(counts) * 
        100
    data$percent_ribo <- Matrix::colSums(counts[rb.genes, , drop = FALSE])/Matrix::colSums(counts) * 
        100
    par(mfrow = c(2, 2))
    hist(data$nFeature_RNA, breaks = 100, prob = TRUE, main = "nFeature_RNA", xlab = "Detected genes")
    hist(data$nCount_RNA, breaks = 100, prob = TRUE, main = "nCount_RNA", xlab = "UMI counts")
    hist(data$percent_mito, breaks = 100, prob = TRUE, main = "Mitochondrial", xlab = "Percent")
    hist(data$percent_ribo, breaks = 100, prob = TRUE, main = "Ribosomal", xlab = "Percent")
}

# featureplot_new：基因或评分的嵌入分布图；返回单图或图列表 ----
featureplot_new <- function(data, reduction = "umap", pt.size = 1, color = "parula", features, raster = NULL, 
    outlier.rm = F) {
    Seurat::DefaultAssay(data) <- "RNA"
    if (color == "parula") {
        col <- colorRampPalette(c("#352A87", "#343EB1", "#1558D9", "#036CE0", "#107AD9", "#1387D3", "#0997D1", 
            "#06A4C9", "#0CADBB", "#23B4A9", "#43BB97", "#6ABE83", "#8FBF73", "#AEBD66", "#CBBB5B", "#E6B94F", 
            "#FDBE3C", "#FBCF2C", "#F5E21E", "#F9FB0E"))(100)
    }
    else if (color == "kovesi") {
        col <- colorRampPalette(c("#1B0084", "#1B099A", "#1C14AE", "#1C20C0", "#1F2DCE", "#233DD6", "#2950CE", 
            "#2E68AB", "#377989", "#3F876A", "#46954D", "#53A036", "#65AB26", "#7BB41A", "#95BE16", "#AFC61C", 
            "#C6CE26", "#DAD636", "#EAE04E", "#F1EC74", "#F8F7AF", "#FFFFFF"))(100)
    }
    else if (color == "ocean.haline") {
        col <- colorRampPalette(c("#2A186C", "#2E1B87", "#2B229F", "#1D35A1", "#11459A", "#0C5294", "#115D8F", 
            "#19678C", "#22718A", "#2B7A89", "#328388", "#388D88", "#3E9886", "#44A284", "#4CAC80", "#55B67A", 
            "#62BF72", "#73C869", "#8AD060", "#A4D75C", "#BEDC62", "#D5E171", "#EAE784", "#FDEF9A"))(100)
    }
    else if (color == "ocean.dense") {
        col <- colorRampPalette(c("#360E24", "#5D1A58", "#723692", "#795CC3", "#7487E0", "#7EB0E4", "#AAD2E3", 
            "#E6F1F1"))(100)
    }
    else if (color == "magma") {
        col <- c(colorRampPalette(c(viridis::magma(323, begin = 0.15)[1]))(10), viridis::magma(323, begin = 0.18))
    }
    else if (color == "white2blue") {
        col <- colorRampPalette(c("#FBF4F8", "#E5E0ED", "#BFC6DD", "#8CADCC", "#4E92BA", "#1871A8", "#085889", 
            "#003758"))(100)
    }
    else if (color == "blue2red") {
        col <- colorRampPalette(c("#0070b2", "#009bc7", "#5ec7dd", "#b8e3ea", "#f3f3f1", "#fccdb9", "#f79676", 
            "#f15e4c", "#da1735"))(100)
    }
    else if (color == "rainbow") {
        col <- colorRampPalette(c("#343391", "#0064af", "#0090cc", "#00b6db", "#01b7c2", "#53c0a3", "#8dcb8a", 
            "#bbd967", "#fbd324", "#f6bd25", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f", "#c02e2f", "#8b2a21"))(100)
    }
    if (length(features) > 1) {
        lapply(features, function(i) {
            Seurat::FeaturePlot(data, features = i, cols = col, reduction = reduction, ncol = 1, pt.size = pt.size, 
                min.cutoff = if (outlier.rm) 
                  "q5"
                else NA, max.cutoff = if (outlier.rm) 
                  "q95"
                else NA, label = F, order = T, raster = raster) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_text(colour = "black", 
                  size = 15), plot.title = ggplot2::element_text(size = 17), panel.border = ggplot2::element_rect(fill = NA, 
                  color = "black", linewidth = 1), legend.position = "none") + ggplot2::labs(x = "UMAP1", 
                y = "UMAP2", title = i)
        })
    }
    else {
        Seurat::FeaturePlot(data, features = features, cols = col, reduction = reduction, ncol = 1, pt.size = pt.size, 
            min.cutoff = if (outlier.rm) 
                "q5"
            else NA, max.cutoff = if (outlier.rm) 
                "q95"
            else NA, label = F, order = T, raster = raster) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
            axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), plot.title = ggplot2::element_text(size = 17), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1), legend.position = "none") + ggplot2::labs(x = "UMAP1", 
            y = "UMAP2")
    }
}

# export_featureplot：图列表拼接与导出 ----
export_featureplot <- function(all_plots, ncol = 3, dpi = 300, output = NULL) {
    plot <- cowplot::plot_grid(plotlist = all_plots, ncol = ncol)
    if (!is.null(output)) 
        ggplot2::ggsave(output, plot, dpi = dpi, width = ncol * 5.5, height = ceiling(length(all_plots)/ncol) * 
            5, limitsize = FALSE)
    plot
}

# dotplot_marker：marker 点图；自定义 marker 时 species = NULL ----
dotplot_marker <- function(data, group.by = "cluster", marker, species = "human", width = 15, height = 6, 
    output = NULL) {
    if (is.null(species)) {
        gene_list <- marker
    }
    else if (species == "human") {
        Tcell = c("CD3D", "CD3E")
        CD8T = c("CD8A", "CD8B")
        CD4T = c("CD4")
        Treg = c("FOXP3", "IL2RA", "IL7R")
        NKcell = c("KLRD1", "NKG7", "NCAM1")
        Bcell = c("CD19", "CD79A", "MS4A1")
        Plasma = c("IGKC", "IGLC2", "TNFRSF17")
        Macro = c("APOC1", "SPP1", "C1QC")
        Mono = c("FCN1", "S100A8", "S100A9")
        cDC = c("FCER1A", "CD207", "XCR1")
        pDC = c("IL3RA", "LILRA4")
        MAST = c("KIT", "MS4A2")
        Neu = c("FCGR3B", "CEACAM8", "CSF3R")
        Fibro = c("COL1A1", "COL1A2")
        Endo = c("PECAM1", "VWF")
        hepatocyte = c("ALB", "CYP3A4", "HNF4A")
        biliary = c("EPCAM", "KRT7", "KRT19")
        gene_list <- list(Tcell = Tcell, CD8T = CD8T, CD4T = CD4T, Treg = Treg, NKcell = NKcell, Bcell = Bcell, 
            Plasma = Plasma, Macro = Macro, Mono = Mono, cDC = cDC, pDC = pDC, MAST = MAST, Neu = Neu, 
            Fibro = Fibro, Endo = Endo, hepatocyte = hepatocyte, biliary = biliary)
    }
    else if (species == "mouse") {
        Tcell = c("Cd3d", "Cd3e")
        CD8T = c("Cd8a", "Cd8b1")
        CD4T = c("Cd4")
        Treg = c("Foxp3", "Il2ra", "Il7r")
        NKcell = c("Klrd1", "Nkg7", "Ncam1")
        Bcell = c("Cd19", "Cd79a", "Ms4a1")
        Plasma = c("Igkc", "Tnfrsf17")
        Macro = c("Apoc1", "Spp1", "C1qc")
        Mono = c("S100a8", "S100a9")
        cDC = c("Fcer1a", "Cd207", "Xcr1")
        pDC = c("Il3ra", "Gm14548")
        MAST = c("Kit", "Ms4a2")
        Neu = c("Csf3r", "Fut4")
        Fibro = c("Col1a1", "Col1a2")
        Endo = c("Pecam1", "Vwf")
        hepatocyte = c("Alb", "Cyp3a11", "Hnf4a")
        biliary = c("Epcam", "Krt7", "Krt19")
        gene_list <- list(Tcell = Tcell, CD8T = CD8T, CD4T = CD4T, Treg = Treg, NKcell = NKcell, Bcell = Bcell, 
            Plasma = Plasma, Macro = Macro, Mono = Mono, cDC = cDC, pDC = pDC, MAST = MAST, Neu = Neu, 
            Fibro = Fibro, Endo = Endo, hepatocyte = hepatocyte, biliary = biliary)
    }
    col <- colorRampPalette(c("#0070b2", "#009bc7", "#b8e3ea", "#f3f3f1", "#fccdb9", "#f15e4c", "#da1735"))(100)
    if (is.list(gene_list)) {
        gene_list <- lapply(gene_list, function(x) intersect(x, rownames(data[["RNA"]])))
        gene_list <- gene_list[lengths(gene_list) > 0]
    }
    else {
        gene_list <- intersect(gene_list, rownames(data[["RNA"]]))
    }
    plot <- Seurat::DotPlot(data, assay = "RNA", scale = T, col.min = -1, group.by = group.by, col.max = 1, 
        features = gene_list) + ggplot2::scale_color_gradientn(colors = col) + ggplot2::theme_bw() + 
        ggplot2::theme(legend.position = "right", legend.box = "vertical", legend.margin = ggplot2::margin(t = 0, 
            unit = "cm"), panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(color = "black", 
            size = 12), legend.text = ggplot2::element_text(size = 12, color = "black"), legend.title = ggplot2::element_text(size = 12, 
            color = "black")) + ggplot2::labs(x = "", y = "", title = "") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
        hjust = 1))
    if (is.null(output)) {
        plot
    }
    else {
        ggplot2::ggsave(output, plot, width = width, height = height)
    }
}

# featureplot_marker：内置人/鼠 marker 面板；按组织修改基因 ----
featureplot_marker <- function(data, reduction = "umap", pt.size = 0.5, outlier.rm = F, species = "human", 
    output = NULL) {
    Seurat::DefaultAssay(data) <- "RNA"
    col <- colorRampPalette(c("#0070b2", "#009bc7", "#5ec7dd", "#b8e3ea", "#f3f3f1", "#fccdb9", "#f79676", 
        "#f15e4c", "#da1735"))(100)
    if (species == "human") {
        select <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "FOXP3", "IL2RA", "IL7R", "KLRD1", "NKG7", 
            "NCAM1", "CD19", "CD79A", "MS4A1", "IGKC", "IGLC2", "TNFRSF17", "APOC1", "SPP1", "C1QC", 
            "FCN1", "S100A8", "S100A9", "FCER1A", "CD207", "XCR1", "IL3RA", "LILRA4", "KIT", "MS4A2", 
            "FCGR3B", "CEACAM8", "CSF3R", "COL1A1", "COL1A2", "PECAM1", "VWF", "ALB", "CYP3A4", "HNF4A", 
            "EPCAM", "KRT7", "KRT19")
    }
    else if (species == "mouse") {
        select <- c("Cd3d", "Cd3e", "Cd8a", "Cd8b1", "Cd4", "Foxp3", "Il2ra", "Il7r", "Klrd1", "Nkg7", 
            "Ncam1", "Cd19", "Cd79a", "Ms4a1", "Igkc", "Tnfrsf17", "Apoc1", "Spp1", "C1qc", "S100a8", 
            "S100a9", "Fcer1a", "Cd207", "Xcr1", "Il3ra", "Gm14548", "Kit", "Ms4a2", "Csf3r", "Fut4", 
            "Col1a1", "Col1a2", "Pecam1", "Vwf", "Alb", "Cyp3a11", "Hnf4a", "Epcam", "Krt7", "Krt19")
    }
    select <- intersect(select, rownames(data))
    all_plots <- lapply(select, function(i) {
        Seurat::FeaturePlot(data, features = i, cols = col, reduction = reduction, ncol = 1, pt.size = pt.size, 
            min.cutoff = if (outlier.rm) 
                "q5"
            else NA, max.cutoff = if (outlier.rm) 
                "q95"
            else NA, label = F, order = T) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
            axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), plot.title = ggplot2::element_text(size = 17), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1), legend.position = "none") + ggplot2::labs(x = "UMAP1", 
            y = "UMAP2", title = i)
    })
    height = round(length(all_plots)/5 + 0.35, 0) * 5
    all_plots <- cowplot::plot_grid(plotlist = all_plots, ncol = 5)
    if (is.null(output)) 
        return(all_plots)
    ggplot2::ggsave(output, all_plots, dpi = 300, width = 25, height = height, limitsize = FALSE)
}

# diff_loadingplot：各群 top marker 的 log2FC 排序图 ----
diff_loadingplot <- function(diff, topn = 15, color = "#560047", output = NULL) {
    if ("avg_log2FC" %in% names(diff)) 
        diff$logfc <- diff$avg_log2FC
    all_plots <- lapply(unique(diff$cluster), function(i) {
        input <- diff[diff$cluster == i, ]
        input <- data.frame(id = input$gene, value = input$logfc)
        input <- na.omit(input)
        input <- input[order(input$value), ]
        input$number <- 1:nrow(input)
        input$number <- input$number/nrow(input)
        input <- input[order(input$value, decreasing = T), ]
        select <- head(input[order(input$value, decreasing = TRUE), ], topn)
        y_min <- min(input$value)
        y_max <- max(input$value)
        top_genes <- select$id
        ggplot2::ggplot(input, ggplot2::aes(x = number, y = value)) + ggplot2::geom_point(size = 1.5, 
            colour = "grey90") + ggplot2::geom_point(data = select, size = 3, color = color) + ggplot2::annotate("text", 
            x = 1.1, y = seq(y_max, y_min, length.out = topn)[1:length(top_genes)], label = top_genes, 
            hjust = 0, col = color, size = 4.5) + ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") + 
            ggplot2::theme_bw() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(), 
            axis.text = ggplot2::element_text(colour = "black", size = 15), axis.text.x = ggplot2::element_blank(), 
            axis.ticks.x = ggplot2::element_blank(), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), plot.margin = grid::unit(c(1, 5, 1, 1), "lines"), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 0.5)) + ggplot2::labs(x = "Genes", y = "log2 fold change", 
            title = i)
    })
    height = round(length(all_plots)/5 + 0.35, 0) * 5
    all_plots <- cowplot::plot_grid(plotlist = all_plots, ncol = 5)
    if (is.null(output)) 
        return(all_plots)
    ggplot2::ggsave(output, all_plots, width = 20, height = height, limitsize = FALSE)
}

# enrichr_cluster_plot：已有富集结果热图；不执行富集检验 ----
enrichr_cluster_plot <- function(data, select = "pvalue", topN = 5, color = "col2", outlier.rm = T) {
    data$plot_value <- -log10(pmax(data[, select], .Machine$double.xmin))
    data <- do.call(rbind, lapply(unique(data$type), function(i) {
        subdata <- data[data$type == i, ]
        head(subdata[order(subdata$plot_value, decreasing = TRUE), ], topN)
    }))
    order <- unique(data$id)
    data <- reshape2::acast(data, id ~ type, value.var = "plot_value", fun.aggregate = mean, na.rm = TRUE)
    data[data == "NaN"] <- 0
    data <- data[order, ]
    if (color == "col1") {
        col <- colorRampPalette(c("#1B0084", "#1B099A", "#1C14AE", "#1C20C0", "#1F2DCE", "#233DD6", "#2950CE", 
            "#2E68AB", "#377989", "#3F876A", "#46954D", "#53A036", "#65AB26", "#7BB41A", "#95BE16", "#AFC61C", 
            "#C6CE26", "#DAD636", "#EAE04E", "#F1EC74", "#F8F7AF", "#FFFFFF"))(50)
    }
    else if (color == "col2") {
        col <- rev(colorRampPalette(c("#360E24", "#5D1A58", "#723692", "#795CC3", "#7487E0", "#7EB0E4", 
            "#AAD2E3", "#E6F1F1"))(50))
    }
    else if (color == "col3") {
        col <- colorRampPalette(c("#352A87", "#343EB1", "#1558D9", "#036CE0", "#107AD9", "#1387D3", "#0997D1", 
            "#06A4C9", "#0CADBB", "#23B4A9", "#43BB97", "#6ABE83", "#8FBF73", "#AEBD66", "#CBBB5B", "#E6B94F", 
            "#FDBE3C", "#FBCF2C", "#F5E21E", "#F9FB0E"))(50)
    }
    if (outlier.rm == T) {
        all_value <- reshape2::melt(data)
        thres <- quantile(all_value$value, probs = c(0.05, 0.95))
        col <- circlize::colorRamp2(seq(from = thres[1], to = thres[2], length.out = length(col)), col)
    }
    ComplexHeatmap::Heatmap(data, col = col, show_row_names = T, show_column_names = T, column_names_rot = 45, 
        cluster_rows = F, show_row_dend = F, cluster_columns = F, show_column_dend = F, rect_gp = grid::gpar(col = "white", 
            lwd = 0.75))
}

# prop_density：按组展示嵌入空间二维密度 ----
prop_density <- function(datafilt, group = "group", coord = "umap") {
    info <- datafilt@meta.data
    input <- data.frame(Seurat::Embeddings(datafilt, coord))
    input$group <- info[, group]
    ggplot2::ggplot(input, ggplot2::aes(x = input[, 1], y = input[, 2])) + ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), 
        geom = "raster", h = 2, contour = FALSE) + ggplot2::geom_density2d(linewidth = 0.1, colour = "#FDAF9199", 
        alpha = 0.3, bins = 15, h = 3.5) + viridis::scale_fill_viridis(option = "magma") + ggplot2::theme_bw() + 
        ggplot2::labs(x = "UMAP1", y = "UMAP2", title = "") + ggplot2::theme(strip.background = ggplot2::element_rect(colour = NA, 
        fill = "grey90"), strip.text.x = ggplot2::element_text(size = 15), axis.title = ggplot2::element_text(colour = "black", 
        size = 15), legend.key.size = grid::unit(0.75, "cm"), legend.text = ggplot2::element_text(size = 12), 
        legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), strip.placement = "outside", 
        axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank()) + ggplot2::scale_x_continuous(expand = c(0, 
        0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::facet_wrap(~group, nrow = 1)
}

# prop_density2：二维密度的另一种配色 ----
prop_density2 <- function(datafilt, group = "group", coord = "umap", color = "col1") {
    if (color == "col1") {
        col <- colorRampPalette(c("#FAF3F7", "#E3BED1", "#CC89AC", "#B55487", "#9E2062"))(100)
    }
    else if (color == "col2") {
        col <- colorRampPalette(c("#FAF3F7", "#FAF3F7", "#CEC6DD", "#A399C3", "#786CAA", "#4D4091"))(100)
    }
    else if (color == "col3") {
        col <- colorRampPalette(c("#FAF3F7", "#FAF3F7", "#D8C5DD", "#B798C3", "#966BA9", "#753E90"))(100)
    }
    info <- datafilt@meta.data
    input <- data.frame(Seurat::Embeddings(datafilt, coord))
    input$group <- info[, group]
    ggplot2::ggplot(input, ggplot2::aes(x = input[, 1], y = input[, 2])) + ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density)), 
        geom = "raster", h = 3, contour = FALSE) + ggplot2::geom_density2d(linewidth = 0.1, colour = "grey50", 
        alpha = 0.3, bins = 15, h = 3) + ggplot2::scale_fill_gradientn(colors = col) + ggplot2::theme_bw() + 
        ggplot2::labs(x = "UMAP1", y = "UMAP2", title = "") + ggplot2::theme(strip.background = ggplot2::element_rect(colour = NA, 
        fill = "grey90"), strip.text.x = ggplot2::element_text(size = 15), axis.title = ggplot2::element_text(colour = "black", 
        size = 15), legend.key.size = grid::unit(0.75, "cm"), legend.text = ggplot2::element_text(size = 12), 
        legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), strip.placement = "outside", 
        axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank()) + ggplot2::scale_x_continuous(expand = c(0, 
        0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::facet_wrap(~group, nrow = 1)
}

# prop_back2back：两组细胞组成的背靠背比例图 ----
prop_back2back <- function(datafilt, group = "group", cluster = "celltype", order = TRUE) {
    info <- datafilt@meta.data
    info <- data.frame(group = info[, group], cluster = info[, cluster])
    input <- data.frame(table(info$group, info$cluster))
    names(input) <- c("group", "cluster", "prop")
    input <- do.call(rbind, lapply(unique(input$group), function(i) {
        subinput <- input[input$group == i, ]
        subinput$prop <- subinput$prop/sum(subinput$prop)
        subinput
    }))
    input$label <- paste0(round(input$prop * 100, 1), "%")
    if (order == TRUE) {
        name <- unique(gtools::mixedsort(input$cluster, decreasing = TRUE))
        input$cluster <- factor(input$cluster, levels = name)
    }
    ggplot2::ggplot(input, ggplot2::aes(x = cluster, fill = group)) + ggplot2::scale_fill_manual(values = c("#A231A1", 
        "#F2D7EE")) + ggplot2::geom_bar(stat = "identity", data = subset(input, group == unique(input$group)[1]), 
        ggplot2::aes(y = prop)) + ggplot2::geom_text(data = subset(input, group == unique(input$group)[1]), 
        ggplot2::aes(y = prop, label = label), size = 5, hjust = -0.1) + ggplot2::geom_bar(stat = "identity", 
        data = subset(input, group == unique(input$group)[2]), ggplot2::aes(y = prop * (-1))) + ggplot2::geom_text(data = subset(input, 
        group == unique(input$group)[2]), ggplot2::aes(y = prop * (-1), label = label), size = 5, hjust = 1.1) + 
        ggplot2::theme_bw() + ggplot2::theme(legend.key.size = grid::unit(1, "cm"), legend.text = ggplot2::element_text(size = 12), 
        legend.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black", 
            size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
            color = "black", linewidth = 1)) + ggplot2::coord_flip() + ggplot2::ylab("Proportion") + 
        ggplot2::xlab("") + ggplot2::ylim(-max(input$prop) - 0.035, max(input$prop) + 0.035)
}

# prop_back2back_lollipop：两组细胞比例 log2 比值图 ----
prop_back2back_lollipop <- function(datafilt, group = "group", group1 = "H", group2 = "L", cluster = "celltype_sig2") {
    info <- datafilt@meta.data
    info <- data.frame(group = info[, group], cluster = info[, cluster])
    input <- data.frame(table(info$group, info$cluster))
    names(input) <- c("group", "cluster", "prop")
    input <- do.call(rbind, lapply(unique(input$group), function(i) {
        subinput <- input[input$group == i, ]
        subinput$prop <- subinput$prop/sum(subinput$prop)
        subinput
    }))
    fc <- do.call(rbind, lapply(unique(input$cluster), function(i) {
        subinput <- input[input$cluster %in% i, ]
        data.frame(id = i, value = log2(subinput$prop[subinput$group == group1]/subinput$prop[subinput$group == 
            group2]))
    }))
    lollipop_plot <- function(input, order = TRUE, decreasing = T) {
        if (order == TRUE) {
            input <- input[order(input$id, decreasing = decreasing), ]
            input$id <- forcats::fct_inorder(input$id)
        }
        ggplot2::ggplot(input, ggplot2::aes(x = id, y = value)) + ggplot2::geom_segment(ggplot2::aes(x = id, 
            xend = id, y = 0, yend = value), color = "lightgray") + ggplot2::geom_point(color = "#5fc2c0", 
            size = 5) + ggplot2::scale_color_viridis_d() + ggplot2::geom_line(ggplot2::aes(id, value, 
            group = 1)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
            ggplot2::geom_hline(yintercept = 0, color = "grey20", linetype = "solid") + ggplot2::theme_bw() + 
            ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black", 
                size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1), legend.key.size = grid::unit(0.75, "cm"), legend.text = ggplot2::element_text(size = 15), 
                legend.title = ggplot2::element_blank(), legend.position = "right") + ggplot2::coord_flip()
    }
    lollipop_plot(input = fc, order = TRUE, decreasing = T)
}

# exp_distribution：表达分布散点及均值 ± SD ----
exp_distribution <- function(datafilt, assay = "RNA", select, group, order = TRUE, decreasing = TRUE, 
    show.point = TRUE, size.point = 2, method = "wilcox.test", ncol = 2, rotate = 45) {
    common_dotbox <- function(input, rotate = 45, order = TRUE, decreasing = F, show.point = TRUE) {
        aggdata1 <- aggregate(input$value, by = list(input$type), FUN = mean)
        aggdata2 <- aggregate(input$value, by = list(input$type), FUN = sd)
        new_data <- data.frame(aggdata1, sd = aggdata2[, 2])
        names(new_data)[1:2] <- c("id", "mean_value")
        if (order == TRUE) {
            levels <- new_data$id[order(new_data$mean_value, decreasing = decreasing)]
            input$type <- factor(input$type, levels = levels)
            new_data$id <- factor(new_data$id, levels = levels)
        }
        if (show.point == TRUE) {
            plot <- ggplot2::ggplot(new_data, ggplot2::aes(x = id, y = mean_value)) + ggrastr::geom_quasirandom_rast(data = input, 
                ggplot2::aes(x = type, y = value), colour = "grey90", size = 1, alpha = 1) + ggplot2::geom_pointrange(ggplot2::aes(ymin = mean_value - 
                sd, ymax = mean_value + sd))
        }
        else {
            plot <- ggplot2::ggplot(new_data, ggplot2::aes(x = id, y = mean_value)) + ggplot2::geom_pointrange(ggplot2::aes(ymin = mean_value - 
                sd, ymax = mean_value + sd))
        }
        plot + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none", 
            axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1)) + 
            ggplot2::labs(x = "", y = "Gene activity", title = "") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate, 
            hjust = 1))
    }
    if (length(select) == 1) {
        data <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[select, ]
        data <- data.frame(as.matrix(data))
        info <- datafilt@meta.data
        input <- data.frame(value = data[, 1], type = as.character(info[, group]))
        input <- na.omit(input)
        plot <- common_dotbox(input = input, rotate = rotate, order = order, decreasing = decreasing, 
            show.point = show.point)
        plot <- plot + ggplot2::ggtitle(select)
        if (length(unique(input$type)) == 2) {
            plot + ggpubr::stat_compare_means(data = input, ggplot2::aes(x = type, y = value), method = method, 
                size = 5, label.x = 1.1)
        }
        else {
            plot
        }
    }
    else {
        common <- intersect(select, rownames(datafilt))
        data <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[common, ]
        data <- data.frame(t(as.matrix(data)))
        allplots <- lapply(colnames(data), function(i) {
            info <- datafilt@meta.data
            input <- data.frame(value = data[, i], type = as.character(info[, group]))
            input <- na.omit(input)
            plot <- common_dotbox(input = input, rotate = rotate, order = order, decreasing = decreasing, 
                show.point = show.point)
            plot <- plot + ggplot2::ggtitle(i)
            if (length(unique(input$type)) == 2) {
                plot + ggpubr::stat_compare_means(data = input, ggplot2::aes(x = type, y = value), method = method, 
                  size = 5, label.x = 1.1)
            }
            else {
                plot
            }
        })
        cowplot::plot_grid(plotlist = allplots, ncol = ncol)
    }
}

# exp_distribution_ecdf：基因表达的累积分布曲线 ----
exp_distribution_ecdf <- function(datafilt, assay = "RNA", select, group, method = "wilcox.test", ncol = 2) {
    genes <- intersect(select, rownames(datafilt[[assay]]))
    expression <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[genes, , drop = FALSE]
    plots <- lapply(genes, function(gene) {
        input <- data.frame(value = as.numeric(expression[gene, ]), type = datafilt[[]][colnames(expression), 
            group])
        ecdf_plot(input, title = gene, method = method)
    })
    cowplot::plot_grid(plotlist = plots, ncol = ncol)
}

# exp_distribution_boxplot：基因表达箱线图 ----
exp_distribution_boxplot <- function(datafilt, assay = "RNA", select, group, method = "wilcox.test", 
    ncol = 2) {
    common_boxplot <- function(input, output = NULL, method, width = 3.5, height = 5) {
        input <- input
        value1 = quantile(input$value, probs = c(0.005), na.rm = T)
        value2 = quantile(input$value, probs = c(0.98), na.rm = T)
        col <- c("#0070B2", "#5CB3DA", "#B8E3EA", "#DA1735", "#F15E4C", "#FF9F99", "#A231A1", "#A37CB7", 
            "#F2D7EE", "#B91372", "#E93B8C", "#ECB2C8", "#FF7149", "#F7AE24", "#FBDD7E", "#679436", "#8BBE53", 
            "#CDE391", "#067D69", "#00A385", "#98D4C6", "#114B5F", "#028090", "#B2DBBF", "#A23E48", "#CD6981", 
            "#FBD0C0", "#788585", "#9CAEA9", "#CCDAD1")
        plot <- ggplot2::ggplot(input, ggplot2::aes(x = type, y = value, fill = type)) + ggplot2::stat_boxplot(geom = "errorbar", 
            width = 0.1, color = "black") + ggplot2::geom_boxplot(alpha = 1, outlier.shape = NA, size = 0.7, 
            width = 0.7) + ggplot2::scale_fill_manual(values = c("#FBD324", "#A37CB7", sample(col))) + 
            ggpubr::stat_compare_means(method = method, size = 5, label.x = 1.1, label.y = 0.99 * value2) + 
            ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none", 
            axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1)) + 
            ggplot2::labs(x = "", y = "Gene activity", title = "")
        if (is.null(output)) {
            plot
        }
        else {
            ggplot2::ggsave(output, plot, width = width, height = height)
        }
    }
    if (length(select) == 1) {
        data <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[select, ]
        data <- data.frame(as.matrix(data))
        info <- datafilt@meta.data
        input <- data.frame(value = data[, 1], type = as.character(info[, group]))
        input <- na.omit(input)
        plot <- common_boxplot(input = input, output = NULL, method = method, width = 3.5, height = 5)
        plot <- plot + ggplot2::ggtitle(select)
        plot
    }
    else {
        common <- intersect(select, rownames(datafilt))
        data <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[common, ]
        data <- data.frame(t(as.matrix(data)))
        allplots <- lapply(colnames(data), function(i) {
            info <- datafilt@meta.data
            input <- data.frame(value = data[, i], type = as.character(info[, group]))
            input <- na.omit(input)
            plot <- common_boxplot(input = input, output = NULL, method = method, width = 3.5, height = 5)
            plot <- plot + ggplot2::ggtitle(i)
            plot
        })
        cowplot::plot_grid(plotlist = allplots, ncol = ncol)
    }
}

# ecdf_plot：通用累积分布图；input 包含 value、type ----
ecdf_plot <- function(input, ylim = c(0, 1), title = NULL, method = "wilcox.test") {
    plot <- ggplot2::ggplot(input, ggplot2::aes(value, color = type)) + ggplot2::stat_ecdf(geom = "step", 
        linewidth = 1) + ggplot2::theme_bw() + ggplot2::labs(x = "Value", y = "Cumulative fraction", 
        color = NULL, title = title) + ggplot2::coord_cartesian(ylim = ylim)
    if (length(unique(input$type)) == 2) {
        test <- ggpubr::compare_means(value ~ type, data = input, method = method)
        plot <- plot + ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = paste0("P = ", 
            signif(test$p[1], 3)))
    }
    plot
}

# heatmap_diff_seurat：各群 top marker 的平均表达 Z-score 热图 ----
heatmap_diff_seurat <- function(datafilt, diffdata, group.by = "seurat_clusters", color = "blue2red", 
    topn = 10, col_rot = 0, label = TRUE, font.size = 9, outlier.rm = T) {
    if ("avg_log2FC" %in% names(diffdata)) 
        diffdata$logfc <- diffdata$avg_log2FC
    if ("p_val_adj" %in% names(diffdata)) 
        diffdata$adjp <- diffdata$p_val_adj
    if (color == "parula") {
        col <- colorRampPalette(c("#352A87", "#343EB1", "#1558D9", "#036CE0", "#107AD9", "#1387D3", "#0997D1", 
            "#06A4C9", "#0CADBB", "#23B4A9", "#43BB97", "#6ABE83", "#8FBF73", "#AEBD66", "#CBBB5B", "#E6B94F", 
            "#FDBE3C", "#FBCF2C", "#F5E21E", "#F9FB0E"))(100)
    }
    else if (color == "kovesi") {
        col <- colorRampPalette(c("#1B0084", "#1B099A", "#1C14AE", "#1C20C0", "#1F2DCE", "#233DD6", "#2950CE", 
            "#2E68AB", "#377989", "#3F876A", "#46954D", "#53A036", "#65AB26", "#7BB41A", "#95BE16", "#AFC61C", 
            "#C6CE26", "#DAD636", "#EAE04E", "#F1EC74", "#F8F7AF", "#FFFFFF"))(100)
    }
    else if (color == "ocean.haline") {
        col <- colorRampPalette(c("#2A186C", "#2E1B87", "#2B229F", "#1D35A1", "#11459A", "#0C5294", "#115D8F", 
            "#19678C", "#22718A", "#2B7A89", "#328388", "#388D88", "#3E9886", "#44A284", "#4CAC80", "#55B67A", 
            "#62BF72", "#73C869", "#8AD060", "#A4D75C", "#BEDC62", "#D5E171", "#EAE784", "#FDEF9A"))(100)
    }
    else if (color == "ocean.dense") {
        col <- colorRampPalette(c("#360E24", "#5D1A58", "#723692", "#795CC3", "#7487E0", "#7EB0E4", "#AAD2E3", 
            "#E6F1F1"))(100)
    }
    else if (color == "magma") {
        col <- c(colorRampPalette(c(viridis::magma(323, begin = 0.15)[1]))(10), viridis::magma(323, begin = 0.18))
    }
    else if (color == "white2blue") {
        col <- colorRampPalette(c("#FBF4F8", "#E5E0ED", "#BFC6DD", "#8CADCC", "#4E92BA", "#1871A8", "#085889", 
            "#003758"))(100)
    }
    else if (color == "blue2red") {
        col <- colorRampPalette(c("#0070b2", "#009bc7", "#5ec7dd", "#b8e3ea", "#f3f3f1", "#fccdb9", "#f79676", 
            "#f15e4c", "#da1735"))(100)
    }
    else if (color == "rainbow") {
        col <- colorRampPalette(c("#343391", "#0064af", "#0090cc", "#00b6db", "#01b7c2", "#53c0a3", "#8dcb8a", 
            "#bbd967", "#fbd324", "#f6bd25", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f", "#c02e2f", "#8b2a21"))(100)
    }
    diffdata <- data.frame(id = diffdata$gene, cluster = diffdata$cluster, logfc = diffdata$logfc, pvalue = diffdata$adjp)
    diffdata <- diffdata[diffdata$pvalue < 0.01 & diffdata$logfc > 0.5, ]
    select <- gtools::mixedsort(as.character(unique(diffdata$cluster)))
    diffdata <- do.call(rbind, lapply(select, function(i) {
        subdata <- diffdata[diffdata$cluster == i, ]
        subdata <- dplyr::top_n(subdata, topn, logfc)
        subdata[order(subdata$pvalue), ]
    }))
    expdata <- as.matrix(Seurat::GetAssayData(datafilt[unique(diffdata$id), ], slot = "data", assay = "RNA"))
    expdata <- data.frame(expdata, check.names = F)
    info <- datafilt@meta.data
    clustinfo <- data.frame(id = rownames(info), cluster = info[, group.by])
    input <- lapply(as.character(unique(clustinfo$cluster)), function(i) {
        select <- clustinfo$id[clustinfo$cluster == i]
        data.frame(rowMeans(expdata[, select, drop = FALSE]), check.names = F)
    })
    input <- do.call(cbind, input)
    colnames(input) <- as.character(unique(clustinfo$cluster))
    input <- input[, gtools::mixedsort(colnames(input))]
    input <- t(scale(t(input)))
    if (outlier.rm == T) {
        all_value <- reshape2::melt(input)
        thres <- quantile(all_value$value, probs = c(0.01, 0.99))
        col <- circlize::colorRamp2(seq(from = thres[1], to = thres[2], length.out = length(col)), col)
    }
    ComplexHeatmap::Heatmap(input, col = col, show_row_names = label, show_column_names = T, column_names_rot = col_rot, 
        cluster_rows = F, show_row_dend = F, cluster_columns = F, show_column_dend = F, row_names_gp = grid::gpar(fontsize = font.size))
}

# prop_plot_hca：堆积比例图；input 包含 Var1、Var2、Freq ----
prop_plot_hca <- function(input, rotate = 45, decreasing = T, species = "human") {
    input$Var2 <- as.character(input$Var2)
    name <- unique(gtools::mixedsort(input$Var2, decreasing = decreasing))
    input$Var2 <- factor(input$Var2, levels = name)
    if (species == "human") {
        col <- c(CD8T = "#B8E3EA", CD4T = "#0070B2", Treg = "#A6D719", NKcell = "#F15E4C", Bcell = "#A231A1", 
            DC = "#D5B81B", Macro = "#028090", Mono = "#B2DBBF", MAST = "#B32226", Neutrophil = "#00A9A3", 
            Fibroblast = "#FBD324", Endothelial = "#CD6981", tumor = "#9CAEA9", None = "#F6F5BD")
    }
    else if (species == "mouse") {
        col <- c(CD8T = "#B8E3EA", CD4T = "#0070B2", Treg = "#A6D719", NKcell = "#F15E4C", Bcell = "#A231A1", 
            DC = "#D5B81B", MoMac = "#B2DBBF", MAST = "#B32226", Neutrophil = "#00A9A3", Fibroblast = "#FBD324", 
            Endothelial = "#CD6981", tumor = "#9CAEA9", None = "#F6F5BD")
    }
    ggplot2::ggplot(input, ggplot2::aes(x = Var1, y = Freq, fill = Var2)) + ggplot2::geom_bar(stat = "identity", 
        position = "fill", width = 0.75) + ggplot2::scale_fill_manual(values = col) + ggplot2::theme_bw() + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.key.size = grid::unit(0.75, "cm"), 
            legend.text = ggplot2::element_text(size = 15), legend.title = ggplot2::element_blank(), 
            axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                size = 15)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::labs(x = "", y = "Proportion", 
        title = "") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate, hjust = 1))
}

# prop_sample_compare：先按样本计算比例，再比较两个实验组 ----
prop_sample_compare <- function(datafilt, sample = "sample", group = "group", cluster = "celltype", output.name = NULL, 
    ncol = 5) {
    metadata <- datafilt[[]]
    counts <- table(metadata[[sample]], metadata[[cluster]])
    input <- as.data.frame(prop.table(counts, margin = 1))
    names(input) <- c("sample", "celltype", "proportion")
    design <- unique(data.frame(sample = metadata[[sample]], group = metadata[[group]]))
    input <- merge(input, design, by = "sample")
    plots <- lapply(gtools::mixedsort(unique(as.character(input$celltype))), function(celltype) {
        values <- input[input$celltype == celltype, ]
        ggplot2::ggplot(values, ggplot2::aes(group, proportion, fill = group)) + ggplot2::geom_boxplot(width = 0.55, 
            outlier.shape = NA) + ggbeeswarm::geom_quasirandom(width = 0.2, size = 2) + ggpubr::stat_compare_means(method = "t.test") + 
            ggplot2::scale_fill_manual(values = c("#FBD324", "#A37CB7")) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") + 
            ggplot2::labs(x = NULL, y = "Cell proportion", title = celltype)
    })
    plot <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
    if (!is.null(output.name)) 
        ggplot2::ggsave(output.name, plot, width = 4 * ncol, height = 5 * ceiling(length(plots)/ncol), 
            limitsize = FALSE)
    plot
}

# heatmap_text：数值矩阵热图，标注超过绝对值阈值的数值 ----
heatmap_text <- function(input, color = "blue2red", order_name = T, cutoff = 0.2, col_rot = 45, cluster_row = T, 
    cluster_col = T, outlier.rm = T) {
    input <- as.matrix(input)
    if (color == "parula") {
        col <- colorRampPalette(c("#352A87", "#343EB1", "#1558D9", "#036CE0", "#107AD9", "#1387D3", "#0997D1", 
            "#06A4C9", "#0CADBB", "#23B4A9", "#43BB97", "#6ABE83", "#8FBF73", "#AEBD66", "#CBBB5B", "#E6B94F", 
            "#FDBE3C", "#FBCF2C", "#F5E21E", "#F9FB0E"))(100)
    }
    else if (color == "kovesi") {
        col <- colorRampPalette(c("#1B0084", "#1B099A", "#1C14AE", "#1C20C0", "#1F2DCE", "#233DD6", "#2950CE", 
            "#2E68AB", "#377989", "#3F876A", "#46954D", "#53A036", "#65AB26", "#7BB41A", "#95BE16", "#AFC61C", 
            "#C6CE26", "#DAD636", "#EAE04E", "#F1EC74", "#F8F7AF", "#FFFFFF"))(100)
    }
    else if (color == "ocean.haline") {
        col <- colorRampPalette(c("#2A186C", "#2E1B87", "#2B229F", "#1D35A1", "#11459A", "#0C5294", "#115D8F", 
            "#19678C", "#22718A", "#2B7A89", "#328388", "#388D88", "#3E9886", "#44A284", "#4CAC80", "#55B67A", 
            "#62BF72", "#73C869", "#8AD060", "#A4D75C", "#BEDC62", "#D5E171", "#EAE784", "#FDEF9A"))(100)
    }
    else if (color == "ocean.dense") {
        col <- colorRampPalette(c("#360E24", "#5D1A58", "#723692", "#795CC3", "#7487E0", "#7EB0E4", "#AAD2E3", 
            "#E6F1F1"))(100)
    }
    else if (color == "magma") {
        col <- c(colorRampPalette(c(viridis::magma(323, begin = 0.15)[1]))(10), viridis::magma(323, begin = 0.18))
    }
    else if (color == "white2blue") {
        col <- colorRampPalette(c("#892C29", "#E2431D", "#FACE1D", "#F0F7F0", "#01ADD8", "#265EAB", "#3C2D96"))(100)
    }
    else if (color == "blue2red") {
        col <- colorRampPalette(rev(c("#892C29", "#E2431D", "#FACE1D", "#F0F7F0", "#01ADD8", "#265EAB", 
            "#3C2D96")))(100)
    }
    else if (color == "rainbow") {
        col <- colorRampPalette(c("#343391", "#0064af", "#0090cc", "#00b6db", "#01b7c2", "#53c0a3", "#8dcb8a", 
            "#bbd967", "#fbd324", "#f6bd25", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f", "#c02e2f", "#8b2a21"))(100)
    }
    if (outlier.rm == T) {
        all_value <- reshape2::melt(input)
        thres <- quantile(all_value$value, probs = c(0.01, 0.99))
        col <- circlize::colorRamp2(seq(from = thres[1], to = thres[2], length.out = length(col)), col)
    }
    if (order_name == T) {
        input <- input[gtools::mixedsort(rownames(input)), ]
        input <- input[, gtools::mixedsort(colnames(input))]
    }
    else {
        input <- input
    }
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(input, col = col, name = "Value", show_row_names = T, 
        show_column_names = T, column_names_rot = col_rot, cluster_rows = cluster_row, show_row_dend = F, 
        cluster_columns = cluster_col, show_column_dend = F, rect_gp = grid::gpar(col = "white", lwd = 0.75), 
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (abs(input[i, j]) > cutoff) 
                grid::grid.text(round(input, 2)[i, j], x, y, gp = grid::gpar(fontsize = 10))
        }, heatmap_legend_param = list(legend_direction = "horizontal", legend_width = grid::unit(6, 
            "cm"))), heatmap_legend_side = "bottom")
}

# seurat_circ_prop：各组的细胞组成环形图 ----
seurat_circ_prop <- function(datafilt, group = "treatment", cluster = "celltype", xlab_order = NULL, 
    ncol = 2) {
    circ_prop <- function(data) {
        data <- data.frame(data)
        names(data)[1:2] <- c("type", "n")
        data$fraction = data$n/sum(data$n)
        data$ymax = cumsum(data$fraction)
        data$ymin = c(0, head(data$ymax, n = -1))
        data = data[order(data$n, decreasing = TRUE), ]
        myLabel = as.vector(data$type)
        myLabel = paste(myLabel, "(", round(data$n/sum(data$n) * 100, 2), "%)", sep = "")
        ggplot2::ggplot(data = data, ggplot2::aes(fill = type, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) + 
            ggplot2::geom_rect(colour = "grey30", show.legend = TRUE) + ggplot2::scale_fill_manual(values = c(col)) + 
            ggplot2::coord_polar(theta = "y") + ggplot2::labs(x = "", y = "", title = "") + ggplot2::xlim(2.5, 
            NA) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_blank()) + 
            ggplot2::theme(axis.ticks = ggplot2::element_blank()) + ggplot2::theme(panel.border = ggplot2::element_blank())
    }
    info <- datafilt@meta.data
    info <- data.frame(group = info[, group], cluster = info[, cluster])
    input <- data.frame(table(info$group, info$cluster))
    names(input) <- c("group", "cluster", "prop")
    number <- length(unique(input$cluster))
    if (number < 16) {
        col <- c("#fbd324", "#bbd967", "#8dcb8a", "#53c0a3", "#01b7c2", "#31A6C6", "#0090cc", "#0070b2", 
            "#784B8F", "#9D478D", "#716299", "#A25F60", "#C5456B", "#D6654E", "#E08D4E", "#EABB63")
    }
    else {
        col <- c("#fbd324", "#bbd967", "#8dcb8a", "#53c0a3", "#01b7c2", "#31A6C6", "#0090cc", "#0070b2", 
            "#784B8F", "#9D478D", "#716299", "#A25F60", "#C5456B", "#D6654E", "#E08D4E", "#EABB63")
        col <- colorRampPalette(col)(number)
    }
    if (is.character(xlab_order) == TRUE) {
        select = xlab_order
    }
    else {
        select = unique(input$group)
    }
    all_plots <- lapply(select, function(i) {
        subinput <- input[input$group == i, ]
        subinput <- subinput[, -1]
        plot <- circ_prop(subinput) + ggplot2::ggtitle(i)
        plot
    })
    cowplot::plot_grid(plotlist = all_plots, ncol = ncol)
}

# sc_pointplot_stat：均值 ± SD/SEM 与组间检验；input 包含 value、type ----
sc_pointplot_stat <- function(input, thres.p = 0.05, display = "sd", raster = FALSE, rotate = 45, xlab_order = TRUE, 
    output = NULL, method = "t.test", y.position = 0.9, width = 7, height = 5) {
    if (is.character(xlab_order) == TRUE) {
        input$type <- factor(input$type, levels = xlab_order)
    }
    if (length(unique(input$type)) < 17) {
        col <- c("grey75", "#c02e2f", "#0064af", "#f6bd25", "#bbd967", "#A23E48", "#53c0a3", "#B8E3EA", 
            "#00b6db", "#0090cc", "#01b7c2", "#8dcb8a", "#fbd324", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f")
    }
    else {
        col <- c("grey75", "#c02e2f", "#0064af", "#f6bd25", "#bbd967", "#A23E48", "#53c0a3", "#B8E3EA", 
            "#00b6db", "#0090cc", "#01b7c2", "#8dcb8a", "#fbd324", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f")
        col <- colorRampPalette(col)(length(unique(input$type)))
    }
    sem <- function(x) {
        return(sd(x)/sqrt(length(x)))
    }
    aggdata1 <- aggregate(input$value, by = list(input$type), FUN = mean)
    if (display == "sd") {
        aggdata2 <- aggregate(input$value, by = list(input$type), FUN = sd)
    }
    else if (display == "sem") {
        aggdata2 <- aggregate(input$value, by = list(input$type), FUN = sem)
    }
    new_data <- data.frame(aggdata1, sd = aggdata2[, 2])
    names(new_data)[1:2] <- c("id", "mean_value")
    if (is.character(xlab_order) == FALSE) {
        if (xlab_order == TRUE) {
            levels <- new_data$id[order(new_data$mean_value, decreasing = TRUE)]
            new_data$id <- factor(new_data$id, levels = levels)
            input$type <- factor(input$type, levels = levels)
        }
    }
    if (raster == FALSE) {
        plot <- ggplot2::ggplot(new_data, ggplot2::aes(x = id, y = mean_value)) + ggbeeswarm::geom_quasirandom(data = input, 
            ggplot2::aes(x = type, y = value, color = type), width = 0.2, size = 3, alpha = 0.1)
    }
    else {
        plot <- ggplot2::ggplot(new_data, ggplot2::aes(x = id, y = mean_value)) + ggrastr::geom_quasirandom_rast(data = input, 
            ggplot2::aes(x = type, y = value, color = type), width = 0.2, size = 3, alpha = 0.1)
    }
    plot <- plot + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_value - sd, ymax = mean_value + sd), 
        width = 0.1) + ggplot2::geom_crossbar(ggplot2::aes(ymin = mean_value, ymax = mean_value), width = 0.35, 
        linewidth = 0.35) + ggplot2::scale_color_manual(values = c(col)) + ggplot2::theme_classic() + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none", axis.text = ggplot2::element_text(colour = "black", 
            size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15)) + ggplot2::labs(x = "", 
        y = "Score", title = "") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate, 
        hjust = 1))
    if (length(unique(input$type)) == 2) {
        plot <- plot + ggpubr::stat_compare_means(data = input, ggplot2::aes(x = type, y = value), method = method, 
            size = 5, label.x = 1.1)
    }
    else {
        stat <- ggpubr::compare_means(value ~ type, data = input, method = method, p.adjust.method = "BH")
        stat <- stat[stat$p.adj < thres.p, ]
        plot <- plot + ggpubr::stat_pvalue_manual(stat, y.position = quantile(input$value, probs = c(y.position)), 
            step.increase = 0.2, label.size = 5, label = "p.adj")
    }
    if (is.null(output)) {
        plot
    }
    else {
        ggplot2::ggsave(output, plot, width = width, height = height)
    }
}

# diff_signal_acrossclust：各亚群的两组平均表达/评分之差 ----
diff_signal_acrossclust <- function(datafilt, assay = "RNA", select, group = "group", group1 = "H", group2 = "L", 
    cluster = "celltype_sig2") {
    lollipop_plot <- function(input, order = TRUE, decreasing = T) {
        if (order == TRUE) {
            input <- input[order(input$id, decreasing = decreasing), ]
            input$id <- forcats::fct_inorder(input$id)
        }
        ggplot2::ggplot(input, ggplot2::aes(x = id, y = value)) + ggplot2::geom_segment(ggplot2::aes(x = id, 
            xend = id, y = 0, yend = value), color = "lightgray") + ggplot2::geom_point(color = "#5fc2c0", 
            size = 5) + ggplot2::scale_color_viridis_d() + ggplot2::geom_line(ggplot2::aes(id, value, 
            group = 1)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
            ggplot2::geom_hline(yintercept = 0, color = "grey20", linetype = "solid") + ggplot2::theme_bw() + 
            ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black", 
                size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1), legend.key.size = grid::unit(0.75, "cm"), legend.text = ggplot2::element_text(size = 15), 
                legend.title = ggplot2::element_blank(), legend.position = "right") + ggplot2::labs(x = "", 
            y = "Mean difference") + ggplot2::coord_flip()
    }
    if (is.character(select) == TRUE) {
        data <- Seurat::GetAssayData(datafilt, slot = "data", assay = assay)[select, ]
        data <- data.frame(as.matrix(data))
        data <- data.frame(id = rownames(data), value = data[, 1])
    }
    else if (is.list(select) == TRUE) {
        datafilt <- Seurat::AddModuleScore(datafilt, select, assay = assay, nbin = 30, ctrl = 100, name = names(select))
        name <- colnames(datafilt@meta.data)
        position <- length(name) - length(names(select)) + 1
        colnames(datafilt@meta.data)[position:length(name)] <- names(select)
        sig <- datafilt@meta.data
        data <- data.frame(id = rownames(sig), value = sig[, names(select)])
    }
    info <- datafilt@meta.data
    info <- data.frame(id = rownames(info), group = info[, group], cluster = info[, cluster])
    input <- merge(data, info, by = "id")
    diffvalue <- do.call(rbind, lapply(unique(input$cluster), function(i) {
        subinput <- input[input$cluster == i, ]
        data.frame(id = i, value = mean(subinput$value[subinput$group == group1]) - mean(subinput$value[subinput$group == 
            group2]))
    }))
    lollipop_plot(input = diffvalue, order = TRUE, decreasing = T)
}

# varhclust_for_hmap：轨迹热图辅助：按趋势聚类变量，返回 id、type ----
varhclust_for_hmap <- function(input, reference, k = 2, nbin = 50, smooth.matrix = T, smooth.method = "loess", 
    rm.singlclust = T) {
    reference <- reference[order(reference$value), ]
    input <- input[, reference$id]
    bin_anno <- data.frame(id = colnames(input), number = 1:ncol(input))
    bin_number <- round(length(bin_anno$id)/nbin, 0)
    bin_anno$bin <- dplyr::ntile(bin_anno$number, n = bin_number)
    input <- data.frame(bin = bin_anno$bin, t(input), check.names = F)
    input <- aggregate(. ~ bin, input, mean)
    input <- t(input[, -1])
    input <- t(scale(t(input)))
    if (smooth.matrix == T) {
        if (smooth.method == "loess") {
            smooth <- function(row) {
                x <- seq_along(row)
                loess_fit <- loess(row ~ x, span = 1)
                predict(loess_fit)
            }
        }
        else if (smooth.method == "savitzky-golay") {
            smooth <- function(row, window_size = 9, poly_order = 1) {
                signal::sgolayfilt(row, p = poly_order, n = window_size)
            }
        }
        input <- t(apply(input, 1, smooth))
    }
    ht <- ComplexHeatmap::Heatmap(input, name = "Z-score", use_raster = TRUE, raster_quality = 2, show_column_names = F, 
        column_names_rot = 0, cluster_rows = T, show_row_dend = T, cluster_columns = F, show_column_dend = F)
    ht <- ComplexHeatmap::draw(ht)
    row_dendrogram <- ComplexHeatmap::row_dend(ht)
    row_hclust <- as.hclust(row_dendrogram)
    row_clusters <- cutree(row_hclust, k = k)
    result <- data.frame(id = names(row_clusters), type = paste0("C", row_clusters))
    if (rm.singlclust == T) {
        number <- data.frame(table(result$type))
        id <- as.character(number$Var1[number$Freq %in% c("1")])
        result <- result[!(result$type %in% id), ]
        result
    }
    else {
        result
    }
}

# heatmap_correlation：按已有参考值排序的趋势热图；不推断拟时序 ----
heatmap_correlation <- function(input, reference, varcluster, select.gene = NULL, color = "rainbow", 
    nbin = 50, smooth.matrix = T, smooth.method = "loess", order.by.ref = T, rownames = F) {
    reference <- reference[order(reference$value), ]
    varcluster <- varcluster[order(varcluster$id), ]
    inputdata <- input[varcluster$id, reference$id, drop = FALSE]
    bin_anno <- data.frame(id = colnames(inputdata), number = 1:ncol(inputdata))
    bin_number <- round(length(bin_anno$id)/nbin, 0)
    bin_anno$bin <- dplyr::ntile(bin_anno$number, n = bin_number)
    inputdata <- data.frame(bin = bin_anno$bin, t(inputdata), check.names = F)
    inputdata <- aggregate(. ~ bin, inputdata, mean)
    inputdata <- t(inputdata[, -1])
    inputdata <- t(scale(t(inputdata)))
    anno <- varcluster
    anno$type <- as.character(anno$type)
    ha_row <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(anno[, 2]), show_legend = FALSE, which = "row", 
        show_annotation_name = F)
    col <- rev(c("#360E24", "#5D1A58", "#723692", "#795CC3", "#7487E0", "#7EB0E4", "#AAD2E3", "#E6F1F1"))
    col = circlize::colorRamp2(seq(from = 0, to = 100, length.out = length(col)), col)
    ha_time = ComplexHeatmap::HeatmapAnnotation(time = seq(from = 0, to = 100, length.out = ncol(inputdata)), 
        col = list(time = col), show_legend = FALSE)
    common <- intersect(rownames(inputdata), select.gene)
    mark_at <- match(common, rownames(inputdata))
    ha_gene <- if (length(common)) 
        ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = mark_at, labels = common))
    else NULL
    if (smooth.matrix == T) {
        if (smooth.method == "loess") {
            smooth <- function(row) {
                x <- seq_along(row)
                loess_fit <- loess(row ~ x, span = 1)
                predict(loess_fit)
            }
        }
        else if (smooth.method == "savitzky-golay") {
            smooth <- function(row, window_size = 9, poly_order = 1) {
                signal::sgolayfilt(row, p = poly_order, n = window_size)
            }
        }
        inputdata <- t(apply(inputdata, 1, smooth))
    }
    if (color == "rainbow") {
        col2 <- c("#343391", "#0064af", "#0090cc", "#00b6db", "#01b7c2", "#53c0a3", "#8dcb8a", "#bbd967", 
            "#fbd324", "#f6bd25", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f", "#c02e2f", "#8b2a21")
    }
    else if (color == "yellow2green") {
        col2 <- rev(c("#FDE725", "#BAD741", "#77C75C", "#35B779", "#339C80", "#328287", "#31688E", "#37457A", 
            "#3D2367", "#440154"))
        col2 <- colorRampPalette(col2)(16)
    }
    else if (color == "ocean") {
        col2 <- c("#081D58", "#21338F", "#1F71B0", "#3BAFC2", "#92D4B8", "#E1F3B1", "#FFFFD9")
        col2 <- colorRampPalette(col2)(16)
    }
    all_value <- reshape2::melt(inputdata, id.vars = c("value"))
    thres <- quantile(all_value$value, probs = c(0.01, 0.99))
    col <- circlize::colorRamp2(seq(from = thres[1], to = thres[2], length.out = 16), col2)
    inputdata <- t(scale(t(inputdata)))
    if (order.by.ref == T) {
        get_peak_position <- function(seq) {
            x <- seq_along(seq)
            loess_fit <- loess(seq ~ x, span = 1)
            fitted_values <- predict(loess_fit)
            peak_position <- which.max(fitted_values)
            return(peak_position)
        }
        pseudo_bulk <- aggregate(inputdata, by = list(type = varcluster$type), FUN = mean)
        rownames(pseudo_bulk) <- pseudo_bulk$type
        pseudo_bulk <- as.matrix(pseudo_bulk[, -1])
        pseudo_peak_positions <- apply(pseudo_bulk, 1, get_peak_position)
        group_order <- order(pseudo_peak_positions)
        sorted_groups <- rownames(pseudo_bulk)[group_order]
        sorted_matrix <- NULL
        id_order <- c()
        for (group in sorted_groups) {
            group_ids <- varcluster$id[varcluster$type == group]
            group_matrix <- inputdata[group_ids, , drop = FALSE]
            dist_matrix <- dist(group_matrix, method = "euclidean")
            hc <- hclust(dist_matrix, method = "ward.D2")
            sorted_group_matrix <- group_matrix[hc$order, ]
            sorted_matrix <- rbind(sorted_matrix, sorted_group_matrix)
            id_order <- c(id_order, rownames(sorted_group_matrix))
        }
        inputdata <- sorted_matrix
        sorted_varcluster <- varcluster[match(id_order, varcluster$id), ]
        sorted_varcluster$type <- forcats::fct_inorder(sorted_varcluster$type)
        anno <- sorted_varcluster
        common <- intersect(rownames(inputdata), select.gene)
        mark_at <- match(common, rownames(inputdata))
        ha_gene <- if (length(common)) 
            ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = mark_at, labels = common))
        else NULL
        ha_row <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(anno[, 2]), show_legend = FALSE, 
            which = "row", show_annotation_name = F)
        if (rownames == T) {
            ComplexHeatmap::draw(ComplexHeatmap::Heatmap(inputdata, col = col, name = "Z-score", left_annotation = ha_row, 
                top_annotation = ha_time, use_raster = TRUE, raster_quality = 2, column_gap = grid::unit(0, 
                  "mm"), row_gap = grid::unit(0, "mm"), border = "#404040", show_row_names = rownames, 
                show_column_names = F, column_names_rot = 0, row_split = anno$type, cluster_rows = F, 
                show_row_dend = T, cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(legend_direction = "horizontal", 
                  legend_width = grid::unit(6, "cm"))), heatmap_legend_side = "bottom")
        }
        else {
            ComplexHeatmap::draw(ComplexHeatmap::Heatmap(inputdata, col = col, name = "Z-score", left_annotation = ha_row, 
                right_annotation = ha_gene, top_annotation = ha_time, use_raster = TRUE, raster_quality = 2, 
                column_gap = grid::unit(0, "mm"), row_gap = grid::unit(0, "mm"), border = "#404040", 
                show_row_names = rownames, row_split = anno$type, show_column_names = F, column_names_rot = 0, 
                cluster_rows = F, show_row_dend = T, cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(legend_direction = "horizontal", 
                  legend_width = grid::unit(6, "cm"))), heatmap_legend_side = "bottom")
        }
    }
    else {
        if (rownames == T) {
            ComplexHeatmap::draw(ComplexHeatmap::Heatmap(inputdata, col = col, name = "Z-score", left_annotation = ha_row, 
                top_annotation = ha_time, use_raster = TRUE, raster_quality = 2, column_gap = grid::unit(0, 
                  "mm"), row_gap = grid::unit(0, "mm"), border = "#404040", show_row_names = rownames, 
                show_column_names = F, column_names_rot = 0, row_split = anno$type, cluster_rows = T, 
                show_row_dend = T, cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(legend_direction = "horizontal", 
                  legend_width = grid::unit(6, "cm"))), heatmap_legend_side = "bottom")
        }
        else {
            ComplexHeatmap::draw(ComplexHeatmap::Heatmap(inputdata, col = col, name = "Z-score", left_annotation = ha_row, 
                right_annotation = ha_gene, top_annotation = ha_time, use_raster = TRUE, raster_quality = 2, 
                column_gap = grid::unit(0, "mm"), row_gap = grid::unit(0, "mm"), border = "#404040", 
                show_row_names = rownames, show_column_names = F, column_names_rot = 0, row_split = anno$type, 
                cluster_rows = T, show_row_dend = T, cluster_columns = F, show_column_dend = F, heatmap_legend_param = list(legend_direction = "horizontal", 
                  legend_width = grid::unit(6, "cm"))), heatmap_legend_side = "bottom")
        }
    }
}

# hmap_cluster_vis：变量模块沿参考值变化的平滑曲线 ----
hmap_cluster_vis <- function(input, reference, varcluster, smooth.SE = F) {
    input <- t(input)
    allplots <- lapply(unique(varcluster$type), function(i) {
        id <- varcluster$id[varcluster$type == i]
        subinput <- scale(input[, id, drop = FALSE])
        meanvalue <- data.frame(id = names(rowMeans(subinput)), meanvalue = unname(rowMeans(subinput)), 
            check.names = F)
        data <- merge(meanvalue, reference, by = "id")
        ggplot2::ggplot(data, ggplot2::aes(value, meanvalue)) + ggplot2::geom_smooth(span = 0.7, method = "loess", 
            se = smooth.SE, color = "#0070B2") + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
            axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                size = 15), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1)) + 
            ggplot2::labs(x = "Pseudotime", y = "Activity") + ggplot2::ggtitle(i)
    })
    cowplot::plot_grid(plotlist = allplots, ncol = 3)
}

# hmap_variable_vis：指定变量沿参考值变化的曲线，可分组 ----
hmap_variable_vis <- function(input, input.group = NULL, reference, select, smooth.SE = F, ncol = 3) {
    input <- t(input)
    col <- c("#c02e2f", "#0064af", "#f6bd25", "#bbd967", "#A23E48", "#53c0a3", "#B8E3EA", "#00b6db", 
        "#0090cc", "#01b7c2", "#8dcb8a", "#fbd324", "#f4a02e", "#ed6f32", "#ea5c2e", "#d5452f")
    if (is.null(input.group)) {
        if (length(select) == 1) {
            subinput <- data.frame(id = rownames(input), datavalue = input[, select], check.names = F)
            data <- merge(subinput, reference, by = "id")
            ggplot2::ggplot(data, ggplot2::aes(value, datavalue)) + ggplot2::geom_smooth(span = 0.7, 
                method = "loess", se = smooth.SE, color = "#0070B2") + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                  size = 15), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1), 
                legend.key.size = grid::unit(1, "cm"), legend.text = ggplot2::element_text(size = 12), 
                legend.title = ggplot2::element_blank()) + ggplot2::labs(x = "Pseudotime", y = "Activity") + 
                ggplot2::ggtitle(select)
        }
        else {
            allplots <- lapply(select, function(i) {
                subinput <- data.frame(id = rownames(input), datavalue = input[, i], check.names = F)
                data <- merge(subinput, reference, by = "id")
                ggplot2::ggplot(data, ggplot2::aes(value, datavalue)) + ggplot2::geom_smooth(span = 0.7, 
                  method = "loess", se = smooth.SE, color = "#0070B2") + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                  axis.text = ggplot2::element_text(colour = "black", size = 15), axis.title = ggplot2::element_text(colour = "black", 
                    size = 15), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1), 
                  legend.key.size = grid::unit(1, "cm"), legend.text = ggplot2::element_text(size = 12), 
                  legend.title = ggplot2::element_blank()) + ggplot2::labs(x = "Pseudotime", y = "Activity") + 
                  ggplot2::ggtitle(i)
            })
            cowplot::plot_grid(plotlist = allplots, ncol = ncol)
        }
    }
    else {
        if (length(select) == 1) {
            subinput <- data.frame(id = rownames(input), datavalue = input[, select], check.names = F)
            subinput <- merge(subinput, input.group, by = "id")
            data <- merge(subinput, reference, by = "id")
            ggplot2::ggplot(data, ggplot2::aes(value, datavalue, color = group)) + ggplot2::geom_smooth(span = 0.7, 
                method = "loess", se = smooth.SE) + ggplot2::scale_color_manual(values = c(col)) + ggplot2::theme_bw() + 
                ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black", 
                  size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                  color = "black", linewidth = 1), legend.key.size = grid::unit(1, "cm"), legend.text = ggplot2::element_text(size = 12), 
                  legend.title = ggplot2::element_blank()) + ggplot2::labs(x = "Pseudotime", y = "Activity") + 
                ggplot2::ggtitle(select)
        }
        else {
            allplots <- lapply(select, function(i) {
                subinput <- data.frame(id = rownames(input), datavalue = input[, i], check.names = F)
                subinput <- merge(subinput, input.group, by = "id")
                data <- merge(subinput, reference, by = "id")
                ggplot2::ggplot(data, ggplot2::aes(value, datavalue, color = group)) + ggplot2::geom_smooth(span = 0.7, 
                  method = "loess", se = smooth.SE) + ggplot2::scale_color_manual(values = c(col)) + 
                  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black", 
                  size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                  color = "black", linewidth = 1), legend.key.size = grid::unit(1, "cm"), legend.text = ggplot2::element_text(size = 12), 
                  legend.title = ggplot2::element_blank()) + ggplot2::labs(x = "Pseudotime", y = "Activity") + 
                  ggplot2::ggtitle(i)
            })
            cowplot::plot_grid(plotlist = allplots, ncol = ncol)
        }
    }
}

# heatmap_dotplot_marker：marker 平均表达热图及表达比例圆圈 ----
heatmap_dotplot_marker <- function(datafilt, group.by = "celltype_sig2", marker, color = "blue2red", 
    dot.scale = 1) {
    if (is.list(marker)) {
        feature <- stack(marker)
        names(feature) <- c("gene", "term")
    }
    else {
        feature <- read.csv(marker, stringsAsFactors = FALSE)
    }
    feature <- feature[!duplicated(feature$gene) & feature$gene %in% rownames(datafilt[["RNA"]]), ]
    dot <- Seurat::DotPlot(datafilt, assay = "RNA", group.by = group.by, features = feature$gene)$data
    mean_expression <- reshape2::acast(dot, features.plot ~ id, value.var = "avg.exp.scaled")
    proportion <- reshape2::acast(dot, features.plot ~ id, value.var = "pct.exp")
    mean_expression <- mean_expression[feature$gene, , drop = FALSE]
    proportion <- proportion[rownames(mean_expression), colnames(mean_expression), drop = FALSE]
    palettes <- list(blue2red = c("#0070B2", "#F3F3F1", "#DA1735"), parula = c("#352A87", "#06A4C9", 
        "#F9FB0E"), ocean.dense = c("#360E24", "#7487E0", "#E6F1F1"), ocean2 = c("#7F256F", "#32B1C0", 
        "#F4F2D9"), rainbow = c("#343391", "#53C0A3", "#FBD324", "#8B2A21"))
    colors <- (grDevices::colorRampPalette(palettes[[color]]))(100)
    annotation <- ComplexHeatmap::rowAnnotation(marker_set = feature$term, show_annotation_name = FALSE, 
        show_legend = TRUE)
    heatmap <- ComplexHeatmap::Heatmap(mean_expression, col = colors, name = "Mean expression (Z)", left_annotation = annotation, 
        row_split = factor(feature$term, levels = unique(feature$term)), cluster_rows = FALSE, cluster_columns = FALSE, 
        cluster_row_slices = FALSE, column_names_rot = 45, column_title = "Circle area: expressing-cell proportion", 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.circle(x, y, r = sqrt(proportion[i, j]/100) * 0.45 * dot.scale * min(grid::unit.c(width, 
                height)), gp = grid::gpar(fill = NA, col = "black"))
        })
    ComplexHeatmap::draw(heatmap)
}

# common_boxplot3：通用散点与均值线；input 包含 value、type ----
common_boxplot3 <- function(input, output = NULL, rotate = 45, order = TRUE, decreasing = T, width = 7, 
    height = 5) {
    input <- input
    col <- c("#0070B2", "#5CB3DA", "#B8E3EA", "#DA1735", "#F15E4C", "#FF9F99", "#A231A1", "#A37CB7", 
        "#F2D7EE", "#B91372", "#E93B8C", "#ECB2C8", "#FF7149", "#F7AE24", "#FBDD7E", "#679436", "#8BBE53", 
        "#CDE391", "#067D69", "#00A385", "#98D4C6", "#114B5F", "#028090", "#B2DBBF", "#A23E48", "#CD6981", 
        "#FBD0C0", "#788585", "#9CAEA9", "#CCDAD1")
    mean_data <- aggregate(input$value, by = list(input$type), FUN = mean)
    names(mean_data)[1:2] <- c("id", "mean_value")
    if (order == TRUE) {
        levels <- mean_data$id[order(mean_data$mean_value, decreasing = decreasing)]
        input$type <- factor(input$type, levels = levels)
        mean_data$id <- factor(mean_data$id, levels = levels)
    }
    if (length(unique(input$type)) < 36) {
        plot <- ggplot2::ggplot(mean_data, ggplot2::aes(x = id, y = mean_value)) + ggplot2::geom_crossbar(ggplot2::aes(ymin = mean_value, 
            ymax = mean_value), width = 0.5, linewidth = 0.5) + ggbeeswarm::geom_quasirandom(data = input, 
            ggplot2::aes(x = type, y = value, fill = type, color = type), width = 0, size = 3, alpha = 0.2) + 
            ggplot2::scale_color_manual(values = c(sample(col))) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
            legend.position = "none", axis.text = ggplot2::element_text(colour = "black", size = 15), 
            axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1)) + ggplot2::labs(x = "", y = "Gene activity", title = "") + 
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate, hjust = 1))
    }
    else {
        plot <- ggplot2::ggplot(mean_data, ggplot2::aes(x = id, y = mean_value)) + ggplot2::geom_crossbar(ggplot2::aes(ymin = mean_value, 
            ymax = mean_value), width = 0.5, linewidth = 0.5) + ggbeeswarm::geom_quasirandom(data = input, 
            ggplot2::aes(x = type, y = value), width = 0, size = 3, alpha = 0.2) + ggplot2::theme_bw() + 
            ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none", axis.text = ggplot2::element_text(colour = "black", 
                size = 15), axis.title = ggplot2::element_text(colour = "black", size = 15), panel.border = ggplot2::element_rect(fill = NA, 
                color = "black", linewidth = 1)) + ggplot2::labs(x = "", y = "Gene activity", title = "") + 
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate, hjust = 1))
    }
    if (is.null(output)) {
        plot
    }
    else {
        ggplot2::ggsave(output, plot, width = width, height = height)
    }
}

# seurat_cytometry：两个指标的密度图及四象限比例 ----
seurat_cytometry <- function(metadata, select_x = "sig1", select_y = "sig2", select_group = "group", 
    median.thres = TRUE, density.bin = 30, density.h = 0.7, xcutoff = 0.5, ycutoff = 0.75, color = "col7") {
    input <- data.frame(x = metadata[[select_x]], y = metadata[[select_y]], group = if (is.null(select_group)) 
        "All"
    else metadata[[select_group]])
    palettes <- list(col1 = c("#F7D0B5", "#CB1B4F", "#1B112B"), col7 = c("#E6F1F1", "#7487E0", "#360E24"), 
        col8 = c("#F4F2D9", "#32B1C0", "#7F256F"), col9 = c("#FFFFD9", "#3BAFC2", "#081D58"), col11 = c("#FBF4F8", 
            "#8CADCC", "#003758"))
    plots <- lapply(unique(input$group), function(group) {
        values <- input[input$group == group, ]
        x_threshold <- if (median.thres) 
            median(values$x)
        else xcutoff
        y_threshold <- if (median.thres) 
            median(values$y)
        else ycutoff
        fractions <- c(mean(values$x <= x_threshold & values$y > y_threshold), mean(values$x <= x_threshold & 
            values$y <= y_threshold), mean(values$x > x_threshold & values$y > y_threshold), mean(values$x > 
            x_threshold & values$y <= y_threshold))
        labels <- data.frame(x = c(-Inf, -Inf, Inf, Inf), y = c(Inf, -Inf, Inf, -Inf), label = paste0(round(100 * 
            fractions, 2), "%"), hjust = c(-0.1, -0.1, 1.1, 1.1), vjust = c(1.5, -0.5, 1.5, -0.5))
        ggplot2::ggplot(values, ggplot2::aes(x, y)) + ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)), 
            geom = "polygon", contour_var = "ndensity", bins = density.bin, h = density.h) + ggplot2::geom_density_2d(color = "black", 
            linewidth = 0.2, contour_var = "ndensity", bins = density.bin, h = density.h) + ggplot2::scale_fill_gradientn(colors = palettes[[color]]) + 
            ggplot2::geom_vline(xintercept = x_threshold, linetype = "dashed") + ggplot2::geom_hline(yintercept = y_threshold, 
            linetype = "dashed") + ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label, 
            hjust = hjust, vjust = vjust), inherit.aes = FALSE, size = 3.5) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + 
            ggplot2::labs(x = select_x, y = select_y, title = group, fill = "Relative density")
    })
    cowplot::plot_grid(plotlist = plots, ncol = min(2, length(plots)))
}

