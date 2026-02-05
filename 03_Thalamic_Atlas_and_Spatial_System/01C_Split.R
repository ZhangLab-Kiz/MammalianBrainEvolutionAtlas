library(Seurat)

group_by <- "orig.ident"

output_dir <- "/data/xiaotx/Th/Figure/02_FeaturePlot/Spilt_by_Species"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

positive_genes <- c(
    "OTX2", "CHRM2", "SERPINF1", "LMO1", "GRIK1", "LHX1",
    "GAD2", "GATA3", "DMBX1", "TMEM145", "ASIC4",
    "TMEM260", "POU2F1", "SYNPR"
)

valid_score_genes <- positive_genes[positive_genes %in% rownames(seu)]

# ==============================================================================
# Calculate Module Score
# ==============================================================================

seu <- AddModuleScore(
    object = seu,
    features = list(valid_score_genes),
    name = "OTX2_Network_Score"
)

score_name <- "OTX2_Network_Score1"

features_list <- c("OTX2", score_name)
colors_list <- c("blue", "red")

species_list <- c("All", unique(as.character(seu[[group_by]][, 1])))

# ==============================================================================
# Loop: Feature Plots by Species
# ==============================================================================

for (sp in species_list) {

    if (sp == "All") {
        sub_seu <- seu
    } else {
        cells_to_keep <- rownames(seu@meta.data)[seu@meta.data[[group_by]] == sp]
        if (length(cells_to_keep) == 0) {
            message(sprintf("Skip: %s (No cells)", sp))
            next
        }
        sub_seu <- subset(seu, cells = cells_to_keep)
    }

    for (feature in features_list) {
        is_score <- feature == score_name

        if (!is_score && !(feature %in% rownames(sub_seu))) {
            next
        }

        for (plot_color in colors_list) {

            color_suffix <- paste0("_", plot_color)
            high_color_param <- plot_color

            if (is_score) {
                file_tag <- "OTX2_Network_Module"

                plot_base <- FeaturePlot(sub_seu,
                    features = feature,
                    reduction = "ref.umap",
                    label = FALSE, raster = FALSE,
                    min.cutoff = 0, max.cutoff = 1
                ) +

                    scale_color_gradient(low = "lightgrey", high = high_color_param, limits = c(0, 1)) +
                    theme(aspect.ratio = 1, plot.title = element_blank())
            } else {
                file_tag <- paste0(feature, "_Expression")

                plot_base <- FeaturePlot(sub_seu,
                    features = feature,
                    reduction = "ref.umap",
                    order = TRUE,
                    label = FALSE, raster = FALSE,
                    max.cutoff = 2
                ) +
                    scale_color_gradient(low = "lightgrey", high = high_color_param) +
                    theme(aspect.ratio = 1, plot.title = element_blank())
            }

            p_pure <- plot_base + NoAxes() + NoLegend()
            p_for_legend <- plot_base + NoAxes()

            ggsave(
                filename = file.path(output_dir, paste0(file_tag, "_", sp, color_suffix, "_pure.png")),
                plot = p_pure,
                device = "png", width = 12, height = 12, units = "in", dpi = 600
            )

            g_obj <- ggplot2::ggplotGrob(p_for_legend)
            legend_index <- which(sapply(g_obj$grobs, function(x) x$name) == "guide-box")

            if (length(legend_index) > 0) {
                legend_content <- g_obj$grobs[[legend_index[1]]]

                legend_filename <- if (is_score) {
                    paste0(file_tag, "_", sp, color_suffix, "_Legend.pdf")
                } else {
                    paste0(file_tag, "_", sp, color_suffix, "_legend.pdf")
                }

                pdf(file.path(output_dir, legend_filename), width = 2, height = 4)
                grid::grid.draw(legend_content)
                dev.off()
            }
        }
    }

    if (sp != "All") {
        rm(sub_seu)
        gc()
    }
}

# ==============================================================================
# Loop: Split UMAP Hightlighting Target Cellgroup
# ==============================================================================

output_dir <- "/data/xiaotx/Th/Figure/01_Umap/Split_by_Species"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

target_celltype <- "Midbrain-derived inhibitory"
highlight_color <- "red"
bg_color <- "grey90"

species_list <- c("All", unique(as.character(seu$orig.ident)))

for (sp in species_list) {

    if (sp == "All") {
        sub_seu <- seu
    } else {
        if (sum(seu$orig.ident == sp) == 0) {
            next
        }
        sub_seu <- subset(seu, subset = orig.ident == sp)
    }

    sub_seu$highlight_group <- ifelse(
        sub_seu$final.celltype == target_celltype,
        "Target",
        "Background"
    )

    has_target <- "Target" %in% unique(sub_seu$highlight_group)

    if (has_target) {
        plot_order <- "Target"
    } else {
        plot_order <- NULL
    }

    my_cols <- c(
        "Target"     = highlight_color,
        "Background" = bg_color
    )

    p_pure <- DimPlot(sub_seu,
        reduction = "ref.umap",
        group.by = "highlight_group",
        cols = my_cols,
        order = plot_order,
        label = FALSE,
        raster = FALSE
    ) +
        NoAxes() +
        NoLegend() +
        theme(
            aspect.ratio = 1,
            plot.title = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        )

    file_name <- sprintf("Umap_Highlight_Midbrain_%s.png", sp)
    output_path <- file.path(output_dir, file_name)

    ggsave(
        filename = output_path,
        plot = p_pure,
        device = "png",
        width = 12, height = 12, units = "in", dpi = 600
    )

    if (sp != "All") {
        rm(sub_seu)
        gc()
    }
}
