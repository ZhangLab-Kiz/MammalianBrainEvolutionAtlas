library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

# ==============================================================================
# Load Mouse ST Data
# ==============================================================================

query <- readRDS("/data/xiaotx/ABC_Atlas/processed/wholebrain/st_whole.rds")

query_subset1 <- subset(query,
    subset = brain_section_label == "C57BL6J-638850.36"
)

query_subset2 <- subset(query,
    subset = brain_section_label == "C57BL6J-638850.44"
)

non_neuronal_classes <- c(
    "30 Astro-Epen",
    "31 OPC-Oligo",
    "32 OEC",
    "33 Vascular",
    "34 Immune"
)

output_dir <- "/data/xiaotx/Th/Figure/05_Mouse_ST_FeaturePlot/Neuron_Only"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

st_objects_list_raw <- list(
    `C57BL6J-638850.36` = query_subset1,
    `C57BL6J-638850.44` = query_subset2
)

positive_genes_mouse <- c(
    "Otx2", "Chrm2", "Serpinf1", "Lmo1", "Grik1", "Lhx1",
    "Gad2", "Gata3", "Dmbx1", "Tmem145", "Asic4",
    "Tmem260", "Pou2f1", "Synpr"
)

score_name_base <- "OTX2_Network_Score"
score_name_meta <- paste0(score_name_base, "1")

# ==============================================================================
# Process Slices: Compute Module Scores
# ==============================================================================

neuron_objects_list <- list()

for (slide_name in names(st_objects_list_raw)) {
    original_obj <- st_objects_list_raw[[slide_name]]

    neuron_obj <- subset(original_obj, subset = !class %in% non_neuronal_classes)

    spatial_coords <- neuron_obj@meta.data[, c("x", "y")]
    colnames(spatial_coords) <- c("spatial_1", "spatial_2")

    neuron_obj <- NormalizeData(neuron_obj)

    neuron_obj[["spatial"]] <- CreateDimReducObject(
        embeddings = as.matrix(spatial_coords),
        key = "spatial_",
        assay = DefaultAssay(neuron_obj)
    )

    valid_genes <- positive_genes_mouse[positive_genes_mouse %in% rownames(neuron_obj)]

    if (length(valid_genes) > 0) {
        neuron_obj <- AddModuleScore(
            object = neuron_obj,
            features = list(valid_genes),
            name = score_name_base,
            ctrl = 5
        )
    } else {
        neuron_obj[[score_name_meta]] <- 0
    }

    neuron_objects_list[[slide_name]] <- neuron_obj
}

# ==============================================================================
# Visualization Loop
# ==============================================================================

features_to_plot <- c("Otx2", score_name_meta)

for (feature in features_to_plot) {
    is_score <- feature == score_name_meta

    for (slide_name in names(neuron_objects_list)) {
        slide_obj <- neuron_objects_list[[slide_name]]

        if (!is_score && !(feature %in% rownames(slide_obj))) {
            message(paste("Feature", feature, "not found in", slide_name))
            next
        }

        cat(sprintf("Processing: %s - %s\n", slide_name, feature))

        if (is_score) {
            file_tag <- "OTX2_Network_Module"

            plot_layer <- scale_color_gradient2(
                low = "blue", mid = "white", high = "red",
                midpoint = 0,
                limits = c(-2, 3),
                oob = scales::squish
            )
        } else {
            file_tag <- paste0(feature, "_Expression")

            plot_layer <- scale_color_gradient(
                low = "lightgrey", high = "red",
                limits = c(0, 6),
                oob = scales::squish
            )
        }

        p_base <- FeaturePlot(
            object = slide_obj,
            features = feature,
            reduction = "spatial",
            pt.size = 0.01,
            order = TRUE,
            raster = FALSE
        ) +
            plot_layer +
            scale_y_reverse() +
            theme_classic() +
            theme(
                aspect.ratio = 0.6,
                plot.title = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank()
            )

        p_pure <- p_base + NoLegend() + NoAxes()

        output_filename <- paste0(file_tag, "_", slide_name, "_Neuron_pure.png")

        ggsave(
            filename = file.path(output_dir, output_filename),
            plot = p_pure,
            device = "png", width = 8, height = 4.8, units = "in", dpi = 600
        )

        p_for_legend <- p_base

        g_legend <- ggplot2::ggplotGrob(p_for_legend)
        legend_index <- which(sapply(g_legend$grobs, function(x) x$name) == "guide-box")

        if (length(legend_index) > 0) {
            legend_grob <- g_legend$grobs[[legend_index[1]]]

            legend_filename <- paste0(file_tag, "_", slide_name, "_Legend.pdf")

            pdf(file.path(output_dir, legend_filename), width = 2, height = 4)
            grid::grid.draw(legend_grob)
            dev.off()
        }
    }
}
