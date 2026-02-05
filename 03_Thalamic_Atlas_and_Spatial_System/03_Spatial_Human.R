setwd("/data/xiaotx/Th")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)

options(future.globals.maxSize = +Inf)

# ==============================================================================
# Load Data
# ==============================================================================

sc_data <- readRDS("/data/xiaotx/Th/Object/Th_final.rds")

ref_macaque <- subset(sc_data, subset = orig.ident == "Macaque")

rm(sc_data)
gc()

spatial_data <- readRDS("/data/xiaotx/Th/ST/spatial_list.rds")

sample_names <- c("GW10", "GW23_L1", "GW23_L2", "GW23_L3", "GW16", "GW23_M2", "GW23_vM", "GW21", "GW23_PV", "GW21_2")

target_sample_name <- "GW23_vM"

target_index <- which(sample_names == target_sample_name)

human_subset1 <- spatial_data[[target_index]]

target_sample_name <- "GW23_L1"

target_index <- which(sample_names == target_sample_name)

human_subset2 <- spatial_data[[target_index]]

# ==============================================================================
# Prepare Objects for Integration
# ==============================================================================

ref_macaque[["RNA"]]$counts <- ref_macaque[["RNA"]]$counts.Macaque
ref_macaque[["RNA"]]$counts.Macaque <- NULL
ref_macaque[["RNA"]]$data.Macaque <- NULL
ref_macaque[["RNA"]]$scale.data.Macaque <- NULL

human_subset1[["RNA"]] <- human_subset1[["Vizgen"]]
DefaultAssay(human_subset1) <- "RNA"

human_subset2[["RNA"]] <- human_subset2[["Vizgen"]]
DefaultAssay(human_subset2) <- "RNA"

obj.list <- list()

obj.list[["Human_ST_GW23vM"]] <- human_subset1
obj.list[["Human_ST_GW23L1"]] <- human_subset2

obj.list[["Macaque_Sc"]] <- ref_macaque

for (name in names(obj.list)) {
  obj.list[[name]]$dataset <- name
}

macaque_obj <- obj.list[["Macaque_Sc"]]
macaque_obj$dataset <- paste0(macaque_obj$dataset, "_batch", macaque_obj$batch)
obj.list[["Macaque_Sc"]] <- macaque_obj

common_genes <- Reduce(intersect, lapply(obj.list, rownames))
cat(sprintf("Found %d common genes for integration.\n", length(common_genes)))

obj.list <- lapply(obj.list, function(obj) {
  obj <- subset(obj, features = common_genes)
  return(obj)
})

gc()

obj.list <- lapply(obj.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
  return(x)
})

vfs <- SelectIntegrationFeatures(object.list = obj.list)

marker_file_path <- "/data/xiaotx/Th/Script/Cell_type_markers.R"
source(marker_file_path)

marker_lists <- list(
  mNeurons, mExN, mInN, mTEx,
  mEN1, mEN2, mMidbrainN,
  mTNR_N, mPretectumN, mTelencephalonN,
  mOPC, mOligo, mAstro, mMicro, mEpendy,
  mVas, mFibro
)
all_markers <- unique(unlist(marker_lists))

marker_genes <- intersect(common_genes, all_markers)
vfs <- unique(c(vfs, marker_genes))

if ("OTX2" %in% common_genes) {
  if (!("OTX2" %in% vfs)) {
    vfs <- union(vfs, "OTX2")
    cat("Info: OTX2 was not in the feature set. It has been added.\n")
  } else {
    cat("Info: OTX2 was already present in the integration features.\n")
  }
} else {
  cat("Warning: OTX2 is not present in both datasets.\n")
}

# ==============================================================================
# Data Integration (CCA)
# ==============================================================================

combined_obj <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
rm(obj.list)
gc()

combined_obj[["RNA"]] <- JoinLayers(combined_obj[["RNA"]])

combined_obj[["RNA"]] <- split(combined_obj[["RNA"]], f = combined_obj$dataset)

combined_obj <- ScaleData(combined_obj, features = vfs, verbose = FALSE)
combined_obj <- RunPCA(combined_obj, features = vfs, verbose = FALSE)

combined_obj <- IntegrateLayers(
  object = combined_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "cca",
  features = vfs,
  dims = 1:30,
  verbose = TRUE
)

combined_obj <- RunUMAP(
  combined_obj,
  reduction = "cca",
  dims = 1:30,
  reduction.name = "cca.umap"
)

# ==============================================================================
# Visualization: UMAP by Dataset
# ==============================================================================

combined_obj$celltype <- ifelse(
  !is.na(combined_obj$final.celltype),
  as.character(combined_obj$final.celltype),
  as.character(combined_obj$clusters)
)

my_colors <- readRDS("/data/xiaotx/tools/celltype_palette_41.rds")

my_colors <- sample(my_colors)

n_celltypes <- length(levels(combined_obj$celltype))

ArchR_Palettes <- readRDS("/data/xiaotx/tools/ArchR_Palettes.RDS")

chosen_palette_name <- "stallion"
my_palette_raw <- ArchR_Palettes[[chosen_palette_name]]

my_palette <- unname(my_palette_raw)

if (is.factor(combined_obj$dataset)) {
  cell_type_names <- levels(combined_obj$dataset)
} else {
  cell_type_names <- sort(unique(as.character(combined_obj$dataset)))
}
num_celltypes <- length(cell_type_names)

if (length(my_palette) < num_celltypes) {
  cat(paste("Extending palette '", chosen_palette_name, "' to ", num_celltypes, " colors.\n", sep = ""))
  my_palette <- colorRampPalette(my_palette)(num_celltypes)
}

plot_colors <- my_palette[1:num_celltypes]
names(plot_colors) <- cell_type_names
p1_pure <- DimPlot(combined_obj,
  label = FALSE, raster = FALSE,
  reduction = "cca.umap", group.by = "dataset",
  cols = plot_colors,
  pt.size = 0.1
) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_p1 <- file.path(output_dir, "04a_Human_ST_Integration_by_dataset.png")
ggsave(
  filename = output_filepath_p1,
  plot = p1_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p1_full <- DimPlot(combined_obj,
  label = TRUE, raster = FALSE,
  reduction = "cca.umap", group.by = "dataset",
  cols = plot_colors, repel = TRUE
) +
  NoAxes() + theme(aspect.ratio = 1, plot.title = element_blank())

g_p1 <- ggplot2::ggplotGrob(p1_full)
legend_index_p1 <- which(sapply(g_p1$grobs, function(x) x$name) == "guide-box")
if (length(legend_index_p1) > 0) {
  legend_p1 <- g_p1$grobs[[legend_index_p1[1]]]
  output_filepath_p1_legend <- file.path(output_dir, "04a_Human_ST_Integration_by_dataset_legend.pdf")
  pdf(output_filepath_p1_legend, width = 5, height = 8)
  grid::grid.draw(legend_p1)
  dev.off()
}

# ==============================================================================
# Visualization: Highlight Target Cells
# ==============================================================================

target_celltypes <- c("Midbrain-derived inhibitory", "IN1")
combined_obj$highlight_individual <- ifelse(
  combined_obj$celltype %in% target_celltypes,
  as.character(combined_obj$celltype),
  "Other"
)
order_levels <- c(target_celltypes, "Other")
color_map <- c(
  "Midbrain-derived inhibitory" = "#E41A1C",
  "IN1" = "#4DAF4A",
  "Other" = "lightgrey"
)

p_highlight_pure <- DimPlot(combined_obj,
  reduction = "cca.umap",
  group.by = "highlight_individual",
  order = order_levels,
  cols = color_map,
  label = FALSE,
  raster = FALSE,
  pt.size = 0.1
) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_highlight <- file.path(output_dir, "04b_Human_ST_Integration_highlight.png")
ggsave(
  filename = output_filepath_highlight,
  plot = p_highlight_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_highlight_full <- DimPlot(combined_obj,
  reduction = "cca.umap",
  group.by = "highlight_individual",
  order = order_levels,
  cols = color_map,
  label = TRUE, repel = TRUE,
  raster = FALSE, pt.size = 0.1
) +
  NoAxes() + theme(aspect.ratio = 1, plot.title = element_blank())

g_highlight <- ggplot2::ggplotGrob(p_highlight_full)
legend_index_highlight <- which(sapply(g_highlight$grobs, function(x) x$name) == "guide-box")
if (length(legend_index_highlight) > 0) {
  legend_highlight <- g_highlight$grobs[[legend_index_highlight[1]]]
  output_filepath_highlight_legend <- file.path(output_dir, "04b_Human_ST_Integration_highlight_legend.pdf")
  pdf(output_filepath_highlight_legend, width = 5, height = 8)
  grid::grid.draw(legend_highlight)
  dev.off()
}

# ==============================================================================
# Visualization: OTX2 Expression
# ==============================================================================

p_otx2_pure <- FeaturePlot(combined_obj,
  reduction = "cca.umap",
  features = "OTX2",
  label = FALSE,
  raster = FALSE,
  pt.size = 0.1
) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_otx2 <- file.path(output_dir, "04b_Human_ST_Integration_OTX2.png")
ggsave(
  filename = output_filepath_otx2,
  plot = p_otx2_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_otx2_full <- FeaturePlot(combined_obj,
  reduction = "cca.umap",
  features = "OTX2",
  label = FALSE,
  raster = FALSE,
  pt.size = 0.1
) +
  NoAxes() + theme(aspect.ratio = 1, plot.title = element_blank())

g_otx2 <- ggplot2::ggplotGrob(p_otx2_full)
legend_index_otx2 <- which(sapply(g_otx2$grobs, function(x) x$name) == "guide-box")
if (length(legend_index_otx2) > 0) {
  legend_otx2 <- g_otx2$grobs[[legend_index_otx2[1]]]
  output_filepath_otx2_legend <- file.path(output_dir, "04b_Human_ST_Integration_OTX2_legend.pdf")
  pdf(output_filepath_otx2_legend, width = 2, height = 4)
  grid::grid.draw(legend_otx2)
  dev.off()
}

spatial_data <- readRDS("/data/xiaotx/Th/ST/spatial_list.rds")

sample_names <- c("GW10", "GW23_L1", "GW23_L2", "GW23_L3", "GW16", "GW23_M2", "GW23_vM", "GW21", "GW23_PV", "GW21_2")

target_samples <- c("GW23_vM", "GW23_L1")
target_celltypes <- c("IN1")

for (target_sample_name in target_samples) {
  target_index <- which(sample_names == target_sample_name)
  cat("Processing sample:", target_sample_name, "at index:", target_index, "\n")
  specific_slice <- spatial_data[[target_index]]

  p_gene <- ImageFeaturePlot(specific_slice,
    features = "OTX2",
    size = 1.0
  ) +
    scale_fill_gradient(low = "lightgrey", high = "blue") +
    theme_classic() +
    NoAxes() +
    theme(
      text = element_blank(),
      aspect.ratio = 1,
      legend.position = "none"
    )

  ggsave(
    filename = file.path(output_dir, paste0("Human_", target_sample_name, "_OTX2.png")),
    plot = p_gene,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )

  specific_slice$highlight_individual <- ifelse(
    specific_slice$clusters %in% target_celltypes,
    as.character(specific_slice$clusters),
    "Other"
  )

  p_cell <- ImageDimPlot(specific_slice,
    group.by = "highlight_individual",
    cols = c("grey90", "red"),
    size = 1.0
  ) +
    theme_classic() +
    NoAxes() +
    theme(
      text = element_blank(),
      aspect.ratio = 1,
      legend.position = "none"
    )

  ggsave(
    filename = file.path(output_dir, paste0("Human_", target_sample_name, "_IN1_Highlight.png")),
    plot = p_cell,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )
}

names(spatial_data) <- sample_names

output_dir <- "/data/xiaotx/Th/Figure/04_Human_ST/Spatial_OTX2_Module"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

target_samples <- c("GW23_vM", "GW23_L1")

positive_genes <- c(
  "OTX2", "CHRM2", "SERPINF1", "LMO1", "GRIK1", "LHX1",
  "GAD2", "GATA3", "DMBX1", "TMEM145", "ASIC4",
  "TMEM260", "POU2F1", "SYNPR"
)

score_name_base <- "OTX2_Network_Score"
score_name_meta <- paste0(score_name_base, "1")

for (i in seq_along(spatial_data)) {
  current_slice <- spatial_data[[i]]
  valid_genes <- positive_genes[positive_genes %in% rownames(current_slice)]

  if (length(valid_genes) > 0) {
    spatial_data[[i]] <- AddModuleScore(
      object = current_slice,
      features = list(valid_genes),
      name = score_name_base,
      ctrl = 5
    )
  } else {
    spatial_data[[i]]@meta.data[[score_name_meta]] <- 0
    warning(paste("Slice", names(spatial_data)[i], "has no valid score genes."))
  }
}

features_to_plot <- c("OTX2", score_name_meta)
fix_color <- "red"

for (target_sample_name in target_samples) {
  if (!target_sample_name %in% names(spatial_data)) {
    next
  }

  cat(sprintf("\nProcessing sample: %s\n", target_sample_name))
  specific_slice <- spatial_data[[target_sample_name]]

  for (feature in features_to_plot) {
    is_score <- feature == score_name_meta

    if (!is_score && !(feature %in% rownames(specific_slice))) {
      cat(sprintf("  Skip: %s\n", feature))
      next
    }

    if (is_score) {
      file_tag <- "OTX2_Network_Module"

      plot_layer <- scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = c(-2, 2),
        oob = scales::squish
      )
    } else {
      file_tag <- paste0(feature, "_Expression")

      plot_layer <- scale_fill_gradient(
        low = "lightgrey",
        high = "red",
        limits = c(0, 6),
        oob = scales::squish
      )
    }

    p_base <- ImageFeaturePlot(specific_slice,
      features = feature,
      size = 1.0
    ) +
      plot_layer +
      theme_classic() +
      theme(
        text = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        aspect.ratio = 1
      )

    p_pure <- p_base + NoLegend()

    ggsave(
      filename = file.path(output_dir, paste0("Human_", target_sample_name, "_", file_tag, "_pure.png")),
      plot = p_pure,
      width = 6, height = 6, dpi = 600, units = "in"
    )

    g_obj <- ggplot2::ggplotGrob(p_base)
    legend_index <- which(sapply(g_obj$grobs, function(x) x$name) == "guide-box")

    if (length(legend_index) > 0) {
      legend_content <- g_obj$grobs[[legend_index[1]]]

      pdf(file.path(output_dir, paste0("Human_", target_sample_name, "_", file_tag, "_Legend.pdf")), width = 2, height = 4)
      grid::grid.draw(legend_content)
      dev.off()
    }
  }
}
