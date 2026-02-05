setwd("/data/xiaotx/Th")

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(data.table)
library(Cairo)
library(httpgd)

set.seed(1111)

options(future.globals.maxSize = +Inf)

# ==============================================================================
# Load Data
# ==============================================================================

th_final <- readRDS("/data/xiaotx/Th/Object/Th_final.rds")

human <- readRDS("/data/xiaotx/Th/BICCN/Human_joinlayer_re_trans.RDS")

human_th <- subset(human, subset = ROIGroupFine == "Thalamus")

non_thalamic_types <- c(
  "Bergmann glia", "Cerebellar inhibitory", "Granule",
  "Deep-layer excitatory", "Upper-layer excitatory",
  "Medium spiny neuron", "Eccentric medium spiny neuron",
  "Hippocampal CA1-3", "Hippocampal CA4", "Hippocampal dentate gyrus",
  "Deep-layer corticothalamic and 6b", "LAMP5-LHX6 and Chandelier",
  "Upper rhombic lip", "Upper-layer intratelencephalic"
  # "MGE interneuron", "CGE interneuron"
)

human_th <- subset(human_th,
                   subset = predicted.id %in% non_thalamic_types,
                   invert = TRUE)

human_th$final.celltype <- human_th$predicted.id

ncol(human_th)

FeaturePlot(human_th, features = "prediction.score.max",
            reduction = "ref.umap")

human_th <- subset(human_th, subset = prediction.score.max > 0.6)

# ==============================================================================
# Data Integration
# ==============================================================================

if (FALSE) {
  macaque <- subset(th_final, subset = Species == "Macaque")

  DimPlot(human_th_only, reduction = "ref.pca", label = TRUE, repel = TRUE)

  DimPlot(macaque, reduction = "ref.pca", label = TRUE, repel = TRUE)

  human_embeddings <- human_th_only@reductions$ref.pca@cell.embeddings
  macaque_embeddings <- macaque@reductions$ref.pca@cell.embeddings

  print(dim(human_embeddings))
  print(dim(macaque_embeddings))

  combined_embeddings <- rbind(human_embeddings, macaque_embeddings)

  merged_obj <- merge(macaque, y = human_th_only)

  ordered_embeddings <- combined_embeddings[colnames(merged_obj), ]

  print(all(rownames(ordered_embeddings) == colnames(merged_obj)))

  ref_pca_loadings <- human_th_only@reductions$ref.pca@feature.loadings
  ref_pca_stdev <- human_th_only@reductions$ref.pca@stdev

  merged_obj@reductions$ref.pca <- CreateDimReducObject(
    embeddings = ordered_embeddings,
    loadings = ref_pca_loadings,
    stdev = ref_pca_stdev,
    key = "refPC_",
    assay = "RNA"
  )

  print(merged_obj)

  merged_obj <- RunUMAP(merged_obj, reduction = "ref.pca", dims = 1:30)
  DimPlot(merged_obj, group.by = "predicted.id")
}

th_test <- th_final

DimPlot(human_th, reduction = "ref.pca", label = TRUE, repel = TRUE)

DimPlot(th_final, reduction = "ref.pca", label = TRUE, repel = TRUE)

human_embeddings <- human_th@reductions$ref.pca@cell.embeddings
th_embeddings <- th_final@reductions$ref.pca@cell.embeddings

print(dim(human_embeddings))
print(dim(th_embeddings))

combined_embeddings <- rbind(human_embeddings, th_embeddings)

merged_obj <- merge(th_final, y = human_th)

ordered_embeddings <- combined_embeddings[colnames(merged_obj), ]

print(all(rownames(ordered_embeddings) == colnames(merged_obj)))

ref_pca_loadings <- human_th@reductions$ref.pca@feature.loadings
ref_pca_stdev <- human_th@reductions$ref.pca@stdev

merged_obj@reductions$ref.pca <- CreateDimReducObject(
  embeddings = ordered_embeddings,
  loadings = ref_pca_loadings,
  stdev = ref_pca_stdev,
  key = "refPC_",
  assay = "RNA"
)

print(merged_obj)

merged_obj[["RNA"]]$counts.Human.1 <- merged_obj[["RNA"]]$counts.2
merged_obj[["RNA"]]$data.Human.1 <- merged_obj[["RNA"]]$data.2
merged_obj[["RNA"]]$scale.data.Human.1 <- merged_obj[["RNA"]]$scale.data.2
merged_obj[["RNA"]]$counts.2 <- NULL
merged_obj[["RNA"]]$data.2 <- NULL
merged_obj[["RNA"]]$scale.data.2 <- NULL

# ==============================================================================
# Visualization and Clustering
# ==============================================================================

print(merged_obj)

merged_try <- merged_obj

merged_obj <- RunUMAP(merged_obj,
                      reduction = "ref.pca",
                      dims = 1:30,
                      reduction.name = "ref.umap",
                      return.model = TRUE,
                      verbose = TRUE,
                      min.dist = 0.25,
                      n.neighbors = 30)

DimPlot(merged_obj, group.by = "final.celltype", reduction = "ref.umap",
        label = FALSE,
        raster = FALSE) +
  theme(aspect.ratio = 1)

merged_obj <- FindNeighbors(merged_obj, reduction = "ref.pca", dims = 1:30)

merged_obj <- FindClusters(merged_obj, resolution = 0.8,
                           random.seed = 1111, algorithm = 4, verbose = TRUE)

saveRDS(merged_obj, file = "/data/xiaotx/Th/Object/Th_final_human.rds")

non_neuron_markers <- c(
  "AQP4",   # Astrocytes
  "MOG",   # Oligodendrocytes
  "PDGFRA", # OPC
  "PTPRC",  # Microglia
  "CFAP126",  # Ependymal Cells
  "CLDN5",  # Vascular Cells
  "COL1A2",     # Fibroblasts
  "SNAP25"
)

FeaturePlot(
  merged_obj,
  features = non_neuron_markers,
  reduction = "ref.umap",
  repel = TRUE,
  label = TRUE,
  order = TRUE,
  ncol = 4
) & theme(aspect.ratio = 1)

DimPlot(merged_obj, group.by = "seurat_clusters", reduction = "ref.umap",
# ==============================================================================
# Cluster Annotation
# ==============================================================================

        label = TRUE,
        raster = TRUE) +
  theme(aspect.ratio = 1)

annotation_map <- c(
  "1" = "Oligodendrocyte",
  "2" = "Oligodendrocyte",
  "3" = "Thalamic excitatory",
  "4" = "Astrocyte",
  "5" = "Thalamic excitatory",
  "6" = "Oligodendrocyte precursor",
  "7" = "Microglia",
  "8" = "Splatter",
  "9" = "Thalamic excitatory",
  "10" = "Thalamic excitatory",
  "11" = "Astrocyte",
  "12" = "Midbrain-derived inhibitory",
  "13" = "Thalamic excitatory",
  "14" = "Oligodendrocyte",
  "15" = "Thalamic excitatory",
  "16" = "Ependymal",
  "17" = "Thalamic excitatory",
  "18" = "Splatter",
  "19" = "Vascular",
  "20" = "Thalamic excitatory",
  "21" = "Fibroblast",
  "22" = "Thalamic excitatory",
  "23" = "Oligodendrocyte",
  "24" = "Thalamic excitatory",
  "25" = "Microglia",
  "26" = "Thalamic excitatory",
  "27" = "CGE/MGE inhibitory",
  "28" = "Thalamic excitatory",
  "29" = "Oligodendrocyte",
  "30" = "Astrocyte",
  "31" = "Oligodendrocyte",
  "32" = "Oligodendrocyte",
  "33" = "Oligodendrocyte",
  "34" = "Oligodendrocyte",
  "35" = "Oligodendrocyte",
  "36" = "Oligodendrocyte",
  "37" = "Oligodendrocyte",
  "38" = "Oligodendrocyte",
  "39" = "Oligodendrocyte",
  "40" = "Oligodendrocyte",
  "41" = "Oligodendrocyte",
  "42" = "Oligodendrocyte",
  "43" = "Oligodendrocyte",
  "44" = "Oligodendrocyte",
  "45" = "Oligodendrocyte",
  "46" = "Oligodendrocyte",
  "47" = "Oligodendrocyte",
  "48" = "Oligodendrocyte",
  "49" = "Oligodendrocyte",
  "50" = "Oligodendrocyte",
  "51" = "Oligodendrocyte",
  "52" = "Oligodendrocyte",
  "53" = "Oligodendrocyte"
)

new_metadata <- data.frame(
  final.celltype = annotation_map[as.character(merged_obj$seurat_clusters)],
  row.names = colnames(merged_obj)
)

merged_obj <- AddMetaData(merged_obj, metadata = new_metadata)

DimPlot(merged_obj, reduction = "ref.umap", group.by = "final.celltype",
        label = TRUE, repel = TRUE) +
  labs(title = "UMAP with Manual Celltype") + theme(aspect.ratio = 1)

# ==============================================================================
# Sub-clustering Inhibitory Neurons
# ==============================================================================

merged_test <- merged_obj

merged_test <- JoinLayers(merged_test)

inhibitory_subset <- subset(merged_test,
                            subset = final.celltype == "CGE/MGE inhibitory")

inhibitory_subset <- RunUMAP(inhibitory_subset,
                             reduction = "ref.pca",
                             dims = 1:30, verbose = FALSE)

inhibitory_subset <- FindNeighbors(inhibitory_subset,
                                   reduction = "ref.pca",
                                   dims = 1:30, verbose = FALSE)

inhibitory_subset <- FindClusters(inhibitory_subset,
                                  resolution = 0.5, verbose = FALSE)

DimPlot(inhibitory_subset, reduction = "umap",
        label = TRUE, group.by = "seurat_clusters") +
  labs(title = "Re-clustered CGE/MGE Neurons") + theme(aspect.ratio = 1)

FeaturePlot(
  inhibitory_subset,
  features = c("SST", "LHX6", "VIP", "LAMP5"),
  reduction = "umap",
  ncol = 2
) & theme(aspect.ratio = 1)

sub_annotation_map <- c(
  "0" = "CGE interneuron",
  "1" = "CGE interneuron",
  "2" = "MGE interneuron",
  "3" = "MGE interneuron",
  "4" = "MGE interneuron",
  "5" = "CGE interneuron",
  "6" = "MGE interneuron"
)

inhibitory_subset$inhibitory.subtype <- plyr::mapvalues(
  x = inhibitory_subset$seurat_clusters,
  from = names(sub_annotation_map),
  to = unname(sub_annotation_map)
)

DimPlot(inhibitory_subset, group.by = "inhibitory.subtype",
        reduction = "umap", label = TRUE) +
  labs(title = "Precise Annotations within Inhibitory Subset") + 
  theme(aspect.ratio = 1)

new_subtypes_with_names <- as.character(inhibitory_subset$inhibitory.subtype)
names(new_subtypes_with_names) <- colnames(inhibitory_subset)

cells_to_update <- colnames(inhibitory_subset)

merged_obj$final.celltype[cells_to_update] <- new_subtypes_with_names

saveRDS(merged_obj, file = "/data/xiaotx/Th/Object/Th_final_human.rds")

## ArchR_Palettes <- readRDS("/data/xiaotx/tools/ArchR_Palettes.RDS")

chosen_palette_name <- "stallion"
my_palette_raw <- ArchR_Palettes[[chosen_palette_name]]

my_palette <- unname(my_palette_raw)

if (is.factor(merged_obj$final.celltype)) {
  cell_type_names <- levels(merged_obj$final.celltype)
} else {
  cell_type_names <- sort(unique(as.character(merged_obj$final.celltype)))
}
num_celltypes <- length(cell_type_names)

if (length(my_palette) < num_celltypes) {
  cat(paste("Extending palette '",
            chosen_palette_name, "' to ", num_celltypes, " colors.\n", sep=""))
  my_palette <- colorRampPalette(my_palette)(num_celltypes)
}

plot_colors <- my_palette[1:num_celltypes]
### names(plot_colors) <- cell_type_names

plot_colors <- readRDS("./Data/ct_col_new.RDS")

species_colors <- readRDS("./Data/species_color_38_new.rds")

names(species_colors) <- gsub("_", " ", names(species_colors))

missing_species <- setdiff(unique(merged_obj$latin_name), names(species_colors))

if (length(missing_species) > 0) {
  cat("以下物种在调色板中缺失:\n")
  print(missing_species)
} else {
  cat("所有物种均已匹配。\n")
}

output_dir <- "/data/xiaotx/Th/Figure/01_Umap"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

p_pure <- DimPlot(merged_obj,
                  group.by = "latin_name",
                  reduction = "ref.umap",
                  label = FALSE,
                  raster = FALSE,
                  cols = species_colors) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_pure <- file.path(output_dir, "Umap_by_species_pure.png")
ggsave(
  filename = output_filepath_pure,
  plot = p_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_pure <- DimPlot(merged_obj,
                  group.by = "final.celltype",
                  reduction = "ref.umap",
                  label = FALSE,
                  raster = FALSE,
                  cols = plot_colors) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_pure <- file.path(output_dir, "Umap_by_final.celltype_pure.png")
ggsave(
  filename = output_filepath_pure,
  plot = p_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_for_legend <- DimPlot(merged_obj,
                        group.by = "final.celltype",
                        reduction = "ref.umap",
                        label = TRUE,
                        repel = TRUE,
                        raster = FALSE,
                        cols = plot_colors) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g <- ggplot2::ggplotGrob(p_for_legend)
legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")

if(length(legend_index) > 0) {
  legend_grob <- g$grobs[[legend_index[1]]]

  output_filepath_legend <- file.path(output_dir,
                                      "Umap_by_final.celltype_legend_new.pdf")
  pdf(output_filepath_legend, width = 5, height = 8)
  grid::grid.draw(legend_grob)
  dev.off()
}

merged_plot <- merged_obj

merged_plot <- JoinLayers(merged_plot)

merged_plot <- NormalizeData(merged_plot)
merged_plot <- ScaleData(merged_plot)

merged_plot$latin_name[merged_plot$species == "Human"] <- "Homo sapiens"

saveRDS(merged_plot, "/data/xiaotx/Th/Object/Final_human_joinlayers.rds")

output_dir <- "/data/xiaotx/Th/Figure/02_FeaturePlot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

p_all_pure <- FeaturePlot(merged_plot, features = "OTX2",
                          reduction = "ref.umap",
                          label = FALSE,
                          raster = FALSE,
                          max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_all_pure <- file.path(output_dir,
                                      "OTX2_Expression_All_pure.png")
ggsave(
  filename = output_filepath_all_pure,
  plot = p_all_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_all_for_legend <- FeaturePlot(merged_plot, features = "OTX2",
                                reduction = "ref.umap",
                                label = FALSE,
                                raster = FALSE,
                                max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_all <- ggplot2::ggplotGrob(p_all_for_legend)
legend_index_all <- which(sapply(g_all$grobs,
                                 function(x) x$name) == "guide-box")
if (length(legend_index_all) > 0) {
  legend_all <- g_all$grobs[[legend_index_all[1]]]
  output_filepath_all_legend <- file.path(output_dir,
                                          "OTX2_Expression_All_legend.pdf")
  pdf(output_filepath_all_legend, width = 2, height = 4)
  grid::grid.draw(legend_all)
  dev.off()
}

human_seurat <- subset(merged_plot, subset = species == "Human")

p_human_pure <- FeaturePlot(human_seurat, features = "OTX2",
                            reduction = "ref.umap",
                            label = FALSE,
                            raster = FALSE,
                            max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_human_pure <- file.path(output_dir,
                                        "OTX2_Expression_Human_pure.png")
ggsave(
  filename = output_filepath_human_pure,
  plot = p_human_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_human_for_legend <- FeaturePlot(human_seurat, features = "OTX2",
                                  reduction = "ref.umap",
                                  label = FALSE,
                                  raster = FALSE,
                                  max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_human <- ggplot2::ggplotGrob(p_human_for_legend)
legend_index_human <- which(sapply(g_human$grobs,
                                   function(x) x$name) == "guide-box")
if (length(legend_index_human) > 0) {
  legend_human <- g_human$grobs[[legend_index_human[1]]]
  output_filepath_human_legend <- file.path(output_dir,
                                            "OTX2_Expression_Human_legend.pdf")
  pdf(output_filepath_human_legend, width = 2, height = 4)
  grid::grid.draw(legend_human)
  dev.off()
}

macaque_seurat <- subset(merged_plot, subset = species == "Macaque")

p_macaque_pure <- FeaturePlot(macaque_seurat, features = "OTX2",
                              reduction = "ref.umap",
                              label = FALSE,
                              raster = FALSE,
                              max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_macaque_pure <- file.path(output_dir,
                                          "OTX2_Expression_Macaque_pure.png")
ggsave(
  filename = output_filepath_macaque_pure,
  plot = p_macaque_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_macaque_for_legend <- FeaturePlot(macaque_seurat, features = "OTX2",
                                    reduction = "ref.umap",
                                    label = FALSE,
                                    raster = FALSE,
                                    max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_macaque <- ggplot2::ggplotGrob(p_macaque_for_legend)
legend_index_macaque <- which(sapply(g_macaque$grobs,
                                     function(x) x$name) == "guide-box")
if(length(legend_index_macaque) > 0) {
  legend_macaque <- g_macaque$grobs[[legend_index_macaque[1]]]
  output_filepath_macaque_legend <- file.path(output_dir,
                                              "OTX2_Macaque_legend.pdf")
  pdf(output_filepath_macaque_legend, width = 2, height = 4)
  grid::grid.draw(legend_macaque)
  dev.off()
}

mouse_seurat <- subset(merged_plot, subset = species == "Mouse")

p_mouse_pure <- FeaturePlot(mouse_seurat, features = "OTX2",
                            reduction = "ref.umap",
                            label = FALSE,
                            raster = FALSE,
                            max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_mouse_pure <- file.path(output_dir,
                                        "OTX2_Expression_Mouse_pure.png")
ggsave(
  filename = output_filepath_mouse_pure,
  plot = p_mouse_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_mouse_for_legend <- FeaturePlot(mouse_seurat, features = "OTX2",
                                  reduction = "ref.umap",
                                  label = FALSE,
                                  raster = FALSE,
                                  max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_mouse <- ggplot2::ggplotGrob(p_mouse_for_legend)
legend_index_mouse <- which(sapply(g_mouse$grobs,
                                   function(x) x$name) == "guide-box")
if (length(legend_index_mouse) > 0) {
  legend_mouse <- g_mouse$grobs[[legend_index_mouse[1]]]
  output_filepath_mouse_legend <- file.path(output_dir,
                                            "OTX2_Expression_Mouse_legend.pdf")
  pdf(output_filepath_mouse_legend, width = 2, height = 4)
  grid::grid.draw(legend_mouse)
  dev.off()
}


p_all_pure <- FeaturePlot(merged_plot, features = "GATA3",
                          reduction = "ref.umap",
                          label = FALSE,
                          raster = FALSE,
                          max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_all_pure <- file.path(output_dir,
                                      "GATA3_Expression_All_pure.png")
ggsave(
  filename = output_filepath_all_pure,
  plot = p_all_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_all_for_legend <- FeaturePlot(merged_plot, features = "GATA3",
                                reduction = "ref.umap",
                                label = FALSE,
                                raster = FALSE,
                                max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_all <- ggplot2::ggplotGrob(p_all_for_legend)
legend_index_all <- which(sapply(g_all$grobs,
                                 function(x) x$name) == "guide-box")
if (length(legend_index_all) > 0) {
  legend_all <- g_all$grobs[[legend_index_all[1]]]
  output_filepath_all_legend <- file.path(output_dir,
                                          "GATA3_Expression_All_legend.pdf")
  pdf(output_filepath_all_legend, width = 2, height = 4)
  grid::grid.draw(legend_all)
  dev.off()
}

human_seurat <- subset(merged_plot, subset = species == "Human")

p_human_pure <- FeaturePlot(human_seurat, features = "GATA3",
                            reduction = "ref.umap",
                            label = FALSE,
                            raster = FALSE,
                            max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_human_pure <- file.path(output_dir,
                                        "GATA3_Expression_Human_pure.png")
ggsave(
  filename = output_filepath_human_pure,
  plot = p_human_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_human_for_legend <- FeaturePlot(human_seurat, features = "GATA3",
                                  reduction = "ref.umap",
                                  label = FALSE,
                                  raster = FALSE,
                                  max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_human <- ggplot2::ggplotGrob(p_human_for_legend)
legend_index_human <- which(sapply(g_human$grobs,
                                   function(x) x$name) == "guide-box")
if (length(legend_index_human) > 0) {
  legend_human <- g_human$grobs[[legend_index_human[1]]]
  output_filepath_human_legend <- file.path(output_dir,
                                            "GATA3_Expression_Human_legend.pdf")
  pdf(output_filepath_human_legend, width = 2, height = 4)
  grid::grid.draw(legend_human)
  dev.off()
}

macaque_seurat <- subset(merged_plot, subset = species == "Macaque")

p_macaque_pure <- FeaturePlot(macaque_seurat, features = "GATA3",
                              reduction = "ref.umap",
                              label = FALSE,
                              raster = FALSE,
                              max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_macaque_pure <- file.path(output_dir,
                                          "GATA3_Expression_Macaque_pure.png")
ggsave(
  filename = output_filepath_macaque_pure,
  plot = p_macaque_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_macaque_for_legend <- FeaturePlot(macaque_seurat, features = "GATA3",
                                    reduction = "ref.umap",
                                    label = FALSE,
                                    raster = FALSE,
                                    max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_macaque <- ggplot2::ggplotGrob(p_macaque_for_legend)
legend_index_macaque <- which(sapply(g_macaque$grobs,
                                     function(x) x$name) == "guide-box")
if(length(legend_index_macaque) > 0) {
  legend_macaque <- g_macaque$grobs[[legend_index_macaque[1]]]
  output_filepath_macaque_legend <- file.path(output_dir,
                                              "GATA3_Macaque_legend.pdf")
  pdf(output_filepath_macaque_legend, width = 2, height = 4)
  grid::grid.draw(legend_macaque)
  dev.off()
}

mouse_seurat <- subset(merged_plot, subset = species == "Mouse")

p_mouse_pure <- FeaturePlot(mouse_seurat, features = "GATA3",
                            reduction = "ref.umap",
                            label = FALSE,
                            raster = FALSE,
                            max.cutoff = 2) +
  NoAxes() +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_mouse_pure <- file.path(output_dir,
                                        "GATA3_Expression_Mouse_pure.png")
ggsave(
  filename = output_filepath_mouse_pure,
  plot = p_mouse_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_mouse_for_legend <- FeaturePlot(mouse_seurat, features = "GATA3",
                                  reduction = "ref.umap",
                                  label = FALSE,
                                  raster = FALSE,
                                  max.cutoff = 2) +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

g_mouse <- ggplot2::ggplotGrob(p_mouse_for_legend)
legend_index_mouse <- which(sapply(g_mouse$grobs,
                                   function(x) x$name) == "guide-box")
if (length(legend_index_mouse) > 0) {
  legend_mouse <- g_mouse$grobs[[legend_index_mouse[1]]]
  output_filepath_mouse_legend <- file.path(output_dir,
                                            "GATA3_Expression_Mouse_legend.pdf")
  pdf(output_filepath_mouse_legend, width = 2, height = 4)
  grid::grid.draw(legend_mouse)
  dev.off()
}


merged_try <- subset(merged_try,
                     subset = final.celltype %in% c("CGE interneuron",
                                                    "MGE interneuron"),
                     invert = TRUE)

merged_try <- RunUMAP(merged_try,
                      reduction = "ref.pca",
                      dims = 1:30,
                      reduction.name = "ref.umap",
                      return.model = TRUE,
                      verbose = TRUE,
                      min.dist = 0.25,
                      n.neighbors = 30)

DimPlot(merged_try, group.by = "final.celltype", reduction = "ref.umap",
        label = FALSE,
        raster = FALSE,
        cols = plot_colors) +
  theme(aspect.ratio = 1)

merged_try <- FindNeighbors(merged_try, reduction = "ref.pca", dims = 1:30)

merged_try <- FindClusters(merged_try, resolution = 0.8,
                           random.seed = 1111, algorithm = 4, verbose = TRUE)

saveRDS(merged_try, file = "/data/xiaotx/Th/Object/Th_final_human_try.rds")

non_neuron_markers <- c(
  "AQP4",   # Astrocytes
  "MOG",   # Oligodendrocytes
  "PDGFRA", # OPC
  "PTPRC",  # Microglia
  "CFAP126",  # Ependymal Cells
  "CLDN5",  # Vascular Cells
  "COL1A2",     # Fibroblasts
  "SNAP25"
)

FeaturePlot(
  merged_try,
  features = non_neuron_markers,
  reduction = "ref.umap",
  repel = TRUE,
  label = TRUE,
  order = TRUE,
  ncol = 4
) & theme(aspect.ratio = 1)



ArchR_Palettes <- readRDS("/data/xiaotx/tools/ArchR_Palettes.RDS")

chosen_palette_name <- "stallion"
my_palette_raw <- ArchR_Palettes[[chosen_palette_name]]

my_palette <- unname(my_palette_raw)

if (is.factor(merged_try$final.celltype)) {
  cell_type_names <- levels(merged_try$final.celltype)
} else {
  cell_type_names <- sort(unique(as.character(merged_try$final.celltype)))
}
num_celltypes <- length(cell_type_names)

if (length(my_palette) < num_celltypes) {
  cat(paste("Extending palette '", chosen_palette_name, "' to ", num_celltypes, " colors.\n", sep=""))
  my_palette <- colorRampPalette(my_palette)(num_celltypes)
}

plot_colors <- my_palette[1:num_celltypes]
names(plot_colors) <- cell_type_names
