library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(patchwork)
library(sceasy)
library(reticulate)
library(ggridges)
library(readr)
library(gridExtra)
library(grid)

Sys.unsetenv("LD_LIBRARY_PATH")

library(reticulate)

use_python("/data/xiaotx/miniconda3/bin/python", required = TRUE)

py_config()

# ==============================================================================
# Load Mouse ST Data
# ==============================================================================

h5ad_file_path <- "/data/xiaotx/ABC_Atlas/whole_brain/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw.h5ad"
query <- sceasy::convertFormat(h5ad_file_path, from = "anndata", to = "seurat")

print(query)

metadata_path <- "/data/xiaotx/ABC_Atlas/whole_brain/metadata/MERFISH-C57BL6J-638850/20241115/views/cell_metadata_with_cluster_annotation.csv"
cell_metadata <- read.csv(metadata_path, row.names = 1)

seurat_cells <- colnames(query)

metadata_cells <- rownames(cell_metadata)

common_cells <- intersect(seurat_cells, metadata_cells)

query <- subset(query, cells = common_cells)

cell_metadata <- cell_metadata[common_cells, ]

query <- AddMetaData(query, metadata = cell_metadata)

# ==============================================================================
# Gene Name Mapping
# ==============================================================================

gene_map_path <- "/data/xiaotx/ABC_Atlas/whole_brain/metadata/MERFISH-C57BL6J-638850/20241115/gene.csv"
gene_map <- read.csv(gene_map_path)

current_gene_ids <- rownames(query)

match_indices <- match(current_gene_ids, gene_map$gene_identifier)

new_gene_names <- gene_map$gene_symbol[match_indices]

na_indices <- is.na(new_gene_names)

counts_matrix <- GetAssayData(query, assay = "RNA", slot = "counts")

rownames(counts_matrix) <- new_gene_names

query <- CreateSeuratObject(counts = counts_matrix, meta.data = query@meta.data)

# ==============================================================================
# Load Macaque Reference
# ==============================================================================

# Macaque SC
th_final <- readRDS("/data/xiaotx/Th/Object/Th_final.rds")

ref_macaque <- subset(th_final, subset = orig.ident == "Macaque")

rm(th_final)
gc()

# ==============================================================================
# Prepare Mouse Subsets
# ==============================================================================

# Mouse ST(4 samples)
query_subset1 <- subset(query,
                        subset = brain_section_label == "C57BL6J-638850.31")

query_subset2 <- subset(query,
                        subset = brain_section_label == "C57BL6J-638850.36")

query_subset3 <- subset(query,
                        subset = brain_section_label == "C57BL6J-638850.40")

query_subset4 <- subset(query,
                        subset = brain_section_label == "C57BL6J-638850.44")

ref_macaque[["RNA"]]$counts <- ref_macaque[["RNA"]]$counts.Macaque
ref_macaque[["RNA"]]$counts.Macaque <- NULL
ref_macaque[["RNA"]]$data.Macaque <- NULL
ref_macaque[["RNA"]]$scale.data.Macaque <- NULL

# ==============================================================================
# Homology Mapping (Mouse to Human Symbols)
# ==============================================================================

homology_table <- read_delim(
  "/data/xiaotx/tools/HMD_HumanPhenotype.rpt",
  delim = "\t", 
  col_names = c("human_symbol", "human_entrez", "mouse_symbol", "mgi_id", "source"),
  col_types = cols(.default = "c")
)

human_mouse_map <- homology_table %>%
  dplyr::select(human_symbol, mouse_symbol) %>%
  dplyr::distinct()

mouse_genes <- rownames(query_subset1)

mouse_genes <- data.frame(mouse_symbol = mouse_genes)

homologs <- mouse_genes %>%
  inner_join(human_mouse_map, by = "mouse_symbol") %>%
  dplyr::rename(MGI.symbol = mouse_symbol, HGNC.symbol = human_symbol)

head(homologs)

homologs_1_to_1 <- homologs %>%
  group_by(MGI.symbol) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  group_by(HGNC.symbol) %>%
  filter(n() == 1) %>%
  ungroup()

process_mouse_to_human <- function(mouse_obj) {
  mouse_counts <- GetAssayData(mouse_obj, layer = "counts")

  mouse_genes_to_keep <- intersect(rownames(mouse_counts), homologs_1_to_1$MGI.symbol)
  mouse_counts_filt <- mouse_counts[mouse_genes_to_keep, ]
  cat(sprintf("  - Kept %d genes with one-to-one orthologs.\n", length(mouse_genes_to_keep)))

  homologs_subset <- filter(homologs_1_to_1, MGI.symbol %in% mouse_genes_to_keep)
  gene_map <- setNames(homologs_subset$HGNC.symbol, homologs_subset$MGI.symbol)

  rownames(mouse_counts_filt) <- gene_map[rownames(mouse_counts_filt)]

  new_human_named_obj <- CreateSeuratObject(
    counts = mouse_counts_filt,
    meta.data = mouse_obj@meta.data
  )

  return(new_human_named_obj)
}

# ==============================================================================
# Data Integration (CCA)
# ==============================================================================

obj.list <- list()

obj.list[["Mouse_ST_31"]] <- process_mouse_to_human(query_subset1)
obj.list[["Mouse_ST_36"]] <- process_mouse_to_human(query_subset2)
obj.list[["Mouse_ST_40"]] <- process_mouse_to_human(query_subset3)
obj.list[["Mouse_ST_44"]] <- process_mouse_to_human(query_subset4)
rm(query_subset1, query_subset2, query_subset3, query_subset4)
gc()

obj.list[["Macaque_Sc"]] <- ref_macaque
rm(ref_macaque)
gc()

for (name in names(obj.list)) {
  obj.list[[name]]$dataset <- name
}

macaque_obj <- obj.list[["Macaque_Sc"]]
macaque_obj$dataset <- paste0(macaque_obj$dataset, "_batch", macaque_obj$batch)
obj.list[["Macaque_Sc"]] <- macaque_obj

rm(macaque_obj, homology_table, human_mouse_map_1_to_1)
gc()

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

combined_obj$celltype <- coalesce(
  combined_obj$final.celltype,
  combined_obj$class
)

my_colors <- readRDS("/data/xiaotx/tools/celltype_palette_41.rds")

my_colors <- sample(my_colors)

set.seed(1111)
new_colors <- grDevices::hcl(h = runif(3, 0, 360), c = 60, l = 65)

my_colors <- c(my_colors, new_colors)

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
  cat(paste("Extending palette '", chosen_palette_name, "' to ", num_celltypes, " colors.\n", sep=""))
  my_palette <- colorRampPalette(my_palette)(num_celltypes)
}

plot_colors <- my_palette[1:num_celltypes]
names(plot_colors) <- cell_type_names

p1 <- DimPlot(combined_obj, label = FALSE, raster = FALSE,
              reduction = "cca.umap", group.by = "dataset",
              cols = plot_colors) +
  NoAxes() + theme(aspect.ratio = 1,
                   legend.position = "none",
                   plot.title = element_blank())

output_dir <- "/data/xiaotx/Th/Figure"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_filepath_p1 <- file.path(output_dir, "05a_Mouse_ST_Integration_by_dataset.png")
ggsave(
  filename = output_filepath_p1,
  plot = p1,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p1_full <- DimPlot(combined_obj, label = TRUE, raster = FALSE,
                   reduction = "cca.umap", group.by = "dataset",
                   cols = plot_colors, repel = TRUE) +
  NoAxes() + theme(aspect.ratio = 1,
                   plot.title = element_blank())

output_filepath_p1_full <- file.path(output_dir, "05a_Mouse_ST_Integration_by_dataset_full.png")
ggsave(
  filename = output_filepath_p1_full,
  plot = p1_full,
  device = "png",
  width = 16,
  height = 12,
  units = "in",
  dpi = 600
)

g <- ggplotGrob(p1_full)
legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")

if(length(legend_index) > 0) {
  legend_p2 <- g$grobs[[legend_index[1]]]

  output_filepath_legend <- file.path(output_dir, "05a_Mouse_ST_Integration_legend_only.pdf")
  pdf(output_filepath_legend, width = 8, height = 8)
  grid.draw(legend_p2)
  dev.off()
}

target_celltypes <- c("Midbrain-derived inhibitory", "20 MB GABA")

combined_obj$highlight_individual <- ifelse(
  combined_obj$celltype %in% target_celltypes,
  as.character(combined_obj$celltype),
  "Other"
)

order_levels <- c(target_celltypes, "Other")

color_map <- c(
  "Midbrain-derived inhibitory" = "#E41A1C",
  "20 MB GABA" = "#377EB8",
  "Other" = "lightgrey"
)

p_highlight_individual <- DimPlot(combined_obj,
                                  reduction = "cca.umap",
                                  group.by = "highlight_individual",
                                  order = order_levels,
                                  cols = color_map,
                                  label = FALSE,
                                  raster = FALSE,
                                  pt.size = 0.1) +
  NoAxes() + theme(aspect.ratio = 1,
                   legend.position = "none",
                   plot.title = element_blank())

p_otx2 <- FeaturePlot(combined_obj,
                      reduction = "cca.umap",
                      features =  "OTX2",
                      label = FALSE,
                      raster = FALSE,
                      pt.size = 0.1) +
  NoAxes() +  theme(aspect.ratio = 1,
  plot.title = element_blank(),
  legend.title = element_blank(),
  legend.text = element_blank(),
  legend.key.height = unit(1.2, "cm")
)

output_filepath_highlight <- file.path(output_dir, "05b_Mouse_ST_Integration_highlight.png")
ggsave(
  filename = output_filepath_highlight,
  plot = p_highlight_individual,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

output_filepath_otx2 <- file.path(output_dir, "05b_Mouse_ST_Integration_OTX2.png")
ggsave(
  filename = output_filepath_otx2,
  plot = p_otx2,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

combined_obj <- FindNeighbors(combined_obj, dims = 1:30)

combined_obj$celltype.new <- as.character(combined_obj$celltype)

combined_obj$celltype.new[combined_obj$celltype %in% c("20 MB GABA", "Midbrain-derived inhibitory")] <- "MB_GABA_combined"

Idents(combined_obj) <- "celltype.new"

combined_obj <- FindSubCluster(combined_obj,
                               cluster = "MB_GABA_combined",
                               resolution = 0.2,
                               graph.name = "RNA_snn",
                               algorithm = 4)

mb_gaba_obj <- subset(combined_obj, celltype.new == "MB_GABA_combined")

mb_gaba_obj <- RunUMAP(mb_gaba_obj, reduction = "cca", dims = 1:30, reduction.name = "cca.umap")

chosen_palette_name <- "stallion"
my_palette_raw <- ArchR_Palettes[[chosen_palette_name]]

my_palette <- unname(my_palette_raw)

if (is.factor(mb_gaba_obj$sub.cluster)) {
  cell_type_names <- levels(mb_gaba_obj$sub.cluster)
} else {
  cell_type_names <- sort(unique(as.character(mb_gaba_obj$sub.cluster)))
}
num_celltypes <- length(cell_type_names)

if (length(my_palette) < num_celltypes) {
  cat(paste("Extending palette '", chosen_palette_name, "' to ", num_celltypes, " colors.\n", sep=""))
  my_palette <- colorRampPalette(my_palette)(num_celltypes)
}

plot_colors <- my_palette[1:num_celltypes]
names(plot_colors) <- cell_type_names

p_cluster <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "sub.cluster",
        label = TRUE, repel = TRUE,
             cols = plot_colors) +
  ggtitle("MB GABA Subclusters") +
  theme(aspect.ratio = 1)

p_celltype <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "celltype",
        label = TRUE, repel = TRUE) +
  ggtitle("MB GABA Subclusters") +
  theme(aspect.ratio = 1)

p_cluster | p_celltype

FeaturePlot(mb_gaba_obj, reduction = "cca.umap",
            features = c("OTX2", "GATA3")) +
  theme(aspect.ratio = 1)

output_dir <- "/data/xiaotx/Th/Figure"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

p_cluster_pure <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "sub.cluster",
        label = FALSE,
        cols = plot_colors,
        pt.size = 1.0) +
  NoLegend() +
  NoAxes() +
  ggtitle("") +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_cluster <- file.path(output_dir, "06a_MB_GABA_by_subcluster.png")
ggsave(
  filename = output_filepath_cluster,
  plot = p_cluster_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_cluster_full <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "sub.cluster",
        label = TRUE, repel = TRUE,
        cols = plot_colors,
        pt.size = 1.0) +
  theme(aspect.ratio = 1)

g_cluster <- ggplot2::ggplotGrob(p_cluster_full)
legend_index_cluster <- which(sapply(g_cluster$grobs, function(x) x$name) == "guide-box")
if(length(legend_index_cluster) > 0) {
  legend_cluster <- g_cluster$grobs[[legend_index_cluster[1]]]
  output_filepath_cluster_legend <- file.path(output_dir, "06a_MB_GABA_by_subcluster_legend.pdf")
  pdf(output_filepath_cluster_legend, width = 5, height = 8)
  grid::grid.draw(legend_cluster)
  dev.off()
}

my_colors <- readRDS("/data/xiaotx/tools/celltype_palette_41.rds")

my_colors <- sample(my_colors)

n_celltypes <- length(levels(mb_gaba_obj$celltype))

p_celltype_pure <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "celltype",
        label = FALSE,
        pt.size = 1.0) +
  NoLegend() +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank()) +
  scale_color_manual(values = my_colors)

output_filepath_celltype <- file.path(output_dir, "06b_MB_GABA_by_celltype.png")
ggsave(
  filename = output_filepath_celltype,
  plot = p_celltype_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_celltype_full <- DimPlot(mb_gaba_obj, reduction = "cca.umap",
        group.by = "celltype",
        label = TRUE, repel = TRUE,
        pt.size = 1.0) +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = my_colors)

g_celltype <- ggplot2::ggplotGrob(p_celltype_full)
legend_index_celltype <- which(sapply(g_celltype$grobs, function(x) x$name) == "guide-box")
if(length(legend_index_celltype) > 0) {
    legend_celltype <- g_celltype$grobs[[legend_index_celltype[1]]]
    output_filepath_celltype_legend <- file.path(output_dir, "06b_MB_GABA_by_celltype_legend.pdf")
    pdf(output_filepath_celltype_legend, width = 5, height = 8)
    grid::grid.draw(legend_celltype)
    dev.off()
}

p_otx2_gaba_pure <- FeaturePlot(mb_gaba_obj, reduction = "cca.umap",
            features = "OTX2",
            pt.size = 1.0) +
  NoLegend() +
  NoAxes() +
  theme(aspect.ratio = 1, plot.title = element_blank())

output_filepath_otx2_gaba <- file.path(output_dir, "06c_MB_GABA_OTX2.png")
ggsave(
  filename = output_filepath_otx2_gaba,
  plot = p_otx2_gaba_pure,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)

p_otx2_gaba_full <- FeaturePlot(mb_gaba_obj, reduction = "cca.umap",
            features = "OTX2",
            pt.size = 1.0) +
  theme(aspect.ratio = 1)

g_otx2 <- ggplot2::ggplotGrob(p_otx2_gaba_full)
legend_index_otx2 <- which(sapply(g_otx2$grobs, function(x) x$name) == "guide-box")
if(length(legend_index_otx2) > 0) {
    legend_otx2 <- g_otx2$grobs[[legend_index_otx2[1]]]
    output_filepath_otx2_gaba_legend <- file.path(output_dir, "06c_MB_GABA_OTX2_legend.pdf")
    pdf(output_filepath_otx2_gaba_legend, width = 2, height = 4)
    grid::grid.draw(legend_otx2)
    dev.off()
}

clusters_of_interest <- c("MB_GABA_combined_1", "MB_GABA_combined_2",
                          "MB_GABA_combined_3", "MB_GABA_combined_7")

parent_meta <- combined_obj@meta.data
final_cell_indices <- which(
  parent_meta$celltype.new == "MB_GABA_combined" &
  parent_meta$sub.cluster %in% clusters_of_interest &
  grepl("^Mouse_ST", parent_meta$dataset)
)
final_cells_to_keep <- rownames(parent_meta[final_cell_indices, ])

st_slides <- unique(combined_obj@meta.data[grepl("^Mouse_ST", combined_obj@meta.data$dataset), "dataset"])
plot_list <- list()

for (slide_id in st_slides) {

  message(paste("--- Processing slide:", slide_id, "---"))

  slide_obj_full <- subset(combined_obj, subset = dataset == slide_id)

  spatial_coords <- slide_obj_full@meta.data[, c("x", "y")]
  colnames(spatial_coords) <- c("spatial_1", "spatial_2")
  slide_obj_full[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(spatial_coords),
    key = "spatial_",
    assay = DefaultAssay(slide_obj_full)
  )

  cells_to_highlight_on_this_slide <- intersect(final_cells_to_keep, colnames(slide_obj_full))

  p <- DimPlot(
    object = slide_obj_full,
    reduction = "spatial",
    cells.highlight = cells_to_highlight_on_this_slide,
    cols.highlight = "red",
    cols = "grey80",
    sizes.highlight = 0.005,
    pt.size = 0.005,
    alpha = 0.3,
    raster = FALSE
  ) +
    scale_y_reverse() +
    ggtitle(paste(slide_id, "20 MB GABA")) +
    NoAxes() +
    theme(
      aspect.ratio = 0.6,
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  plot_list[[slide_id]] <- p
}

if (length(plot_list) > 0) {
  final_plot <- wrap_plots(plot_list, ncol = 2)
} else {
  message("Warning: No plots were generated.")
}

output_dir <- "/data/xiaotx/Th/Figure"

output_filepath <- file.path(output_dir, "05_Mouse_ST.png")

ggsave(
  filename = output_filepath,
  plot = final_plot,
  device = "png",
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)
