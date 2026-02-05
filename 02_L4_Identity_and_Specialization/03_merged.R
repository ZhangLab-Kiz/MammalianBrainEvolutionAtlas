# ==============================================================================
# Step 03: Merge Reference and Query
# ==============================================================================

setwd("/home/xiaotx/WholeBrain")
library(Seurat)
library(ggplot2)

DIR_OBJ <- "./Query/Query_object"
DIR_FIG <- "./Query/Query_Figure"

# ==============================================================================
# Load Data
# ==============================================================================

if (!exists("human_ref")) {
  cat("Loading Human Reference...\n")
  human_ref <- readRDS(file.path(DIR_OBJ, "Human_sobj_cortex.RDS"))
}

if (!exists("seu_query")) {
  cat("Loading Mapped Query...\n")
  seu_query <- readRDS(file.path(DIR_OBJ, "All_Query_Mapped.rds"))
}

# ==============================================================================
# Align Reductions
# ==============================================================================

# Align PCA
query_pca_key <- Key(seu_query[["ref.pca"]])
n_dims <- ncol(Embeddings(seu_query, "ref.pca"))
human_emb <- Embeddings(human_ref, "pca")[, 1:n_dims]
colnames(human_emb) <- paste0(query_pca_key, 1:n_dims)
human_ref[["ref.pca"]] <- CreateDimReducObject(
  embeddings = human_emb,
  key = query_pca_key,
  assay = "RNA"
)

# Align UMAP
query_umap_key <- Key(seu_query[["ref.umap"]])
human_umap_emb <- Embeddings(human_ref, "umap")
colnames(human_umap_emb) <- paste0(query_umap_key, 1:ncol(human_umap_emb))
human_ref[["ref.umap"]] <- CreateDimReducObject(
  embeddings = human_umap_emb,
  key = query_umap_key,
  assay = "RNA"
)

# Unify Metadata Columns
if ("cell" %in% colnames(human_ref@meta.data)) {
  human_ref$final.celltype <- human_ref$cell
} else {
  stop("Error: 'cell column not found in Human Reference!")
}

query_cols <- colnames(seu_query@meta.data)
if ("predicted.predicted.id" %in% query_cols) {
  seu_query$final.celltype <- seu_query$predicted.predicted.id
  cat("Used 'predicted.predicted.id' for Query.\n")
} else {
  stop("Error: Could not find prediction column in Query!")
}

# ==============================================================================
# Merge Objects
# ==============================================================================

obj_integrated <- merge(
  x = human_ref,
  y = seu_query,
  merge.dr = c("ref.umap", "ref.pca")
)

obj_integrated[["RNA"]]$scale.data.1 <- NULL
obj_integrated[["RNA"]] <- JoinLayers(obj_integrated[["RNA"]])

rm(human_ref, seu_query)
gc()

# ==============================================================================
# Quality Control
# ==============================================================================

# Assign score 0.99 to Human reference cells
obj_integrated$predicted.predicted.id.score[obj_integrated$orig.ident == "Human"] <- 0.99
obj_filtered <- subset(obj_integrated, subset = predicted.predicted.id.score > 0.5)
Idents(obj_filtered) <- "final.celltype"
umap_coords <- Embeddings(obj_filtered, reduction = "ref.umap")

k_param <- 20
knn_res <- RANN::nn2(data = umap_coords, query = umap_coords, k = k_param + 1)
mean_dist <- rowMeans(knn_res$nn.dists[, 2:(k_param + 1)])
names(mean_dist) <- rownames(umap_coords) 
obj_filtered$umap_local_density_score <- mean_dist
threshold <- quantile(mean_dist, probs = 0.998) 
cells_to_keep <- mean_dist < threshold
n_removed <- sum(!cells_to_keep)
obj_filtered <- subset(obj_filtered, cells = names(cells_to_keep)[cells_to_keep])

# ==============================================================================
# Save Final Integrated Atlas
# ==============================================================================
save_path <- file.path(DIR_OBJ, "All_Query_Mapped_Final_0128.rds")
saveRDS(obj_filtered, save_path)
