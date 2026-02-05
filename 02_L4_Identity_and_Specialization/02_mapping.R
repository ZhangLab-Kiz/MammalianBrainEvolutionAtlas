# ==============================================================================
# Step 02: Cross-Species Mapping
# ==============================================================================

setwd("/home/xiaotx/WholeBrain")
library(Seurat)
library(future)

plan("multicore", workers = 8)
options(future.globals.maxSize = 200 * 1024^3)

DIR_OBJ <- "./Query/Query_object"
DIR_FIG <- "./Query/Query_Figure"

FILE_REF <- file.path(DIR_OBJ, "Human_sobj_cortex.RDS")
FILE_QUERY <- file.path(DIR_OBJ, "All_Query_Species_Raw.rds")

# ==============================================================================
# Load Data
# ==============================================================================
human_ref <- readRDS(FILE_REF)

seu_query <- readRDS(FILE_QUERY)

# ==============================================================================
# Process Query
# ==============================================================================

seu_query <- NormalizeData(seu_query,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

seu_query[["RNA"]] <- JoinLayers(seu_query[["RNA"]])

# ==============================================================================
# Find Anchors & Map
# ==============================================================================

anchors <- FindTransferAnchors(
  reference = human_ref,
  query = seu_query,
  dims = 1:30,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  verbose = TRUE
)

seu_query <- MapQuery(
  anchorset = anchors,
  reference = human_ref,
  query = seu_query,
  refdata = list(predicted.id = "cell"),
  reference.reduction = "pca",
  reduction.model = "umap",
  verbose = TRUE
)

rm(anchors, human_ref)
gc()

save_path <- file.path(DIR_OBJ, "All_Query_Mapped.rds")

saveRDS(seu_query, save_path)
