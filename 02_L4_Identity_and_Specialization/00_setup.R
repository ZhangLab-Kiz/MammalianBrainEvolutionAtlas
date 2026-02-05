# ==============================================================================
# Step 00: Setup and Config
# ==============================================================================
setwd("/home/xiaotx/WholeBrain")
library(Seurat)
library(tidyverse)

target_file <- "singlecell/AllCells_include_human_umap_7.29.RDS"
OUTPUT_DIR <- "./Layer_object"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Cell Type Filtering
# ==============================================================================
seu_obj <- readRDS(target_file)

target_types <- c("Upper-layer excitatory", "Deep-layer excitatory")

if (!"predicted.id" %in% colnames(seu_obj@meta.data)) {
  stop("Column 'predicted.id' not found in meta.data")
}

seu_obj <- subset(seu_obj, subset = predicted.id %in% target_types)

gc()

print(table(seu_obj$predicted.id))

save_path_step1 <- file.path(OUTPUT_DIR, "CellType_Filtered.rds")
saveRDS(seu_obj, save_path_step1)


# ==============================================================================
# Species Filtering
# ==============================================================================

seu_obj <- readRDS(save_path_step1)

exclude_abbrev <- c(
  "Cal", "Deer", "Eri", "Glis", "Hys",
  "Pag", "Pbv", "Pro", "Rhin", "Rhiz", "Rous",
  "Scap", "Tam"
)

seu_obj <- subset(seu_obj,
  subset = orig.ident %in% exclude_abbrev,
  invert = TRUE
)
gc()

seu_obj <- JoinLayers(seu_obj)

seu_obj[["RNA"]] <- split(seu_obj[["RNA"]], f = seu_obj$orig.ident)

save_path_step2 <- file.path(OUTPUT_DIR, "Species_Filtered_Final.rds")
saveRDS(seu_obj, save_path_step2)

seu_obj <- readRDS(save_path_step2)

# ==============================================================================
# Split Data
# ==============================================================================

INPUT_FILE <- "./Layer_object/Species_Filtered_Final.rds"
DIR_OBJ <- "./Query/Query_object"
DIR_FIG <- "./Query/Query_Figure"

if (!dir.exists(DIR_OBJ)) dir.create(DIR_OBJ, recursive = TRUE)
if (!dir.exists(DIR_FIG)) dir.create(DIR_FIG, recursive = TRUE)

# ==============================================================================
# Load Master Data
# ==============================================================================
all_obj <- readRDS(INPUT_FILE)

print(table(all_obj$orig.ident))

seu_query <- subset(all_obj, subset = orig.ident != "Human")

print(table(seu_query$orig.ident))

query_save_path <- file.path(DIR_OBJ, "All_Query_Species_Raw.rds")
saveRDS(seu_query, query_save_path)

rm(seu_query)
gc()