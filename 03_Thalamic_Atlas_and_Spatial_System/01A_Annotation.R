setwd("/data/xiaotx/Th")

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(data.table)
library(Cairo)
library(httpgd)
library(ggridges)

# ==============================================================================
# Load Data
# ==============================================================================

th <- readRDS("/data/xiaotx/Th/Object/Th.rds")

set.seed(1111)

# ==============================================================================
# Species Mapping
# ==============================================================================

species_map <- c(
  "Cangshu" = "Phodopus sungorus",
  "Cat" = "Felis catus",
  "Cav" = "Cavia porcellus",
  "Chin" = "Chinchilla lanigera",
  "Colobus" = "Colobus guereza",
  "Deer" = "Cervus nippon",
  "Dog" = "Canis lupus familiaris",
  "Glis" = "Graphiurus kelleni",
  "Goat" = "Capra hircus",
  "Heter" = "Heterocephalus glaber",
  "Hippo" = "Hipposideros armiger",
  "Hys" = "Hystrix brachyura",
  "Lemur" = "Lemur catta",
  "Macaque" = "Macaca mulatta",
  "Marmoset" = "Callithrix jacchus",
  "Marmot" = "Marmota himalayana",
  "Meles" = "Meles meles",
  "Mouse" = "Mus musculus",
  "Pag" = "Paguma larvata",
  "Pig" = "Sus scrofa",
  "Pro" = "Procyon lotor",
  "Rabbit" = "Oryctolagus cuniculus",
  "Rat" = "Rattus norvegicus",
  "Rhin" = "Rhinolophus blythi",
  "Rhiz" = "Rhizomys sinensis",
  "Rous" = "Rousettus leschenaultii",
  "Scap" = "Scaptonyx fusicaudus",
  "Tam" = "Tamias sibiricus",
  "Ferret" = "Mustela putorius furo",
  "Cattle" = "Bos taurus",
  "Ger" = "Meriones unguiculatus",
  "Cal" = "Callosciurus erythraeus",
  "Rdog" = "Nyctereutes procyonoides",
  "Afox" = "Vulpes lagopus",
  "Shrew" = "Suncus murinus",
  "Eri" = "Atelerix albiventris",
  "Tup" = "Tupaia chinensis",
  "Yak" = "Bos grunniens",
  "Pbv" = "Petaurus breviceps papuanus"
)

th@meta.data$latin_name <- species_map[as.character(th$orig.ident)]

# ==============================================================================
# Initial Filtering
# ==============================================================================

non_thalamic_types <- c(
  "Bergmann glia", "Cerebellar inhibitory", "Granule",
  "Deep-layer excitatory", "Upper-layer excitatory",
  "Medium spiny neuron", "Eccentric medium spiny neuron",
  "Hippocampal CA1-3", "Hippocampal CA4", "Hippocampal dentate gyrus"
  # "MGE interneuron", "CGE interneuron"
)

th_only <- subset(th, subset = predicted.id %in% non_thalamic_types, invert = TRUE)

th_filtered <- subset(th_only, subset = prediction.score.max > 0.6)

# ==============================================================================
# Dimensionality Reduction and Clustering
# ==============================================================================

th_filtered <- RunUMAP(th_filtered, reduction = "ref.pca",
                       dims = 1:30, verbose = TRUE)

th_filtered <- FindNeighbors(th_filtered, reduction = "ref.pca", dims = 1:30)

th_filtered <- FindClusters(th_filtered, resolution = 0.8,
                            random.seed = 1111, algorithm = 4, verbose = TRUE)

th_filtered <- subset(th_filtered, subset = ! (seurat_clusters %in% c("28")))

th_filtered <- RunUMAP(th_filtered, reduction = "ref.pca", dims = 1:30,
                       min.dist = 0.25, n.neighbors = 30, verbose = TRUE)

th_filtered <- FindNeighbors(th_filtered, reduction = "ref.pca",
                             dims = 1:30, verbose = TRUE)

th_filtered <- FindClusters(th_filtered, resolution = 0.8, algorithm = 4)

# ==============================================================================
# Cluster Annotation
# ==============================================================================

annotation_map <- c(
  "1" = "Oligodendrocyte",
  "2" = "Oligodendrocyte",
  "3" = "Astrocyte",
  "4" = "Thalamic excitatory",
  "5" = "Oligodendrocyte precursor",
  "6" = "Thalamic excitatory",
  "7" = "Thalamic excitatory",
  "8" = "Oligodendrocyte",
  "9" = "Thalamic excitatory",
  "10" = "Splatter",
  "11" = "Microglia",
  "12" = "Thalamic excitatory",
  "13" = "Microglia",
  "14" = "Thalamic excitatory",
  "15" = "Astrocyte",
  "16" = "Midbrain-derived inhibitory",
  "17" = "Thalamic excitatory",
  "18" = "Splatter",
  "19" = "Oligodendrocyte",
  "20" = "Thalamic excitatory",
  "21" = "Ependymal",
  "22" = "Splatter",
  "23" = "Vascular",
  "24" = "Fibroblast",
  "25" = "Astrocyte",
  "26" = "Thalamic excitatory",
  "27" = "Thalamic excitatory",
  "28" = "Thalamic excitatory",
  "29" = "CGE/MGE inhibitory",
  "30" = "Thalamic excitatory",
  "31" = "Oligodendrocyte",
  "32" = "Thalamic excitatory",
  "33" = "Astrocyte",
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
  "53" = "Oligodendrocyte",
  "54" = "Oligodendrocyte",
  "55" = "Oligodendrocyte"
)

new_metadata <- data.frame(
  manual.celltype = annotation_map[as.character(th_filtered$seurat_clusters)],
  row.names = colnames(th_filtered)
)

th_filtered <- AddMetaData(th_filtered, metadata = new_metadata)

th_final <- th_filtered

# ==============================================================================
# Sub-clustering Inhibitory Neurons
# ==============================================================================

inhibitory_subset <- subset(th_filtered,
                            subset = manual.celltype == "CGE/MGE inhibitory")

inhibitory_subset <- RunUMAP(inhibitory_subset, reduction = "ref.pca", dims = 1:30, verbose = FALSE)

inhibitory_subset <- FindNeighbors(inhibitory_subset, reduction = "ref.pca", dims = 1:30, verbose = FALSE)
inhibitory_subset <- FindClusters(inhibitory_subset, resolution = 0.5, verbose = FALSE)

sub_annotation_map <- c(
  "0" = "CGE interneuron",
  "1" = "MGE interneuron",
  "2" = "CGE interneuron",
  "3" = "MGE interneuron",
  "4" = "MGE interneuron",
  "5" = "MGE interneuron"
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


th_final$final.celltype <- as.character(th_final$manual.celltype)

new_subtypes_with_names <- as.character(inhibitory_subset$inhibitory.subtype)
names(new_subtypes_with_names) <- colnames(inhibitory_subset)

cells_to_update <- colnames(inhibitory_subset)

th_final$final.celltype[cells_to_update] <- new_subtypes_with_names

# ==============================================================================
# Phenotype Data Integration
# ==============================================================================

phenotype_file <- "/data/xiaotx/Th/Data/Phenotype.csv"
phenotype_data <- read.csv(phenotype_file, header = TRUE)

cols_to_add <- c(
  "orig.ident", "Genus", "Family", "Order", "Class", "Weight", "BW", "BWW", 
  "EQ", "BSR", "SL", "GL", "VL", "DM", "soc", "qixi", "xixing", "feeding", 
  "hearing", "olfaction", "touch", "vision", "life_span", "length"
)
phenotype_data_subset <- phenotype_data[, cols_to_add]

original_meta <- th_final@meta.data

merged_meta <- left_join(original_meta, phenotype_data_subset, by = c("species" = "orig.ident"))

rownames(merged_meta) <- rownames(original_meta)

th_final@meta.data <- merged_meta

print(head(th_final@meta.data[, (which(colnames(th_final@meta.data) == "species")):ncol(th_final@meta.data)]))

th_final <- RunUMAP(th_final,
                    reduction = "ref.pca",
                    dims = 1:30,
                    reduction.name = "ref.umap",
                    return.model = TRUE,
                    verbose = TRUE,
                    min.dist = 0.25,
                    n.neighbors = 30)

saveRDS(th_final, file = "/data/xiaotx/Th/Object/Th_final.rds")