# re-analyze human brain atlas

library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(anndata)
library(reticulate)
use_python("/home/zhangc/miniconda3/bin/python")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
source("functions.R")

# set this option when analyzing large datasets
options(future.globals.maxSize = 1e+09)

##### load human data from https://www.science.org/doi/10.1126/science.add7046
Human_Neruon_sc = read_h5ad("../snRNAseq_data/Human/Human_neuron.h5ad") 

Human_Neruon_sobj <- CreateSeuratObject(counts = t(Human_Neruon_sc$X), project = "Human", meta.data = Human_Neruon_sc$obs,
                                 min.cells = 0, min.features = 500)

# make gene name unique
Human_Neruon_sobj = Human_Neruon_sobj[!duplicated(Human_Neruon_sc$var[rownames(Human_Neruon_sobj), "Gene"]), ]

Human_Neruon_sobj = reNameGene(Human_Neruon_sobj, Human_Neruon_sc$var[rownames(Human_Neruon_sobj), "Gene"])



#
Human_NN_sc = read_h5ad("../snRNAseq_data/Human/Human_Non_neuron.h5ad")

Human_NN_sobj <- CreateSeuratObject(counts = t(Human_NN_sc$X), project = "Human", meta.data = Human_NN_sc$obs,
                                        min.cells = 0, min.features = 500)
Human_NN_meta_data = Human_NN_sobj@meta.data

# make gene name unique
Human_NN_sobj = Human_NN_sobj[!duplicated(Human_NN_sc$var[rownames(Human_NN_sobj), "Gene"]), ]
Human_NN_sobj = reNameGene(Human_NN_sobj, Human_NN_sc$var[rownames(Human_NN_sobj), "Gene"])


rm(Human_NN_sc, Human_Neruon_sc)
gc()

### merge neuron and non-neuron cells
Human_sobj = merge(Human_NN_sobj, Human_Neruon_sobj)
saveRDS(Human_sobj, "Human_sobj_refenence.RDS")

rm(Human_NN_sobj, Human_Neruon_sobj)
gc()

Human_sobj = subset(Human_sobj, tissue %in% c("cerebellum", "cerebral cortex", "hippocampal formation",
                                              "hypothalamus", "midbrain", "thalamic complex"))

Human_sobj$species = "Human"

# load macaque data
Macaque_sc = read_h5ad("Macaque.merge.h5ad")

Macaque_sc <- CreateSeuratObject(counts = t(Macaque_sc$X), project = "Macaque", meta.data = Macaque_sc$obs,
                               min.cells = 30, min.features = 500)

Macaque_sc$species = "Macaque"

# use overlaped genes
overlap.gene = intersect(rownames(Human_sobj), rownames(Macaque_sc))

#
Macaque_sc = NormalizeData(Macaque_sc[overlap.gene, ])
Macaque_sc = FindVariableFeatures(Macaque_sc)
Macaque_sc = ScaleData(Macaque_sc)
Macaque_sc = RunPCA(Macaque_sc)

#
Human_sobj = NormalizeData(Human_sobj[overlap.gene, ])
Human_sobj = FindVariableFeatures(Human_sobj)
Human_sobj = ScaleData(Human_sobj)
Human_sobj = RunPCA(Human_sobj)
Human_sobj = RunUMAP(Human_sobj, dims = 1:50, umap.method = "umap-learn", return.model = TRUE)

DimPlot(Human_sobj, group.by = c("supercluster_term"),  label = T)
DimPlot(subset(Human_sobj, cell_type == "neuron"), group.by = c("supercluster_term"),  label = T)
FeaturePlot(Human_sobj, features = c("DRD1", "DRD2", "CUX2", "GAD1"))
saveRDS(Human_sobj, "Human_sobj_refenence.RDS")

Human_sobj = readRDS("Human_sobj_refenence.RDS")
ElbowPlot(Human_sobj, ndims = 50, reduction = "pca")
ElbowPlot(Macaque_sc, ndims = 50, reduction = "pca")

#  project macaca data to human 
human.anchors <- FindTransferAnchors(reference = Human_sobj, query = Macaque_sc, dims = 1:20, # dims = 1:50,
                                     reduction = "rpca",
                                     reference.reduction = "pca")

predictions <- TransferData(anchorset = human.anchors, refdata = Human_sobj$supercluster_term, dims = 1:20)

Macaque_sc <- AddMetaData(Macaque_sc, metadata = predictions)
Human_sobj$supercluster_term = factor(Human_sobj$supercluster_term)
Macaque_sc$predicted.id = factor(Macaque_sc$predicted.id, levels = levels((Human_sobj$supercluster_term)))
DimPlot(Macaque_sc, group.by = c("predicted.id"),label = T)


# integration and project cell types
Macaque_sc <- IntegrateEmbeddings(anchorset = human.anchors, reference = Human_sobj, query = Macaque_sc,
                              new.reduction.name = "ref.pca")

Macaque_sc <- ProjectUMAP(Macaque_sc, query.reduction = "ref.pca", reference = Human_sobj, 
                      reference.reduction = "pca", reduction.model = "umap") 


Macaque_sc$predict_confidence = ifelse(rowMax(as.matrix(predictions[,-1])) < 0.5, "low", "high")
DimPlot(Macaque_sc, group.by = c("predict_confidence"), reduction = "ref.umap",label = T)

plot(density(rowMax(as.matrix(predictions[,-1]))))



human_p = FetchData(Human_sobj, vars = c("umap_1", "umap_2", "supercluster_term", "species"))
Macaque_p = FetchData(Macaque_sc, vars = c("refUMAP_1", "refUMAP_2", "predicted.id", "species"))
colnames(Macaque_p) = colnames(human_p)
overlap_p = rbind(human_p, Macaque_p)
ggplot(overlap_p, aes(umap_1, umap_2, color = species)) +
  geom_point(size=0.1, shape=20, alpha = 0.1) +
  theme_void()


saveRDS(Macaque_sc, "Macaque_sc_CellTyping.RDS")
Macaque_sc = readRDS("Macaque_sc_CellTyping.RDS")

Macaque_sc = subset(Macaque_sc, predict_confidence == "high")     

DimPlot(Macaque_sc, group.by = c("predicted.id"), label = T)


# remove doublets
# Scrublet
scrublet = reticulate::import("scrublet")
scr <- scrublet$Scrublet(counts_matrix = t(Macaque_sc@assays$RNA@layers$counts),
                         expected_doublet_rate = 0.12)
scr_out = scr$scrub_doublets()

Macaque_sc$doulet_scores = scr_out[[1]]
Macaque_sc$predicted_doulets = Macaque_sc$doulet_scores > 0.25
DimPlot(Macaque_sc, group.by = c("predicted_doulets"),  reduction = "ref.umap")
plot(density(Macaque_sc$doulet_scores))

Macaque_sc = subset(Macaque_sc, predicted_doulets == FALSE)

Macaque_sc$cell_type = as.character(Macaque_sc$predicted.id)
Macaque_sc$cell_type[!Macaque_sc$predicted.id %in% non_neuron] <- "neuron" 
saveRDS(Macaque_sc, "Macaque_sc_CellTyping1.RDS")


DimPlot(Macaque_sc, reduction = "ref.umap", group.by = "cell_type", label = T)&
  scale_color_manual(values = as.character(c(pals::polychrome(30), pals::kelly(20))))

DimPlot(Macaque_sc, reduction = "ref.umap", group.by = "predicted.id", label = T)&
  scale_color_manual(values = as.character(c(pals::polychrome(30), pals::kelly(20))))

# 
Marker = c("INA", "SLC17A6", "SLC17A7", "SLC32A1", "PTPRC", "CLDN5", "ACTA2", "LUM", "PDGFRA", "SOX10", "PLP1", "AQP4", "FOXJ1", "TTR")
FeaturePlot(Macaque_sc, reduction = "ref.umap", features = Marker,
            slot = "data", cols = c("gray90", "red")) & NoAxes()


# refine Miscellaneous & Splatter 
DimPlot(subset(Macaque_sc, predicted.id == "Miscellaneous"), group.by = "Organ", label = T)&
  scale_color_manual(values = as.character(c(pals::polychrome(30), pals::kelly(20))))
FeaturePlot(subset(Macaque_sc, predicted.id == "Splatter"& prediction.score.max > 0.9), "prediction.score.max")
FeaturePlot(subset(Macaque_sc, predicted.id == "Miscellaneous"& prediction.score.max > 0.9), "prediction.score.max")

Macaque_sc = subset(Macaque_sc, predicted.id == "Splatter" & prediction.score.max < 0.8, invert = T)
Macaque_sc = subset(Macaque_sc, predicted.id == "Miscellaneous" & prediction.score.max < 0.8, invert = T)

DimPlot(Macaque_sc, reduction = "ref.umap", group.by = "predicted.id", split.by = "predicted.id", label = T, ncol=8)&
  scale_color_manual(values = as.character(c(pals::polychrome(30), pals::kelly(20))))

saveRDS(Macaque_sc, "Macaque_sc_CellTyping2.RDS")
