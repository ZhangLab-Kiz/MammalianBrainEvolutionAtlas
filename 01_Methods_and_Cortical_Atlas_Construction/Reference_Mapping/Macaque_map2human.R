#####Map macaque to human atlas
library(Seurat)
library(SeuratObject)
library(anndata)
library(reticulate)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(harmony)
use_python("/bin/python")


Macaque_sc = read_h5ad("./Macaque.h5ad")
Macaque_sc <- CreateSeuratObject(counts = t(Macaque_sc$X), project = "Macaque", meta.data = Macaque_sc$obs)
Macaque_out = Macaque_sc
Macaque_sc$species = "Macaque"
#QC
MT <- c("COX1","COX2","COX3","ND1","ND2","ND3","ND4L","ND4","ND5","ND6","CYTB","ATP6","ATP8")
Macaque_sc[["percent.mt"]] = colSums(x = GetAssayData(object = Macaque_sc, assay = "RNA", layer = "counts")[MT, , drop = FALSE]) / Macaque_sc[["nCount_RNA"]] * 100

Macaque_sc <- subset(Macaque_sc,subset = percent.mt <25)
Macaque_sc <- subset(Macaque_sc,subset = nCount_RNA >500)
scrublet = reticulate::import("scrublet")
scr <- scrublet$Scrublet(counts_matrix = t(Macaque_sc@assays$RNA@layers$counts),
                         expected_doublet_rate = 0.12)
scr_out = scr$scrub_doublets()

Macaque_sc$doulet_scores = scr_out[[1]]
Macaque_sc$predicted_doulets = Macaque_sc$doulet_scores > 0.25
Macaque_sc = subset(Macaque_sc, predicted_doulets == FALSE)
#ref human
Human_sobj = readRDS("Human_sobj_refenence.RDS")
overlap.gene = intersect(rownames(Human_sobj), rownames(Macaque_sc))
#basic analysis
Macaque_sc = NormalizeData(Macaque_sc[overlap.gene, ])
Macaque_sc = FindVariableFeatures(Macaque_sc)
Macaque_sc = ScaleData(Macaque_sc)
Macaque_sc = RunPCA(Macaque_sc)
Macaque_sc = RunHarmony(Macaque_sc,reduction = "pca",group.by.vars = "ID",reduction.save = "harmony")

#project to Human umap
human.anchors <- FindTransferAnchors(reference = Human_sobj, query = Macaque_sc, dims = 1:20, # dims = 1:50,
                                     reduction = "rpca",
                                     reference.reduction = "pca")

predictions <- TransferData(anchorset = human.anchors, refdata = Human_sobj$supercluster_term, dims = 1:20)
Macaque_sc <- AddMetaData(Macaque_sc, metadata = predictions)
Human_sobj$supercluster_term = factor(Human_sobj$supercluster_term)
Macaque_sc$predicted.id = factor(Macaque_sc$predicted.id, levels = levels((Human_sobj$supercluster_term)))
Macaque_sc <- IntegrateEmbeddings(anchorset = human.anchors, reference = Human_sobj, query = Macaque_sc,
                              new.reduction.name = "ref.pca")
Macaque_sc <- ProjectUMAP(Macaque_sc, query.reduction = "ref.pca", reference = Human_sobj, 
                      reference.reduction = "pca", reduction.model = "umap") 

human_p = FetchData(Human_sobj, vars = c("umap_1", "umap_2", "supercluster_term", "species"))
Macaque_p = FetchData(Macaque_sc, vars = c("refUMAP_1", "refUMAP_2", "predicted.id", "species"))
colnames(Macaque_p) = colnames(human_p)
overlap_p = rbind(human_p, Macaque_p)
ggplot(overlap_p, aes(umap_1, umap_2, color = species)) +
  geom_point(size=0.1, shape=20, alpha = 0.1) +
  theme_void()

non_neuron = c("Astrocyte", "Oligodendrocyte", "Microglia", "Oligodendrocyte precursor", "Committed oligodendrocyte precursor",
               "Vascular", "Ependymal", "Fibroblast", "Choroid plexus", "Bergmann glia")
Macaque_sc$SuperCellType =  ifelse(Macaque_sc$predicted.id %in% non_neuron, 
                                   "NonNeuron", "Neuron")

Macaque_sc$cell_type = as.character(Macaque_sc$predicted.id)
Macaque_sc$cell_type[!Macaque_sc$predicted.id %in% non_neuron] <- "neuron" 

#save out
Macaque_out = Macaque_out[,colnames(Macaque_sc)]
Macaque_out@meta.data = Macaque_sc@meta.data[rownames(Macaque_out@meta.data), ]

Macaque_out = NormalizeData(Macaque_out)
VariableFeatures(Macaque_out) = VariableFeatures(Macaque_sc)
Macaque_out = ScaleData(Macaque_out)
Macaque_out = RunPCA(Macaque_out)

Macaque_out = RunUMAP(Macaque_out, dims = 1:30, 
                      umap.method = "umap-learn", return.model = TRUE)


saveRDS(Macaque_sc, "Macaque_sc_ref.RDS")