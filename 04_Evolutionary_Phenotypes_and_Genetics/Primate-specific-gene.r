
library(Seurat)
library(dplyr)
library(tidyverse)


spe_vector <- species_vector <- c(
  "Suncus_murinus",
  "Hipposideros_armiger",
  "Felis_catus",
  "Meles_meles",
  "Mustela_putorius_furo",
  "Nyctereutes_procyonoides",
  "Vulpes_lagopus",
  "Canis_lupus_familiaris",
  "Bos_taurus",
  "Capra_hircus",
  "Sus_scrofa",
  "Mus_musculus",
  "Rattus_norvegicus",
  "Meriones_unguiculatus",
  "Phodopus_sungorus",
  "Cavia_porcellus",
  "Chinchilla_lanigera",
  "Heterocephalus_glaber",
  "Oryctolagus_cuniculus",
  "Lemur_catta",
  "Callithrix_jacchus",
  "Macaca_mulatta",
  "Homo Sapiens"
)

spe_df <- data.frame(Species = spe_vector, order = c(rep("Non_Primates"),19),rep("Primates",4))

###Seurat

rds <- readRDS("ct.rds") ## ct.rds is a process Seurat object

rds.meta <- rds@meta.data

result.meta <- left_join(rds.meta, spe_df, by = "Species")
head(result.meta)
rownames(result.meta) <- rownames(rds.meta)
rds@meta.data <- result.meta

Idents(rds) <- "order"
rds <- JoinLayers(rds)

Seurat_res <- FindMarkers(
  object = rds,
  ident.1 = "Primates",
  ident.2 = "Non_Primates",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

Seurat_res_sign <- subset(Seurat_res, p_val_adj < 0.05 & avg_log2FC > 1 & pct.2 < 0.2)

###Pseudo-bulk

gene_expr <- readRDS("gene_expr.rds") # gene_expr.rds is a matrix 

expr_mat = t(gene_expr)
expr_mat = expm1(expr_mat)[,spe_df$Species]

Pseudo_res = apply(expr_mat, 1, function(g){
    p = wilcox.test(g[c(20:23)], g[-c(20:23)], alternative = "greater", exact = T)$p.value
    fc = log2(mean(g[c(20:23)])/(mean(g[-c(20:23)]) + mean(g[c(20:23)]*0.1)))
    c(p, fc)
})
  
Pseudo_res = as.data.frame(t(Pseudo_res))
colnames(Pseudo_res) = c("p", "logFC")
Pseudo_res$gene = rownames(Pseudo_res)

Pseudo_res_sign <- subset(Pseudo_res, p < 0.01 & logFC > 1.5)
