##Other species map to macacque
###Each species was processed using the same workflow/pipeline.
library(dplyr)
source("/project/luyz/work/singlecell/00.scripts/functions.R")
library(Seurat)
library(SeuratObject)
library(anndata)
library(ggplot2)
library(reticulate)
library(harmony)
use_python("/project/luyz/miniconda3/envs/seurat/bin/python")


# find RBH between the species of interest and Macaque
find_RBH = function(spe){
  S1 = read.table(paste0("../00.Blastp/Ref_Macaque/Ref_Macaque_",spe,"/",spe,"_Macaque.BLASTP.txt"))
  S2 = read.table(paste0("../00.Blastp/Ref_Macaque/Ref_Macaque_",spe,"/Macaque_",spe,".BLASTP.txt"))
  
  # filter the low quality matches
  S1 = S1[S1$V11 < 1e-6, ]
  S2 = S2[S2$V11 < 1e-6, ]
  
  # find the best match
  S1 = S1[order(S1$V11), ]
  S2 = S2[order(S2$V11), ]
  
  S1 = S1 %>% distinct(V1, .keep_all = TRUE)
  S2 = S2 %>% distinct(V2, .keep_all = TRUE)
  
  Forword_best = unique(paste(S1$V1, S1$V2))
  Reverse_best = unique(paste(S2$V2, S2$V1))
  
  # RBH
  RBH = c(Forword_best, Reverse_best)
  RBH = data.frame(table(RBH))
  RBH = RBH[RBH$Freq == 2, ]
  
  df = data.frame(stringr::str_split_fixed(RBH$RBH, " ", 2))
  
  # convert ensembl ID to Gene SYBOL
  macaque_id = read.csv("/project/luyz/work/singlecell/27.Macaca/Macaque_gene_infor.csv")
  rownames(macaque_id) = macaque_id$gene_ids
  df$X2 = macaque_id[df$X2, 'X']
  return(na.omit(df))
}

# integrate
map2macaque = function(spe, h5ad_path, doublet_remove = TRUE, vars2regress = NA, output_dir = NA){
  
  if(is.na(output_dir)){
    output_dir = spe
  }
  
  if(!dir.exists(file.path(output_dir))){
    dir.create(file.path(output_dir))
  }
  
  output_plots = paste0(output_dir, "/plots")
  if(!dir.exists(file.path(output_plots))){
     dir.create(file.path(output_plots))
  }
   
  # find the RBH
  RBH_spe = find_RBH(spe)
  
  obj_sc = read_h5ad(h5ad_path)
  obj_sc <- CreateSeuratObject(counts = t(obj_sc$X), project = spe, meta.data = obj_sc$obs,
                               min.cells = 0, min.features = 0) 
  MT <- c("COX1","COX2","COX3","ND1","ND2","ND3","ND4L","ND4","ND5","ND6","CYTB","ATP6","ATP8")
  mt_path <- paste0("../00.Blastp/Ref_Macaque/Ref_Macaque_",spe,"/",spe, ".MT.list")
  if (file.exists(mt_path)) {
    mtlist <- readLines(mt_path)
    MT.all <- unique(c(MT,mtlist))
  } else {
    warning(paste("File", mt_path, "not found."))
    MT.all <- unique(c(MT))  
  }
  obj_sc[["percent.mt"]] = colSums(x = GetAssayData(object = obj_sc, assay = "RNA", layer = "counts")[intersect(MT.all, rownames(obj_sc)), , drop = FALSE]) / obj_sc[["nCount_RNA"]] * 100
  #VlnPlot(obj_sc,features = c("nFeature_RNA", "nCount_RNA","percent.mt"),pt.size=-1)
  RBH_spe$X1 = stringr::str_replace_all(RBH_spe$X1, "_", "-")
  obj_sc <- obj_sc[RBH_spe$X1,]
  # remove doublets by Scrublet
if(doublet_remove == TRUE){
    pre_filter_cell_count <- ncol(obj_sc)
    scrublet = reticulate::import("scrublet")
    scr <- scrublet$Scrublet(counts_matrix = t(obj_sc@assays$RNA@layers$counts),
                             expected_doublet_rate = 0.12)
    scr_out = scr$scrub_doublets()
    
    obj_sc$doulet_scores = scr_out[[1]]
    obj_sc$predicted_doulets = obj_sc$doulet_scores > 0.25
    
    doublets_plot = ggplot(obj_sc@meta.data, aes(x=doulet_scores)) + 
      geom_histogram(binwidth = .01, fill = "navy") +
      geom_vline(xintercept = 0.25, color = "orange") +
      theme_classic()
    
    obj_sc = subset(obj_sc, predicted_doulets == F) 
    post_filter_cell_count <- ncol(obj_sc)
    removed_cells <- pre_filter_cell_count - post_filter_cell_count
    cell_counts <- data.frame(
        Stage = c("Before doublet removal", "After doublet removal","Removed by doublets filter"),
        Cell_Count = c(pre_filter_cell_count, post_filter_cell_count,removed_cells)
    )
    
    write.csv(cell_counts, file = "removed_cells_doublets.csv", row.names = FALSE)
}
  # transfer genes to macaque genes
  spe2macaque = RBH_spe[match(rownames(obj_sc), RBH_spe$X1), "X2"]
  obj_sc = reNameGene(obj_sc, spe2macaque)


 #filter by mt count 
  pre_filter_cell_count <- ncol(obj_sc)
  obj_sc <- subset(obj_sc, subset = percent.mt < 25)
  post_filter_cell_count_mt <- ncol(obj_sc)
  removed_cells_mt <- pre_filter_cell_count - post_filter_cell_count_mt
  mt_cell_counts <- data.frame(
    Stage = c("Before mt filter", "After mt filter", "Removed by mt filter"),
    Cell_Count = c(pre_filter_cell_count, post_filter_cell_count_mt, removed_cells_mt)
  )
  write.csv(mt_cell_counts, file = "removed_cells_mt.csv", row.names = FALSE)

 #filter by nFeature count
  obj_sc <- subset(obj_sc, subset = nFeature_RNA > 500)
  post_filter_cell_count_nFeature <- ncol(obj_sc)
  removed_cells_nFeature <- post_filter_cell_count_mt - post_filter_cell_count_nFeature
  nFeature_cell_counts <- data.frame(
    Stage = c("Before nFeature_RNA filter", "After nFeature_RNA filter", "Removed by nFeature_RNA filter"),
    Cell_Count = c(post_filter_cell_count_mt, post_filter_cell_count_nFeature, removed_cells_nFeature)
  )
  write.csv(nFeature_cell_counts, file = "removed_cells_nFeature.csv", row.names = FALSE)

  obj_sc = NormalizeData(obj_sc)
  obj_sc = FindVariableFeatures(obj_sc)
  if(is.na(vars2regress[1])){
    obj_sc = ScaleData(obj_sc)
  } else {
    obj_sc = ScaleData(obj_sc, vars.to.regress = vars2regress)
  }
  
  obj_sc = RunPCA(obj_sc)
  #run harmony
  obj_sc <- RunHarmony(obj_sc,reduction = "pca",group.by.vars = "ID",reduction.save = "harmony",ncores = 20,max_iter = 10)
  #pca harmony plots
  p = DimPlot(obj_sc, group.by = "ID", reduction = "pca", raster=FALSE, pt.size = 0.1) 
  ggsave( paste0(output_plots, "/ID.pca.png"), plot = p,
          dpi = 600, width = 15, height = 8)
  p = DimPlot(obj_sc, group.by = "ID", reduction = "harmony", raster=FALSE, pt.size = 0.1) 
  ggsave( paste0(output_plots, "/ID.harmony.png"), plot = p,
          dpi = 600, width = 15, height = 8)
  obj_sc[["pca"]] = obj_sc[["harmony"]]
  obj_sc[["harmony"]] = NULL

  #  project to macaque UMAP
  Macaque_sc = readRDS("/project/luyz/work/singlecell/27.Macaca/Macaca_result/Macaque_sc_ref.RDS")
  Idents(Macaque_sc) = "predicted.id"
  
  macaque.anchors <- FindTransferAnchors(reference = Macaque_sc, query = obj_sc, dims = 1:30,
                                         reduction = "rpca",
                                       reference.reduction = "pca")
  predictions <- TransferData(anchorset = macaque.anchors, refdata = Macaque_sc$predicted.id, dims = 1:30)
  
  obj_sc <- AddMetaData(obj_sc, metadata = predictions)
  obj_sc$predicted.id = factor(obj_sc$predicted.id, levels = levels((Macaque_sc$predicted.id)))
  obj_sc$predict_confidence = ifelse(predictions$prediction.score.max < 0.5, "low", "high")
  
  Idents(obj_sc) = "predicted.id"
  obj_sc <- IntegrateEmbeddings(anchorset = macaque.anchors, reference = Macaque_sc, query = obj_sc,
                                new.reduction.name = "ref.pca")
  
  obj_sc <- ProjectUMAP(obj_sc, query.reduction = "ref.pca", reference = Macaque_sc,
                        umap.method = "umap-learn",
                        reference.reduction = "pca", reduction.model = "umap") 
  
  # filter low predictive cells
  obj_sc = subset(obj_sc, prediction.score.max > 0.5)
  
  # add cell types in meta.data
  non_neuron = c("Astrocyte", "Oligodendrocyte", "Microglia", "Oligodendrocyte precursor", "Vascular",
                 "Ependymal", "Fibroblast", "Choroid plexus", "Bergmann glia","Committed oligodendrocyte precursor")

  obj_sc$super_cell_type = ifelse(obj_sc$predicted.id %in% non_neuron, 
                                    "NonNeuron", "Neuron")
  
  obj_sc$cell_type = as.character(obj_sc$predicted.id)
  obj_sc$cell_type[!obj_sc$predicted.id %in% non_neuron] <- "Neuron" 
  
  
  # plots
  Macaque_p = FetchData(Macaque_sc, vars = c("umap_1", "umap_2", "predicted.id", "orig.ident"))
  obj_p = FetchData(obj_sc, vars = c("refUMAP_1", "refUMAP_2", "predicted.id", "orig.ident"))
  colnames(obj_p) = colnames(Macaque_p)
  overlap_p = rbind(Macaque_p, obj_p)
  p = ggplot(overlap_p, aes(umap_1, umap_2, color = orig.ident)) +
    geom_point(size=0.3, shape=20, alpha = 0.1) +
    theme_void()
  ggsave( paste0(output_plots, "/integrated.umap.png"), plot = p,
          dpi = 600, width = 6, height = 5)
  
  p = DimPlot(obj_sc, group.by = "predicted.id", reduction = "ref.umap", raster=FALSE, pt.size = 0.1) &
    scale_color_manual(values = as.character(c(pals::polychrome(30))))
  ggsave( paste0(output_plots, "/predicted.id.umap.png"), plot = p,
          dpi = 600, width = 15, height = 8)
  p1 = DimPlot(obj_sc, group.by = "ID", reduction = "ref.umap", raster=FALSE, pt.size = 0.1) 
  p2 = DimPlot(obj_sc, group.by = "Organ", reduction = "ref.umap", raster=FALSE, pt.size = 0.1) 
  p = p1+p2
  ggsave(paste0(output_plots,"/Organ.umap.png"), plot = p,
          dpi = 600, width = 15, height = 8)
  #cell pct plot cell number plot
  p = plot_cell_pct(obj_sc, "predicted.id", "ID",  as.character(c(pals::polychrome(30), pals::kelly(20)))) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
  ggsave( paste0(output_plots, "/sample.ID.cellPercent.umap.pdf"), plot = p,dpi = 600, width = 10, height = 5)
  ggsave( paste0(output_plots, "/sample.ID.cellPercent.umap.png"), plot = p,dpi = 600, width = 10, height = 5)
  p = plot_cell_number(obj_sc, "predicted.id", "ID",  as.character(c(pals::polychrome(30), pals::kelly(20)))) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
  ggsave( paste0(output_plots, "/sample.ID.cellNumber.umap.pdf"), plot = p,dpi = 600, width = 10, height = 5)
  ggsave( paste0(output_plots, "/sample.ID.cellNumber.umap.png"), plot = p,dpi = 600, width = 10, height = 5)
  #
  p = ggplot(obj_sc@meta.data, aes(x=prediction.score.max)) + 
    geom_histogram(binwidth = .01, fill = "navy") + theme_classic()
  ggsave( paste0(output_plots, "/predicted.maxscore.distribution.pdf"), plot = p,
          width = 4, height = 4)
  
  ggsave( paste0(output_plots, "/doublets.score.distribution.pdf"), plot = doublets_plot,
          width = 4, height = 4)
  obj_sc$organ = obj_sc$Organ
  all_pro <- obj_sc@meta.data %>%
  group_by(orig.ident,organ, batch, predicted.id) %>%
  tally(name = "count") %>%
  group_by(organ, batch) %>%
  mutate(proportion = count / sum(count)) 

  p = ggplot(all_pro, aes(x = batch, y = proportion, fill = predicted.id)) +
  geom_bar(stat = "identity", position = "stack") +  
  facet_grid(organ ~ orig.ident) +  
  labs(x = "Batch", y = "Proportion of Predicted ID", title = "Proportion of Predicted ID by Batch, Organ, and Orig.ident") +
  theme_minimal() +
  scale_color_manual(values = pals::alphabet2()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(paste0(output_plots,"/batch-organ-cell-pro.png"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_plots,"/batch-organ-cell-pro.pdf"), plot = p, width = 8, height = 6, dpi = 300)

  # save file
  saveRDS(obj_sc, paste0(output_dir, "/AllCells.RDS"))
  saveRDS(subset(obj_sc, super_cell_type == "Neuron"), paste0(output_dir, "/Neuron.RDS"))
  saveRDS(subset(obj_sc, super_cell_type == "NonNeuron"), paste0(output_dir, "/NonNeuron.RDS"))
  
  write.csv(RBH_spe, paste0(output_dir, "/RBH.gene.csv"), quote = F, row.names = F)
  write.csv(obj_sc@meta.data, paste0(output_dir, "/meta_data.csv"), quote = F)
  
  write.csv(AverageExpression(Macaque_sc, assays = "RNA")$RNA, paste0(output_dir, "/avgExp.Macaque.cellType.csv"), quote = F)
  write.csv(AverageExpression(obj_sc, assays = "RNA")$RNA, paste0(output_dir, "/avgExp.", spe,".cellType.csv"), quote = F)
  
  
  return(obj_sc)
}


## Rscript
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3){
  stop("Usage: Rscript Cell_Map_Macaque.R species h5ad_path output_dir")
} 
map2macaque(spe = args[1], h5ad_path = args[2], output_dir = args[3])


