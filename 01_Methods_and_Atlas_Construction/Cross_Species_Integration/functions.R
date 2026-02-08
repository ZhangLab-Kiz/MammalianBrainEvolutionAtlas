#

# rename genes in Seurat object
reNameGene = function(obj, newname){
  rownames(obj@assays$RNA@features@.Data) = newname
  obj@assays$RNA@layers$counts@Dimnames[[1]] = newname
  obj
}

# transfer h5ad to Seurat object
toSeuratObj = function(anndata, gene, spe){
  # select the genes in RBH
  mouse_sc <- CreateSeuratObject(counts = t(anndata$X[,gene]), project = spe, meta.data = anndata$obs,
                                 min.cells = 30, min.features = 500)
  
  mouse_sc = NormalizeData(mouse_sc)
  mouse_sc = FindVariableFeatures(mouse_sc)
  mouse_sc = ScaleData(mouse_sc)
  mouse_sc = RunPCA(mouse_sc, features = VariableFeatures(mouse_sc))
  mouse_sc = RunUMAP(mouse_sc, dims = 1:20)
  # Scrublet
  scrublet = reticulate::import("scrublet")
  scr <- scrublet$Scrublet(counts_matrix = t(mouse_sc@assays$RNA@layers$counts),
                           expected_doublet_rate = 0.12)
  scr_out = scr$scrub_doublets()
  
  mouse_sc$doulet_scores = scr_out[[1]]
  mouse_sc$predicted_doulets = mouse_sc$doulet_scores > 0.25
  
  return(mouse_sc)
}

# plot cell percent 
plot_cell_pct = function(obj, x, fill, cols = NA){
  obj = obj@meta.data
  if(is.na(cols)[1]) {cols = scales::hue_pal()(length(unique(obj[[fill]])))}
  ggplot(obj, aes(get(x), y = 1, fill = get(fill))) +
    geom_bar(position = "fill", stat = "identity") + 
    cowplot::theme_cowplot() + xlab(x) + ylab("Percentage") +
    scale_fill_manual(values = cols)
}

plot_cell_number = function(obj, x, fill, cols = NA){
  obj = obj@meta.data
  if(is.na(cols)[1]) {cols = scales::hue_pal()(length(unique(obj[[fill]])))}
  ggplot(obj, aes(get(x), fill = get(fill))) +
    geom_bar() + 
    cowplot::theme_cowplot() + xlab(x) + ylab("# cells") +
    scale_fill_manual(values = cols)
}

