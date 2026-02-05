# Oligo Vs Brain Weight
spe = read.csv("sub_spe_caperdata.csv")
df <- spe[,c("Species", "BW")]

BW <- df %>%
  pivot_wider(names_from = Species,  values_from = BW)

bw_values <- as.numeric(BW)
names(bw_values) <- colnames(BW)



plot_expr_BW_cell = function(gene, cells){
  expr_cell = lapply(cells, function(cell){
    expr = avg_expr_list[[cell]][, gene]
    expr = expr / max(expr)
    expr_df = as.data.frame(cbind(expr[names(bw_values)], log10(bw_values+1)))
    colnames(expr_df) = c("expr", "BW")
    expr_df$cell = cell
    expr_df$species = rownames(expr_df)
    expr_df
  })
  expr_cell = do.call(rbind, expr_cell)
  
  r_df <- expr_cell %>%
    group_by(cell) %>%
    summarise(r = cor(expr, BW)) %>%
    mutate(lab = sprintf("r = %.2f", r))
  expr_cell <- expr_cell %>%
    left_join(r_df, by = "cell")
  
  ggplot(expr_cell,
         aes(BW, expr)) +
    
    geom_smooth(
      aes(colour = cell),   
      method = "lm", se = TRUE, fill = "grey90"
    ) +
    geom_point(aes(colour = cell)) +
    scale_linetype_manual(
      name = "Pearson r",       
      values = rep("solid", nrow(r_df))  
    ) +
    guides(
      linetype = guide_legend(order = 2, override.aes = list(size = 1)),  
      colour   = guide_legend(order = 1)                                
    ) +
    scale_color_manual(values = ct_col_new) +
    cowplot::theme_cowplot()+
    theme(
      legend.position      = c(0.98, 0.02),  
      legend.justification = c(1, 0),         
      legend.background    = element_rect(fill = alpha("white", 0.6), colour = NA))
}

#	MBP, PLP, CNP, MAG
plot_expr_BW_cell("MBP", c("Oligodendrocyte"))


expr_cor_BW = function(){
  x = avg_expr_list[[3]][names(bw_values), ]
  res = apply(x, 2, function(y){
    pearson_test <- cor.test((as.numeric(y)), log(bw_values+1), method = "pearson")
    spearman_test <- cor.test((as.numeric(y)), log(bw_values+1), method = "spearman")
    mean_expr = mean(as.numeric(y))
    max_expr = max(as.numeric(y))
    pearson_cor = as.numeric(pearson_test$estimate)
    pearson_p = pearson_test$p.value
    spearman_cor = as.numeric(spearman_test$estimate)
    spearman_p = spearman_test$p.value
    c(mean_expr, max_expr, pearson_cor, pearson_p, spearman_cor, spearman_p)
  })
  res = as.data.frame(t(res))
  colnames(res) = c("mean_expr", "max_expr", "pearson_r", "pearson_p", "spearman_r", "spearman_p")
  res$pearson_padj = p.adjust(res$pearson_p, method = "fdr")
  res$spearman_padj = p.adjust(res$spearman_p, method = "fdr")
  res
}
expr_cor_BW_gene = expr_cor_BW()

plot_expr_BW_cell("LRIG1", c("Oligodendrocyte", "Upper-layer excitatory", "Deep-layer excitatory", "MGE interneuron", "CGE interneuron", "Astrocyte", "Microglia"))
plot_expr_BW_cell("LRIG1", c("Oligodendrocyte"))

expr_cor_BW_gene_p = na.omit(expr_cor_BW_gene[expr_cor_BW_gene$pearson_padj < 0.1 &  (expr_cor_BW_gene$spearman_r) > 0.3 & (expr_cor_BW_gene$pearson_r) > 0.3 & expr_cor_BW_gene$max_expr > 1, ])


Oligo_expr = avg_expr_list[["Oligodendrocyte"]]
Oligo_expr = Oligo_expr[, rownames(expr_cor_BW_gene_p)]
Oligo_expr = Oligo_expr[names(bw_values[order(bw_values)]), ]
spe_ann = data.frame(BrainWeight = log10(df$BW))
rownames(spe_ann) = df$Species

pdf("../plots/Oligo_cor_BrainWeight_genes.heatmap.pdf", width = 6, height = 12)
pheatmap::pheatmap(t(log(Oligo_expr[,sort(colnames(Oligo_expr))]+1)), 
                   color = colorRampPalette(c("lightblue2","white","red2"))(100),
                   border_color = NA,
                   annotation_col = spe_ann,
                   scale = "row", cluster_cols = F, cluster_rows = F)
dev.off()

# GO enrichment
library(enrichR)
available_dbs <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "KEGG_2016", "KEGG_2019_Human", "Disease_Perturbations_from_GEO_up", 
         "WikiPathways_2024", "WikiPathways_2024_Human", "Panther_2016", "SynGO_2024", "Reactome_Pathways_2024",
         "GWAS_Catalog_2025","UK_Biobank_GWAS_v1", "ClinVar_2025", "OMIM_Disease", "OMIM_Expanded")

enriched_results <- enrichr(rownames(expr_cor_BW_gene_p), dbs)
plotEnrich(enriched_results[[11]], showTerms = 20, numChar = 80, 
           y = "Count", orderBy = "P.value")

Oligo_BW_GO.plot = as.data.frame(enriched_results[[1]])
ggplot(Oligo_BW_GO.plot[c(1:4, 6:10, 14),], aes(fct_reorder(Term, -log10(P.value)), -log10(P.value)))+
  geom_bar(stat = "identity", fill = "red3", width = 0.8) + coord_flip() + xlab("")+
  cowplot::theme_cowplot()
ggsave("../plots/Oligo_cor_BrainWeight_genes.GO.pdf", width = 8, height = 4)


## validate expression changes by Nature 2021: https://www.nature.com/articles/s41586-021-03465-8
Nature2021_Glia = readRDS("../data/Nature2021/sample.combined_glia_integration.RDS")
Nature2021_Glia = UpdateSeuratObject(Nature2021_Glia)
DimPlot(Nature2021_Glia)

DefaultAssay(Nature2021_Glia) = "RNA"
Nature2021_Glia = NormalizeData(Nature2021_Glia)

Nature2021_oligo = AverageExpression(subset(Nature2021_Glia, subclass_label == "Oligo"), group.by = "orig.ident", layer = "data")$RNA
Nature2021_oligo = log(Nature2021_oligo[rownames(Nature2021_oligo) %in% colnames(Oligo_expr), c(3,2,1)]+1)
pdf("../plots/Oligo_cor_BrainWeight_Nature2021.heatmap.pdf", width = 3, height = 10)
pheatmap::pheatmap(Nature2021_oligo[sort(rownames(Nature2021_oligo)), ], 
                   color = colorRampPalette(c("lightblue2","white","red2"))(100),
                   border_color = NA,
                   scale = "row", cluster_cols = F, cluster_rows = F)
dev.off()


### Micorglia
expr_cor_BW_MG = function(){
  x = avg_expr_list[[9]][names(bw_values), ]
  res = apply(x, 2, function(y){
    pearson_test <- cor.test((as.numeric(y)), log10(bw_values+1), method = "pearson")
    spearman_test <- cor.test((as.numeric(y)), log10(bw_values+1), method = "spearman")
    mean_expr = mean(as.numeric(y))
    max_expr = max(as.numeric(y))
    pearson_cor = as.numeric(pearson_test$estimate)
    pearson_p = pearson_test$p.value
    spearman_cor = as.numeric(spearman_test$estimate)
    spearman_p = spearman_test$p.value
    c(mean_expr, max_expr, pearson_cor, pearson_p, spearman_cor, spearman_p)
  })
  res = as.data.frame(t(res))
  colnames(res) = c("mean_expr", "max_expr", "pearson_r", "pearson_p", "spearman_r", "spearman_p")
  res$pearson_padj = p.adjust(res$pearson_p, method = "fdr")
  res$spearman_padj = p.adjust(res$spearman_p, method = "fdr")
  res
}
expr_cor_BW_MG_gene = expr_cor_BW_MG()

plot_expr_BW_cell("SLC4A7", c("Microglia"))

expr_cor_BW_MG_gene_p = na.omit(expr_cor_BW_MG_gene[expr_cor_BW_MG_gene$pearson_padj < 0.1 & (expr_cor_BW_MG_gene$pearson_r) > 0.5 & expr_cor_BW_MG_gene$max_expr > 1, ])


MG_expr = avg_expr_list[["Microglia"]]
MG_expr = MG_expr[, rownames(expr_cor_BW_MG_gene_p)]
MG_expr = MG_expr[names(bw_values[order(bw_values)]), ]


spe_ann = data.frame(BrainWeight = log10(df$BW))
rownames(spe_ann) = df$Species

pdf("../plots/MG_cor_BrainWeight_genes.heatmap.pdf", width = 6, height = 5)
pheatmap::pheatmap(t(log(MG_expr[, sort(colnames(MG_expr))]+1)), 
                   color = colorRampPalette(c("lightblue2","white","red2"))(100),
                   border_color = NA,
                   annotation_col = spe_ann,
                   scale = "row", cluster_cols = F, cluster_rows = F)
dev.off()


dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "KEGG_2016", "KEGG_2019_Human", "Disease_Perturbations_from_GEO_up", 
         "WikiPathways_2024", "WikiPathways_2024_Human", "Panther_2016", "SynGO_2024", "Reactome_Pathways_2024",
         "GWAS_Catalog_2025","UK_Biobank_GWAS_v1", "ClinVar_2025", "OMIM_Disease", "OMIM_Expanded")

enriched_results <- enrichr(rownames(expr_cor_BW_MG_gene_p), dbs)
plotEnrich(enriched_results[[1]], showTerms = 20, numChar = 80, 
           y = "Count", orderBy = "P.value")

## validate by nature 2021 data
Nature2021_Glia = readRDS("../data/Nature2021/sample.combined_glia_integration.RDS")
Nature2021_Glia = UpdateSeuratObject(Nature2021_Glia)
DimPlot(Nature2021_Glia)

DefaultAssay(Nature2021_Glia) = "RNA"
Nature2021_Glia = NormalizeData(Nature2021_Glia)

Nature2021_MG = AverageExpression(subset(Nature2021_Glia, subclass_label == "Micro-PVM"), group.by = "orig.ident", layer = "data")$RNA
Nature2021_MG = log(Nature2021_MG[rownames(Nature2021_MG) %in% colnames(MG_expr), c(3,2,1)]+1)
pdf("../plots/MG_cor_BrainWeight_Nature2021.heatmap.pdf", width = 3, height = 5)
pheatmap::pheatmap(Nature2021_MG[sort(rownames(Nature2021_MG)), ], 
                   color = colorRampPalette(c("lightblue2","white","red2"))(100),
                   border_color = NA,
                   scale = "row", cluster_cols = F, cluster_rows = F)
dev.off()


