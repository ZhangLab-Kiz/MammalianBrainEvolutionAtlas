library(anndata)
library(ggplot2)
library(reticulate)
library(tidyr)
library(pals)
library(devtools)
library(pathview)
library(stringr)
library(dplyr)
library(purrr)

spe = read.csv("sub_spe_caperdata.csv")$Species
pheno <- read.csv("EQ_adj.csv")
df <- pheno[pheno$Species %in% spe, c("Species", "EQ", "BW", "Weight")]
rownames(df) = df$Species
df = df[,-1]
pheatmap(df[order(df$EQ),], scale = "column", cluster_rows = F)

df <- pheno[pheno$Species %in% spe, c("Species", "EQ")]

EQ <- df %>%
  pivot_wider(names_from = Species,  values_from = EQ)

eq_values <- as.numeric(EQ)
names(eq_values) <- colnames(EQ)


ggplot(pheno, aes(Species, EQ))+
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_data <- pheno %>%
  pivot_longer(
    cols = c(EQ, BW, Weight), 
    names_to = "Measurement",
    values_to = "Value"      
  )

spe_color = readRDS("species_color_new.rds")

plot_data <- pheno

plot_data$Species <- factor(plot_data$Species,
                            levels = plot_data$Species[order(plot_data$EQ)])

ggplot(pheno[pheno$Species %in% spe, c("Species", "EQ")], aes(x = fct_reorder(Species, EQ), y = EQ, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip()+
  cowplot::theme_cowplot() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = spe_color) 
ggsave("../plots/EQ_spe.barplot.pdf", width = 5, height = 5)

ggplot(plot_data, aes(x = fct_reorder(Species, -log10(Weight)), y = log10(Weight), fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip()+
  cowplot::theme_cowplot() +
  labs(y = "log10(Body Weight), gram") +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = spe_color) 
ggsave("../plots/BodyWeight_spe.barplot.pdf", width = 5, height = 5)

ggplot(plot_data, aes(x = Species, y = BW, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip()+
  cowplot::theme_cowplot() +
  labs(y = "Brain Weight, gram") +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = spe_color) 
ggsave("../plots/BrainWeight_spe.barplot.pdf", width = 5, height = 5)



avg_expr_list <- readRDS("Gene_avg_data_20spe_spename_coexp_expm1.RDS")
names(avg_expr_list)
rownames(avg_expr_list[[1]])
colnames(avg_expr_list[[1]])
head(avg_expr_list[[1]])
nrow(avg_expr_list[[1]])


expr_cor = lapply(avg_expr_list, function(x){
  x = x[names(eq_values), ]
  res = apply(x, 2, function(y){
    pearson_test <- cor.test(as.numeric(y), eq_values, method = "pearson")
    spearman_test <- cor.test(as.numeric(y), eq_values, method = "spearman")
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
})

expr_cor_ct = expr_cor
for (i in names(expr_cor_ct)){
  expr_cor_ct[[i]]$gene = rownames(expr_cor_ct[[i]])
  expr_cor_ct[[i]]$cell = i
}
expr_cor_ct = do.call(rbind, expr_cor_ct)

ggplot(expr_cor_ct, aes(x = pearson_r, color = cell)) +
  geom_density( alpha = 0.6) + 
  theme_minimal()

# filter
expr_cor_ct_p = na.omit(expr_cor_ct[expr_cor_ct$pearson_padj < 0.05 & abs(expr_cor_ct$pearson_r) > 0.5 & expr_cor_ct$max_expr > 0.5, ])

expr_cor_ct_p = expr_cor_ct_p[order(expr_cor_ct_p$pearson_padj), ]
expr_cor_ct_p = expr_cor_ct_p[!expr_cor_ct_p$cell %in% c("Midbrain-derived inhibitory", "Eccentric medium spiny neuron"), ]

write.csv(expr_cor_ct_p, "expr_cor_ct_p.csv", quote = F)


ggplot(expr_cor_ct_p, aes(factor(cell, levels = levels(expr_cor_ct_p_plot$cell)), fill = cell))+
  geom_bar() +
  scale_fill_manual(values = ct_col_new)+
  cowplot::theme_cowplot()
ggsave("../plots/EQ_p.number.barplot.pdf", width = 6, height = 1)



ID_gene = read.csv("ClinGene intellectual disability gene.csv")
ID_gene$gene = sub("HGNC.*", "", ID_gene$Gene)
intersect(ID_gene$gene, expr_cor_ct_p$gene)

plot_expr_eq = function(gene, cell){
  expr = avg_expr_list[[cell]][, gene]
  expr_df = as.data.frame(cbind(expr[names(eq_values)], eq_values))
  colnames(expr_df) = c("expr", "EQ")
  expr_df$species = rownames(expr_df)
  ggplot(expr_df, aes(expr, EQ, species)) +
    #geom_text_repel(aes(label = species), size = 5, box.padding = 0)+
    geom_smooth(method = "lm", se = TRUE, color = "red3")+
    geom_point() +
    cowplot::theme_cowplot(font_size = 30) +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste0("R = ", round(sqrt(..r.squared..), digits = 2), ",  P = ", round(..p.value.., digits = 3))),
      parse = F, size = 12,
      label.x = "left",
      label.y = "top"
    ) +
    ggtitle(paste(cell, gene, sep = " - "))
}

pdf("example_cor_EQ_gene.scatter.pdf", 6, 6)

plot_expr_eq("MOK", "Thalamic excitatory")
plot_expr_eq("CALU", "Hippocampal CA1-3")
plot_expr_eq("ROBO1", "Hippocampal CA4")
plot_expr_eq("TMEM87B", "Deep-layer excitatory")
plot_expr_eq("RBFOX1", "Upper-layer excitatory")
plot_expr_eq("MYRIP", "CGE interneuron")
plot_expr_eq("COG6", "MGE interneuron")
plot_expr_eq("DISC1", "Hippocampal dentate gyrus")
plot_expr_eq("INTU", "Hippocampal dentate gyrus")
plot_expr_eq("AKNAD1", "Medium spiny neuron")
plot_expr_eq("PPP3CB", "Medium spiny neuron")
plot_expr_eq("MAP1A", "Microglia")
plot_expr_eq("APCDD1", "Oligodendrocyte precursor")
plot_expr_eq("NTNG1", "Oligodendrocyte precursor")
plot_expr_eq("KCNH8", "Oligodendrocyte")

plot_expr_eq("KMT2C", "Upper-layer excitatory")
plot_expr_eq("COG6", "Upper-layer excitatory")
plot_expr_eq("SATB1", "Upper-layer excitatory")
plot_expr_eq("PAK1", "Upper-layer excitatory")
plot_expr_eq("PAK3", "Upper-layer excitatory")

dev.off()


plot_expr_eq_cell = function(gene, cells){
  expr_cell = lapply(cells, function(cell){
    expr = avg_expr_list[[cell]][, gene]
    expr_df = as.data.frame(cbind(expr[names(eq_values)], eq_values))
    colnames(expr_df) = c("expr", "EQ")
    expr_df$cell = cell
    expr_df$species = rownames(expr_df)
    expr_df
  })
  expr_cell = do.call(rbind, expr_cell)
  
  r_df <- expr_cell %>%
    group_by(cell) %>%
    summarise(r = cor(expr, EQ)) %>%
    mutate(lab = sprintf("r = %.2f", r))
  expr_cell <- expr_cell %>%
    left_join(r_df, by = "cell")
  
  ggplot(expr_cell,
         aes(expr, EQ)) +
   
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
pdf("example_cor_EQ_COG6.scatter1.pdf", 6, 6)
plot_expr_eq_cell("COG6", c("Upper-layer excitatory","Deep-layer excitatory","CGE interneuron","MGE interneuron",
                             "Hippocampal CA1-3", "Oligodendrocyte", "Hippocampal dentate gyrus"))
dev.off()
pdf("example_cor_EQ_KATNBL1.scatter1.pdf", 6, 6)
plot_expr_eq_cell("KATNBL1", c("Upper-layer excitatory","Deep-layer excitatory"))
dev.off()

plot_expr_eq_cell("SNX9", c("Upper-layer excitatory","Deep-layer excitatory","CGE interneuron","MGE interneuron",
                               "Hippocampal dentate gyrus", "Medium spiny neuron","Astrocyte"))


# GO
library(clusterProfiler)
data_GO <- compareCluster(
  gene~cell, 
  data=expr_cor_ct_p, 
  universe = colnames(avg_expr_list[[1]]),
  fun="enrichGO", 
  keyType = "SYMBOL",
  ont = "ALL",
  OrgDb='org.Hs.eg.db',
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

dotplot(data_GO, showCategory=5,font.size = 8)


library(enrichR)
available_dbs <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "KEGG_2016", "KEGG_2019_Human", "Disease_Perturbations_from_GEO_up", 
         "WikiPathways_2024", "WikiPathways_2024_Human", "Panther_2016", "SynGO_2024", "Reactome_Pathways_2024",
         "GWAS_Catalog_2025","UK_Biobank_GWAS_v1", "ClinVar_2025", "OMIM_Disease", "OMIM_Expanded")

enriched_results <- enrichr(unique(expr_cor_ct_p[expr_cor_ct_p$cell_class == "Neurons", "gene"]), dbs)
enriched_results <- enrichr(unique(expr_cor_ct_p[expr_cor_ct_p$cell_class == "Glia", "gene"]), dbs)

enriched_results <- enrichr(unique(expr_cor_ct_p$gene), dbs)
plotEnrich(enriched_results[[1]], showTerms = 20, numChar = 80, 
           y = "Count", orderBy = "P.value")

enriched_results_df = do.call(rbind, enriched_results)
select_terms = c("Axonogenesis (GO:0007409)","Axon guidance", "Ras signaling pathway", 
                 "Calcium signaling pathway","Neural Crest Cell Migration During Development WP4564",
                "MICROCEPHALY",
                 "INTELLECTUAL DISABILITY")

select_terms_gene = enriched_results_df[enriched_results_df$Term %in% select_terms, c("Term", "Genes")]

library(tidyverse)

tidy_gene_term <- select_terms_gene %>% 
  as_tibble() %>%                       
  select(Term, Genes) %>%               
  mutate(gene = str_split(Genes, ";")) %>% 
  unnest(gene) %>%                      
  select(Term, gene) %>%                
  distinct()  

# GWAS
df = read.table("gwas-intelligence.tsv", sep = "\t", header = T)
gene_column <- na.omit(df$MAPPED_GENE)
gene_list <- strsplit(gene_column, "\\s+-\\s+")
all_genes <- unlist(gene_list)
gene_set <- unique(all_genes)
gene_set <- gene_set[gene_set != ""]
gene_set = unique(gene_set)
gene_set = gene_set[gene_set %in% colnames(avg_expr_list[[1]])] # 1051

unique(intersect(gene_set, expr_cor_ct_p$gene))

celltype_summary <- expr_cor_ct_p %>%
  group_by(cell) %>%
  summarise(
    total_genes = n_distinct(gene),
    in_gene_set = sum(unique(gene) %in% gene_set),
  )
celltype_summary$expected_n = celltype_summary$total_genes / 14933 * 1051
celltype_summary$log2obs_exp = log2((celltype_summary$in_gene_set + 0.1) / celltype_summary$expected_n)

ggplot(celltype_summary, aes(cell, log2obs_exp))+
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

expr_cor_ct_p$cell_class = "Neurons"
expr_cor_ct_p[expr_cor_ct_p$cell %in% c("Microglia", "Vascular","Fibroblast"), "cell_class"] <- "Non-neural cells"
expr_cor_ct_p[expr_cor_ct_p$cell %in% c("Astrocyte","Oligodendrocyte","Oligodendrocyte precursor", "Bergmann glia", "Ependymal"), "cell_class"] <- "Glia"
expr_cor_ct_p$cell_class = factor(expr_cor_ct_p$cell_class, levels = c("Non-neural cells", "Glia", "Neurons"))

cellclass_summary <- expr_cor_ct_p %>%
  group_by(cell_class) %>%
  summarise(
    total_genes = n_distinct(gene),
    in_gene_set = sum(unique(gene) %in% gene_set),
  )
cellclass_summary$expected_n = cellclass_summary$total_genes / 14933 * 1051
cellclass_summary$log2obs_exp = log2((cellclass_summary$in_gene_set + 0.1) / cellclass_summary$expected_n)

ggplot(cellclass_summary, aes(cell_class, log2obs_exp))+
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# cortical folding

pheno_CF = pheno$Group
names(pheno_CF) = pheno$Species
CF_1 = names(pheno_CF[pheno_CF == 1])[-14]
CF_0 = names(pheno_CF[pheno_CF == 0])

DEG_CF = lapply(avg_expr_list, function(x){
  res = apply(x, 2, function(y){
    deg <- wilcox.test(y[CF_1], y[CF_0])$p.value
    logFC = log2(mean(y[CF_1] + 0.01) / mean(y[CF_0] + 0.01))
    mean_expr = mean(y)
    c(mean_expr, logFC, deg)
  })
  res = as.data.frame(t(res))
  colnames(res) = c("mean_expr", "logFC", "p")
  res$padj = p.adjust(res$p, method = "fdr")
  res
})

for (i in names(DEG_CF)){
  DEG_CF[[i]]$gene = rownames(DEG_CF[[i]])
  DEG_CF[[i]]$cell = i
}
DEG_CF = do.call(rbind, DEG_CF)

DEG_CF_p = na.omit(DEG_CF[DEG_CF$mean_expr > 0.1 & DEG_CF$padj < 0.1 & abs(DEG_CF$logFC) > log2(2), ])
DEG_CF_p = DEG_CF_p[DEG_CF_p$cell != "Midbrain-derived inhibitory", ]
write.csv(DEG_CF_p, "DEG_CF_p.csv", quote = F)


Up_EQ_cor_gene = expr_cor_ct_p[expr_cor_ct_p$cell == "Upper-layer excitatory" & expr_cor_ct_p$pearson_p > 0, "gene"]
Up_DEG_gene = DEG_CF_p[DEG_CF_p$cell == "Upper-layer excitatory" & DEG_CF_p$logFC > 0, "gene"]

library(ggvenn)
ggvenn(list(A = Up_EQ_cor_gene, B = Up_DEG_gene),
       fill_color = c("#66C2A5", "#FC8D62"), fill_alpha = 1, auto_scale = T, text_size = 6,
       stroke_color = "black",
       show_percentage = F)
ggsave("../plots/Overlap_UpLayer_EQ_Gouhui.venn.pdf", width = 5, height = 4)


Deep_EQ_cor_gene = expr_cor_ct_p[expr_cor_ct_p$cell == "Deep-layer excitatory" & expr_cor_ct_p$pearson_p > 0, "gene"]
Deep_DEG_gene = DEG_CF_p[DEG_CF_p$cell == "Deep-layer excitatory" & DEG_CF_p$logFC > 0, "gene"]

ggvenn(list(A = Deep_EQ_cor_gene, B = Deep_DEG_gene),
       fill_color = c("#66E2D5", "orange"), fill_alpha = 1, auto_scale = T, text_size = 6,
       stroke_color = "black",
       show_percentage = F)
ggsave("../plots/Overlap_DeepLayer_EQ_Gouhui.venn.pdf", width = 5, height = 4)


up_DEG_hp = avg_expr_list[["Upper-layer excitatory"]][, DEG_CF_p[ DEG_CF_p$logFC > 0 & DEG_CF_p$cell == "Upper-layer excitatory", "gene"]]
up_DEG_hp = t(up_DEG_hp)
deep_DEG_hp = avg_expr_list[["Deep-layer excitatory"]][, DEG_CF_p[ DEG_CF_p$logFC > 0 & DEG_CF_p$cell == "Deep-layer excitatory", "gene"]]
deep_DEG_hp = t(deep_DEG_hp)

cortex_DEG_hp = rbind(up_DEG_hp, deep_DEG_hp)

my_col <- colorRampPalette(pals::parula(20))(100)
breaks <- seq(-1, 1, length.out = 100)

pdf("../plots/CF_DEG_Cortex.heatmap.pdf", width = 5, height = 6)
pheatmap(cortex_DEG_hp,
         scale = "row",
         cluster_rows = F,
         clustering_method = "ward.D2",
         show_rownames = F,
         color = my_col,
         gaps_row = nrow(up_DEG_hp),
         breaks = breaks)
dev.off()


# overlap with GWAS genes
phyper_test = function(m, n){
  N = ncol(avg_expr_list[[1]])
  k <- length(intersect(m, n))
  m <- length(m)
  n <- length(n)
  
  OR = (k * (N - m - n + k)) / ((m - k) * (n - k))
  p = phyper(k - 1, m, N - m, n, lower.tail = FALSE)
  c(OR, p)
}

df = read.table("gwas-intelligence.tsv", sep = "\t", header = T)
gene_column <- na.omit(df$MAPPED_GENE)
gene_list <- strsplit(gene_column, "\\s+-\\s+")
all_genes <- unlist(gene_list)
intelligence_gene_set <- unique(all_genes)
intelligence_gene_set <- intelligence_gene_set[intelligence_gene_set != ""]
intelligence_gene_set = unique(intelligence_gene_set)
intelligence_gene_set = intelligence_gene_set[intelligence_gene_set %in% colnames(avg_expr_list[[1]])] # 1051

df = read.table("gwas-cortical_thickness.tsv", sep = "\t", header = T, quote = "\\|")
gene_column <- na.omit(df$MAPPED_GENE)
gene_list <- strsplit(gene_column, "\\s+-\\s+")
all_genes <- unlist(gene_list)
cortical_gene_set <- unique(all_genes)
cortical_gene_set <- cortical_gene_set[cortical_gene_set != ""]
cortical_gene_set = unique(cortical_gene_set)
cortical_gene_set = cortical_gene_set[cortical_gene_set %in% colnames(avg_expr_list[[1]])] # 696

readGWAS = function(x){
  df = read.table(x, sep = "\t", header = T, quote = "\\|")
  gene_column <- na.omit(df[,15])
  gene_list <- strsplit(gene_column, "\\s+-\\s+|,\\s+")
  all_genes <- unlist(gene_list)
  Dementia_gene_set <- unique(all_genes)
  Dementia_gene_set <- Dementia_gene_set[Dementia_gene_set != ""]
  Dementia_gene_set = unique(Dementia_gene_set)
  Dementia_gene_set = Dementia_gene_set[Dementia_gene_set %in% colnames(avg_expr_list[[1]])]
}
AD_gene_set = readGWAS("gwas-Alzheimer.tsv")  # 4585
Schizophrenia_gene_set = readGWAS("gwas-Schizophrenia.tsv")  # 2172
Microcephaly_gene_set = unique(read.table("Microcephaly_gene.txt")$V1) 
Microcephaly_gene_set = Microcephaly_gene_set[Microcephaly_gene_set %in% colnames(avg_expr_list[[1]])] # 72

INTELL_enrich = lapply(unique(DEG_CF_p$cell), function(x){
  m = DEG_CF_p[DEG_CF_p$cell == x, "gene"]
  res = phyper_test(m, intelligence_gene_set)
  res1 = phyper_test(m, cortical_gene_set)
  res2 = phyper_test(m, AD_gene_set)
  c(x, res, res1, res2)
})
INTELL_enrich = do.call(rbind, INTELL_enrich)
INTELL_enrich = as.data.frame(INTELL_enrich)
rownames(INTELL_enrich) = INTELL_enrich$V1
INTELL_enrich = INTELL_enrich[,-1]
INTELL_enrich = apply(INTELL_enrich, 2, as.numeric)
rownames(INTELL_enrich) = unique(DEG_CF_p$cell)
INTELL_enrich = INTELL_enrich[,c(2,4,6)]
colnames(INTELL_enrich) = c("Intelligence", "Cortical thickness", "Alzheimer disease")

p.mat = INTELL_enrich[levels(expr_cor_ct_p_plot$cell),]
star.mat <- matrix("", nrow = nrow(p.mat), ncol = ncol(p.mat))
star.mat[p.mat < 0.01 & p.mat >= 0.001] <- "*"
star.mat[p.mat < 0.001 & p.mat >= 0.0001] <- "**"
star.mat[p.mat < 0.0001] <- "***"

pdf("../plots/Gouhui_GWAS_enrich.heatmap.pdf", 4, 5)
pheatmap(-log10(INTELL_enrich[levels(expr_cor_ct_p_plot$cell),]), cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("gray", "white", "pink","red2"))(100), border_color = NA,
         display_numbers = star.mat,  
         number_color = "black",          
         fontsize_number = 14)
dev.off()
