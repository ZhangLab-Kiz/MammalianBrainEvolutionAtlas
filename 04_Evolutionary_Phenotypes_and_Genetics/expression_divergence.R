# molecular distance
library(dplyr)


ALL = readRDS("AllCells_include_Human_7.22.RDS")

ALL$Species[ALL$orig.ident == "Human"] <- "Homo_sapiens"
Avg_Expr_Species = AverageExpression(ALL, group.by = "Species")
Avg_Expr_Species = Avg_Expr_Species$RNA


ALL$Species[ALL$orig.ident == "Human"] <- "Homo_sapiens"
DimPlot(ALL, reduction = "ref.umap", group.by = "Species")
DimPlot(ALL, reduction = "ref.umap", group.by = "predicted.id", label = T)
DimPlot(ALL, reduction = "ref.pca", group.by = "Species")

ref_pca_data <- ALL@reductions$ref.pca@cell.embeddings

#### Cell Phylogenies
metadata <- ALL@meta.data

pca_metadata <- cbind(metadata, ref_pca_data)
pca_median <- pca_metadata[,c("Species", "predicted.id", paste0("refpca_", 1:30))] %>%
  group_by(Species, predicted.id) %>%
  dplyr::summarise(across(starts_with("refpca"), median))


ggplot(pca_median, aes(refpca_1, refpca_2, color = predicted.id, fill = predicted.id)) +
  geom_point(stroke =1 ) +
  facet_wrap(.~predicted.id)


pca_values = as.data.frame(pca_median)
rownames(pca_values) = paste(pca_values$predicted.id, pca_values$Species)

cell_ann = pca_values[,1:2]
pca_values = pca_values[, -c(1:2)]



spe_col = readRDS("species_color_new.rds")
ct_col = readRDS("ct_col_new.RDS")


ann_colors = list(Species = spe_col,
                  predicted.id = ct_col)


pdf("cell_type_species_hclust.pdf", 20, 13)
pheatmap::pheatmap(t(pca_values[, 1:30]), color = pals::parula(25),border_color = NA, cluster_rows = F, show_colnames = F, 
                   annotation_col = cell_ann, clustering_method = "ward.D", fontsize_col = 6, clustering_distance_cols = "correlation",
                   annotation_colors = ann_colors)
dev.off()


pdf("cell_type_species_correlation.pdf", 18, 15)
pheatmap::pheatmap(cor(t(pca_values[, 1:30])), color = pals::parula(25), border_color = NA, show_rownames = F, show_colnames = F,
                   annotation_col = cell_ann, annotation_row = cell_ann, clustering_method = "ward.D", fontsize_col = 6, 
                   clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", legend = T,
                   cellwidth = .5,  
                   cellheight = .5,
                   treeheight_row = 30,  
                   treeheight_col = 30,
                   annotation_colors = ann_colors)
dev.off()

## transcriptomic divergence, compare with Human

ct = unique(ALL$predicted.id)
spe = unique(ALL$Species)

spe_cor_human = lapply(ct, function(x){
  print(x)
  cor_human = c()
  for(i in spe[-38]){
    res = NA
    tryCatch({
      avg_expr = AverageExpression(subset(ALL, predicted.id == x & Species %in% c(i, "Homo_sapiens")), group.by = c("Species"))$RNA
      res = cor(avg_expr[,1], avg_expr[,2])
    },
    error = function(e) {return(NA)})
    cor_human = c(cor_human, res)
  }
  names(cor_human) = spe[-38]
  cor_human
})

saveRDS(spe_cor_human, "spe_cor_human.RDS")
names(spe_cor_human) = ct

spe_cor_human1 = do.call(cbind, spe_cor_human)
pdf("../plots/expr_cor_with_human_all.pdf", 10, 6)
pheatmap::pheatmap(t(spe_cor_human1), scale = "none", na_col = "gray", cluster_rows = FALSE, cluster_cols = FALSE,
                   border_color = NA,
                   color = colorRampPalette(pals::parula(8))(100))
dev.off()


spe_order = c("Macaca_mulatta", "Callithrix_jacchus", "Colobus_guereza", "Lemur_catta",
              "Mus_musculus", "Rattus_norvegicus", "Meriones_unguiculatus", "Phodopus_sungorus", "Cavia_porcellus", "Chinchilla_lanigera", "Heterocephalus_glaber","Oryctolagus_cuniculus",
              "Canis_lupus_familiaris", "Vulpes_lagopus", "Nyctereutes_procyonoides", "Mustela_putorius_furo", "Meles_meles", "Felis_catus",
              "Bos_taurus", "Capra_hircus", "Sus_scrofa",
              "Hipposideros_armiger",
              "Suncus_murinus")
custom_breaks <- seq(0.1, 0.4, length.out = 100)

pdf("../plots/expr_cor_with_human_selected_col1.pdf", 10, 6)
pheatmap::pheatmap(t(1-spe_cor_human1[spe_order,
                                    c(8, 3, 1,12,21, 2,5,6,9,14,17,18,20,13, 19,4,7,10,11,15,16)]), scale = "none",
                   na_col = "gray", cluster_rows = FALSE, cluster_cols = FALSE,
                   breaks = custom_breaks,
                   border_color = NA,
                   color = colorRampPalette(pals::parula(20))(100))
dev.off()

cell_dis_order = colnames(spe_cor_human1)[c(8, 3, 1,12,21, 2,5,6,9,14,17,18,20,13, 19,4,7,10,11,15,16)]


# species-wise divergence
spe_cor_MEAN = lapply(ct, function(x){
  expr_mat = t(Gene_avg_counts[[x]])
  expr_mat = expr_mat[,-c(3,9:13,15,25,26,30,32,34,35,37:38)]
  expr_dist = 1-cor(expr_mat)
  expr_dist[upper.tri(expr_dist)]
})

spe_cor_MEAN_mat = do.call(cbind, spe_cor_MEAN)
colnames(spe_cor_MEAN_mat) = ct
spe_cor_MEAN_mat = reshape2::melt(spe_cor_MEAN_mat)

ggplot(na.omit(spe_cor_MEAN_mat), aes(reorder(Var2, value, FUN = median, decreasing = T), value,  fill = Var2))+
  geom_boxplot(outliers=F, alpha = 1) + cowplot::theme_cowplot() +
  geom_jitter(size = 0.1, shape=20)+
  xlab("") + ylab("Transcriptomic divergence \n (1- Pearson correlation)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values = ct_col)
ggsave("cell_type_cross_spe_1minusCor_all2all.pdf", width = 10, height = 6)

library(gghalves)

ggplot(na.omit(spe_cor_MEAN_mat), aes(reorder(Var2, value, FUN = mean, decreasing = T), value, fill = Var2, color = Var2))+
  geom_half_violin(side = "l", position = "dodge") + 
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3, 
    shape = 21,
    fill = "red3",
    color = "white",
    position = position_nudge(x = -0.12) 
  ) +
  geom_half_point(side = "r", size = 0.5, shape = 20) + 
  cowplot::theme_cowplot() +
  xlab("") + ylab("Transcriptomic divergence \n (1- Pearson correlation)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values = ct_col_new) +
  scale_color_manual(values = ct_col_new)
ggsave("cell_type_cross_spe_1minusCor_all2all_violin.pdf", width = 10, height = 6)



ggplot(na.omit(spe_cor_MEAN_mat), aes(reorder(cell_class, value, FUN = mean, decreasing = T), value, fill = cell_class))+
  geom_violin() + 
  cowplot::theme_cowplot() +
  stat_summary(fun = median, geom = "point", colour = "red1", size = 2) +
  xlab("") + ylab("Transcriptomic divergence \n (1- Pearson correlation)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none") +
  scale_fill_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62"))
ggsave("../plots/cell_type_cross_spe_1minusCor_all2all_3class.violin.pdf", width = 3, height = 4)

wilcox.test(na.omit(spe_cor_MEAN_mat[spe_cor_MEAN_mat$cell=="Neurons", 3, drop=T]), na.omit(spe_cor_MEAN_mat[spe_cor_MEAN_mat$cell!="Neurons", 3, drop=T])) # p-value < 2.2e-16




## scatterplot
# Non-neuron
avg_expr = AverageExpression(subset(ALL, predicted.id == "Oligodendrocyte" & Species %in% c("Suncus_murinus", "Homo_sapiens")), group.by = c("Species"))$RNA
avg_expr1 = as.data.frame(avg_expr)
p1 = ggplot(log(avg_expr1+1), aes(`Suncus-murinus`, `Homo-sapiens`)) +
  geom_point(size=0.5) +
  theme_classic()+ggtitle(paste0("Pearson's cor: ", cor(avg_expr[,1], avg_expr[,2])))

avg_expr = AverageExpression(subset(ALL, predicted.id == "Microglia" & Species %in% c("Suncus_murinus", "Homo_sapiens")), group.by = c("Species"))$RNA
avg_expr1 = as.data.frame(avg_expr)
p2 = ggplot(log(avg_expr1+1), aes(`Suncus-murinus`, `Homo-sapiens`)) +
  geom_point(size=0.5) +
  theme_classic()+ggtitle(paste0("Pearson's cor: ", cor(avg_expr[,1], avg_expr[,2])))


# Neuron
avg_expr = AverageExpression(subset(ALL, predicted.id == "Upper-layer excitatory" & Species %in% c("Suncus_murinus", "Homo_sapiens")), group.by = c("Species"))$RNA
avg_expr1 = as.data.frame(avg_expr)
p3 = ggplot(log(avg_expr1+1), aes(`Suncus-murinus`, `Homo-sapiens`)) +
  geom_point(size=0.5) +
  theme_classic()+ggtitle(paste0("Pearson's cor: ", cor(avg_expr[,1], avg_expr[,2])))

avg_expr = AverageExpression(subset(ALL, predicted.id == "Medium spiny neuron" & Species %in% c("Suncus_murinus", "Homo_sapiens")), group.by = c("Species"))$RNA
avg_expr1 = as.data.frame(avg_expr)
p4 = ggplot(log(avg_expr1+1), aes(`Suncus-murinus`, `Homo-sapiens`)) +
  geom_point(size=0.5) +
  theme_classic()+ggtitle(paste0("Pearson's cor: ", cor(avg_expr[,1], avg_expr[,2])))


p1|p2|p3|p4
ggsave("../plots/Expr_Divergence_Suncus_murinus_Homo_sapiens_Oligodendrocyte_Microglia_Upperlayerexcitatory_Mediumspinyneuron.pdf", width = 15, height = 4)



