# cell percentage
library(ape)
library(ggtree)
library(dplyr)
library(tidyr)

cpt = read.csv("celltyperatio_matrix.csv", row.names = 1)
phe = read.csv("phenodata.csv", row.names = 1)
pheno <- read.csv("EQ_adj.csv", row.names = 1)
pheno = pheno[rownames(phe)[-1], ]
pheno = cbind(pheno, phe[-1, 5:7])

phe = pheno
phe$log_Weight = log10(phe$Weight)
phe$log_Brain.Weight = log10(phe$BW)

ph = pheatmap(cpt, scale = "none", clustering_method = "ward.D2")

row_tree <- ph$tree_row        # 这是 hclust 对象
plot(row_tree) 

ph$tree_row$labels[ph$tree_row$order]

pheatmap(phe[ph$tree_row$labels[ph$tree_row$order], ], scale = "column", cluster_rows = F)

cpt = t(t(cpt)/as.numeric(cpt[1,]))
hc <- hclust(dist(cpt), method = "ward.D2")
sorted_labels <- hc$labels[hc$order]
plot(hc, hang = -1, horiz = TRUE)

par(mar = c(5, 4, 4, 12) + 0.1) 
pdf("cell_percentage_hclust.pdf",6,8)
plot(rev(as.dendrogram(hc)) ,         
     horiz = TRUE)
dev.off()

pheatmap(phe[sorted_labels, ], cluster_rows = F,scale = "column")


cpt = read.csv("celltyperatio_matrix.csv", row.names = 1)
cell_order = names(sort(colMeans(cpt)))
cpt.m = melt(as.matrix(cpt))

cell_color = readRDS("~/proj/Brain_scRNA/Scripts/cell_color.rds")
names(cell_color) = str_replace_all(names(cell_color), " ", ".")
names(cell_color) = str_replace_all(names(cell_color), "-", ".")
cpt.m$Var1 = factor(cpt.m$Var1, levels = rev(sorted_labels))
cpt.m$Var2 = factor(cpt.m$Var2, levels = cell_order)
ggplot(cpt.m, aes(Var1, value, fill = Var2)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = as.character(c(pals::cols25(25)))) 

ct_col_new1 = ct_col_new
names(ct_col_new1) = str_replace_all(names(ct_col_new1), " ", ".")
names(ct_col_new1) = str_replace_all(names(ct_col_new1), "-", ".")
ggplot(cpt.m, aes(Var1, value, fill = Var2)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  #scale_fill_manual(values = as.character((pals::alphabet(26)[1:22]))) +
  scale_fill_manual(values = ct_col_new1) +
  scale_x_discrete(position = "top") +  
  theme(axis.text.y = element_text(margin = margin(r = 10)))  +
  cowplot::theme_cowplot()
ggsave("cell_percentage_species.heatmap.pdf", width = 10, height = 7)


phe$species = rownames(phe)
phe.m = melt(phe[,c(10,6,8, 11)])
phe.m$species = factor(phe.m$species, levels = rev(sorted_labels))

cpt$Neuron = rowSums(cpt[, c(3,4,5,6,9,10:14,16,19:21)])
cpt$Glia = rowSums(cpt[, c(1,2,15,17,18)])
cpt$Non_Neruonal = rowSums(cpt[, c(7,8,22)])

cpt.m = melt(as.matrix(cpt[,c("Neuron", "Glia", "Non_Neruonal")]))
cpt.m$Var1 = factor(cpt.m$Var1, levels = rev(sorted_labels))
ggplot(cpt.m, aes(Var1, value, fill = Var2)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = as.character(c(pals::cols25(25)))) 

ggplot(cpt.m, aes(Var1, value, fill = Var2)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  #scale_fill_manual(values = as.character((pals::alphabet(26)[1:22]))) +
  #scale_fill_manual(values = cell_color) +
  scale_x_discrete(position = "top") +  
  theme(axis.text.y = element_text(margin = margin(r = 10)))  +
  cowplot::theme_cowplot()


ggplot()+
  geom_point(data = phe.m[phe.m$variable == "log_Brain.Weight", ], aes(species, y = 1, size = value, color = value)) +
  scale_color_viridis()+
  geom_text(data = phe.m[phe.m$variable == "Cortical.convolutions", ], aes(species, y = 2, label = ifelse(value == 1, "✓", "✗")), 
            size = 4, fontface = "bold", color = ifelse(phe.m[phe.m$variable == "Cortical.convolutions", ]$value == 1, "red3", "black"),
            show.legend = FALSE) +
  geom_text(data = phe.m[phe.m$variable == "Hibernation", ], aes(species, y = 3, label = ifelse(value == 1, "✓", "✗")), 
            size = 4, fontface = "bold", color = ifelse(phe.m[phe.m$variable == "Hibernation", ]$value == 1, "red3", "black"),
            show.legend = FALSE) +
  coord_flip(ylim = c(0,4)) +
  theme_void()
ggsave("phentoy_spe.clust_percent.dotplot.pdf", device = cairo_pdf, width = 1.5, height = 7)




cpt = read.csv("celltyperatio_matrix.csv", row.names = 1)
phe_cor = lapply(colnames(phe), function(x){
  phe1 = phe[,x]
  res = apply(cpt, 2, function(y){
    pearson_test <- cor.test(as.numeric(y), phe1, method = "pearson")
    pearson_cor = as.numeric(pearson_test$estimate)
    pearson_p = pearson_test$p.value
    c(pearson_cor, pearson_p)
  })
  res = as.data.frame(t(res))
  colnames(res) = c( "pearson_r", "pearson_p")
  res$pearson_padj = p.adjust(res$pearson_p, method = "fdr")
  res$phenotype = x
  res$cell = rownames(res)
  res
})

phe_cor = do.call(rbind, phe_cor)

phe_cor$stars <- cut(phe_cor$pearson_padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

ggplot(phe_cor, aes(x=cell, y=phenotype, fill= pearson_r)) +
  geom_tile() + scale_fill_gradientn(colours = c("lightblue", "white", "red")) +
  geom_text(aes(label=stars), vjust = 1, color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="pearson_r") + #geom_vline(xintercept=1.5, size=1.5, color="grey50") + 
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))






