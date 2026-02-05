# some plot for brain evolution

library(Seurat)
library(reticulate)
use_python("/home/zhangc/miniconda3/bin/python")


NN_cell_percent = read.csv("../Nonneuron_cell_type_proportions_by_species.csv", row.names = 1)
pheatmap::pheatmap(NN_cell_percent)

Neuron_cell_percent = read.csv("../Neuron_cell_type_proportions_by_species.csv", row.names = 1)
pheatmap::pheatmap(Neuron_cell_percent, scale = "row")


####### all cell
ALL = readRDS("../AllCells-del-1.9.RDS")
all_metadata = ALL@meta.data
write.csv(all_metadata, "../AllCells.metadata.csv")

all_metadata1 = all_metadata[sample(rownames(all_metadata), 1000000), ]

cols = as.character(c(pals::alphabet2(26), pals::polychrome(20)))


# qc
ggplot(all_metadata1, aes(batch, nFeature_RNA, fill = batch)) +
  geom_boxplot(outliers = F) + 
  cowplot::theme_cowplot() + xlab("species") + ylab("# gene") +
  geom_hline(yintercept=2000, linetype="dashed", color = "red")+
  scale_fill_manual(values = cols)+
  facet_wrap(organ~orig.ident, ncol = 18)

ggplot(all_metadata1, aes(batch, nFeature_RNA, fill = batch)) +
  geom_violin() + 
  geom_boxplot(fill = "white", width = 0.2, outliers = F) +
  cowplot::theme_cowplot() + xlab("species") + ylab("# gene") +
  geom_hline(yintercept=2000, linetype="dashed", color = "red")+
  scale_fill_manual(values = cols)+
  facet_wrap(organ~orig.ident, ncol = 18)

ggplot(all_metadata1, aes(batch, log10(nCount_RNA+1), fill = batch)) +
  geom_violin() + 
  geom_boxplot(fill = "white", width = 0.2, outliers = F) +
  cowplot::theme_cowplot() + xlab("species") + ylab("# log10(counts)") +
  geom_hline(yintercept=log10(5000), linetype="dashed", color = "red")+
  scale_fill_manual(values = cols)+
  facet_wrap(organ~orig.ident, ncol = 18)


ggplot(all_metadata1, aes(batch, nCount_RNA, fill = batch)) +
  geom_violin() + 
  geom_boxplot(fill = "white", width = 0.2, outliers = F) +
  cowplot::theme_cowplot() + xlab("species") + ylab("# log10(counts)") + ylim(0, 30000)+
  #geom_hline(yintercept=log10(5000), linetype="dashed", color = "red")+
  scale_fill_manual(values = cols)+
  facet_wrap(organ~orig.ident, ncol = 18)


# plot cell number

ggplot(all_metadata, aes(orig.ident, fill = organ)) +
  geom_bar() + 
  cowplot::theme_cowplot() + xlab("Species") + ylab("# cells") +
  scale_fill_manual(values = cols) +
  coord_flip()

ggplot(all_metadata, aes(factor(cell_type, levels = names(sort(table(cell_type)))), fill = orig.ident)) +
  geom_bar() + 
  cowplot::theme_cowplot() + xlab("cell type") + ylab("# cells") +
  scale_fill_manual(values = cols) +
  coord_flip()

ggplot(all_metadata, aes(factor(predicted.id, levels = names(sort(table(predicted.id)))), fill = organ)) +
  geom_bar() + 
  cowplot::theme_cowplot() + xlab("cell type") + ylab("# cells") +
  scale_fill_manual(values = cols) +
  coord_flip()

ggplot(all_metadata, aes(batch,  fill = batch)) +
  geom_bar() +
  cowplot::theme_cowplot() + xlab("") + ylab("# cell (after QC)") +
  geom_hline(yintercept=10000, linetype="dashed", color = "red")+
  scale_fill_manual(values = cols)+
  facet_wrap(organ~orig.ident, ncol = 18)


ggplot(all_metadata, aes(factor(Organ, levels = names(sort(table(Organ)))), fill = predicted.id)) +
  geom_bar() + 
  cowplot::theme_cowplot() + xlab("cell type") + ylab("# cells") +
  scale_fill_manual(values = cols) +
  coord_flip()



all_nongene_meta = all_nongene@meta.data
ggplot(all_nongene_meta, aes(factor(orig.ident, levels = names(sort(table(orig.ident)))), fill = predicted.id)) +
  geom_bar() + 
  cowplot::theme_cowplot() + xlab("cell type") + ylab("# cells") +
  scale_fill_manual(values = ct_col_new) +
  coord_flip()




all_metadata1_counts = table(all_metadata1$orig.ident, all_metadata1$predicted.id, all_metadata1$organ)
all_metadata1_counts = reshape2::melt(all_metadata1_counts)
ggplot(all_metadata1_counts, aes(Var1, Var2, fill = Var3)) +
  geom_point(aes(size = log10(value+1), color = Var3)) +
  cowplot::theme_cowplot() + xlab("") + ylab("# cell") +
  scale_fill_manual(values = cols)

all_metadata1_counts$g = paste(all_metadata1_counts$Var1, all_metadata1_counts$Var2)
ggplot(all_metadata1_counts, aes(x=0, y=value, group = Var3,  color = Var3, fill = Var3)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar("y", start = 0) +
  facet_grid(Var1~., scales = "free_y")# + theme_void()



# plot cell percentage
all_metadata1_counts = table(all_metadata1$orig.ident, all_metadata1$organ)
all_metadata1_counts = reshape2::melt(all_metadata1_counts)
ggplot(all_metadata1_counts, aes(x=0, y=value, group = Var2, fill = Var2)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar("y", start = 0) +
  facet_wrap(Var1~., scales = "free_y", ncol=6) +
  scale_fill_manual(values = cols) + theme_void()




p2 = ggplot(all_metadata1, aes(batch, y = 1, fill = super_cell_type)) +
  geom_bar(position = "fill", stat = "identity") + 
  cowplot::theme_cowplot() + xlab("Batch") + ylab("Percentage") +
  scale_fill_manual(values = cols) &
  facet_wrap(organ~orig.ident, ncol = 18)

p3 = ggplot(all_metadata1, aes(batch, y = 1, fill = cell_type)) +
  geom_bar(position = "fill", stat = "identity") + 
  cowplot::theme_cowplot() + xlab("Batch") + ylab("Percentage") +
  scale_fill_manual(values = cols) &
  facet_wrap(organ~orig.ident, ncol = 18)
p2|p3

ggplot(all_metadata1, aes(batch, y = 1, fill = predicted.id)) +
  geom_bar(position = "fill", stat = "identity") + 
  cowplot::theme_cowplot() + xlab("Batch") + ylab("Percentage") +
  scale_fill_manual(values = cols) +
  coord_polar("y", start=0) +
  facet_wrap(organ~orig.ident, ncol = 18)




# plot1cell

ALL1 = ALL[, sample(colnames(ALL), 300000)]
DimPlot(ALL1, reduction = "ref.umap")

Sys.setenv(RETICULATE_PYTHON = "/home/zhangc/miniconda3/bin/python")
reticulate::use_python("/home/zhangc/miniconda3/bin/python", required = TRUE)

ALL1_rmHuman = subset(ALL1, orig.ident != "Human")
ALL1_rmHuman = RunUMAP(ALL1_rmHuman, reduction = "ref.pca",
                       umap.method = "umap-learn",
                       dims = 1:20, min.dist = .3, n.neighbors = 20L)

DimPlot(ALL1_rmHuman, reduction = "umap") & 
  scale_color_manual(values = as.character(pals::alphabet2(22)))

DimPlot(ALL1_rmHuman, reduction = "ref.umap")& 
  scale_color_manual(values = ct_col_new)

ct_col_new = as.character(pals::alphabet())[-c(5,9,16,24)][c(11,12, 4,14:18, 5:8, 21,22,20, 3,9, 13,2, 10, 1, 19)]
names(ct_col_new) = levels((ALL1_rmHuman$predicted.id))
saveRDS(ct_col_new, "ct_col_new.RDS")


library(plot1cell)

prepare_circlize_data1 <- function(
    seu_obj, 
    scale =0.8
){
  celltypes<-levels(seu_obj)
  cell_colors <- ct_col_new[celltypes]
  data_plot <- get_metadata(seu_obj, color = cell_colors, reductions = "ref.umap", coord_scale = scale)
  data_plot <- cell_order(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}

plot_circlize1 = function (data_plot, do.label = T, contour.levels = c(0.2, 0.3), 
          pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
          col.use = NULL, label.cex = 0.5, repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- as.character(ct_col_new[celltypes])

  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 
                                                            0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) - 
                                                                                                              1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = 0.25, col = "black", 
                             facing = "bending.inside", niceFacing = T)
                 circos.axis(labels.cex = 0.15, col = "black", labels.col = "black")
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  
  sample_cell = sample(600000,1000)
  points(data_plot$x[sample_cell], data_plot$y[sample_cell], pch = 19, col = alpha(data_plot$Colors[sample_cell], 
                                                         0.2), cex = pt.size)
  #points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
   #                                                                      0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

ALL1_rmHuman = subset(ALL1, orig.ident != "Human")
circ_data <- prepare_circlize_data1(ALL1_rmHuman, scale = 0.8)

set.seed(1234)
cluster_colors<-cols[1:length(levels(ALL1_rmHuman))]
batch_colors<-rand_color(length(names(table(ALL1_rmHuman$batch))))
rep_colors<- spe_col[names(table(ALL1_rmHuman$Species))]
organ_colors<-rand_color(length(names(table(ALL1_rmHuman$Organ))))
celltype_colors<-ct_col_new[names(table(ALL1_rmHuman$predicted.id))]

plot_circlize1(circ_data, do.label = F, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', contour.levels = c(0.2, 0.3),
               kde2d.n = 1000, repel = F, label.cex = 0.6)


plot_circlize1(circ_data, do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', contour.levels = c(0.2, 0.3),
              kde2d.n = 1000, repel = T, label.cex = 0.6)

add_track(circ_data, group = "batch", colors = batch_colors, track_num = 2) 
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) 
add_track(circ_data, group = "Organ",colors = organ_colors, track_num = 4) 


all_nongene = readRDS("AllCells-NoGenedata.RDS")
DimPlot(all_nongene, reduction = "ref.umap")


plot(circ_data$x, circ_data$y, 
     axes = FALSE,  bty  = "n",  xlab = "",    ylab = "",
     pch = 19, col = alpha(circ_data$Colors,  0.3), cex = 0.1)


names(organ_colors) = names(table(ALL1_rmHuman$Organ))
circ_data$organ_color = organ_colors[circ_data$Organ]
png("../plots/plot1cell_organ.png", width = 5000, height = 5000, res = 600)
plot(circ_data$x, circ_data$y, 
     axes = FALSE,  bty  = "n",  xlab = "",    ylab = "",
     pch = 19, col = alpha(circ_data$organ_color,  0.3), cex = 0.1)
dev.off()

names(rep_colors) = names(table(ALL1_rmHuman$Species))
circ_data$rep_colors = rep_colors[circ_data$Species]
png("../plots/plot1cell_species.png", width = 5000, height = 5000, res = 600)
plot(circ_data$x, circ_data$y, 
     axes = FALSE,  bty  = "n",  xlab = "",    ylab = "",
     pch = 19, col = alpha(circ_data$rep_colors,  0.3), cex = 0.1)
dev.off()



names(batch_colors) = names(table(ALL1_rmHuman$batch))
circ_data$batch_colors = batch_colors[circ_data$batch]
png("../plots/plot1cell_batch.png", width = 5000, height = 5000, res = 600)
plot(circ_data$x, circ_data$y, 
     axes = FALSE,  bty  = "n",  xlab = "",    ylab = "",
     pch = 19, col = alpha(circ_data$batch_colors,  0.3), cex = 0.1)
dev.off()

pdf("../plots/plot1cell_colorLegend.pdf", width = 7, height = 27)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "", ylab = "", axes = FALSE)
legend("center", 
       legend = names(organ_colors),  
       col = organ_colors,            
       pch = 19,                      
       title = "Organs")
legend("left", 
       legend = names(rep_colors),  
       col = rep_colors,            
       pch = 19,  ncol = 1,              
       title = "Species")
legend("right", 
       legend = names(batch_colors),  
       col = batch_colors,            
       pch = 19,                      
       title = "batch")
dev.off()


DimPlot(ALL1, reduction = "ref.umap", pt.size = 0.01, group.by = "orig.ident",  raster=FALSE, shuffle = T) +
  scale_color_manual(values = rep_colors) +
  NoAxes()

DimPlot(ALL1, reduction = "ref.umap", group.by = "cell_type", pt.size = .01, shuffle = T, raster=FALSE) +
  scale_color_manual(values = as.character(c(pals::cols25(10)))) +
  NoAxes()

DimPlot(ALL1, reduction = "ref.umap", group.by = "Organ", pt.size = .01, shuffle = T, raster=FALSE) +
  scale_color_manual(values = organ_colors) +
  NoAxes()

DimPlot(ALL1, reduction = "ref.umap", group.by = "batch", pt.size = .01, shuffle = T, raster=FALSE) +
  scale_color_manual(values = as.character(c(pals::brewer.dark2(5)))) +
  NoAxes()


FeaturePlot(ALL1, reduction = "umap", features = c("SLC17A6", "INA"), blend = T,pt.size = .1, raster=FALSE) +
  #scale_color_continuous(type = "viridis") +
  NoAxes()


#

all_metadata1_table = table(all_metadata1$super_cell_type, all_metadata1$orig.ident)
all_metadata1_table = t(all_metadata1_table) / colSums(all_metadata1_table)
pheatmap::pheatmap(all_metadata1_table, color = pals::parula(25), border_color = NA, cluster_cols = F)

#

all_metadata1_table = table(all_metadata1$predicted.id, all_metadata1$orig.ident)
all_metadata1_table = t(all_metadata1_table) / colSums(all_metadata1_table)
pheatmap::pheatmap(all_metadata1_table, color = pals::parula(25), border_color = NA, cluster_cols = F)


Clust_brain_region = function(x){
  all_metadata1_CC = all_metadata1[all_metadata1$organ == x, ]
  all_metadata1_table = table(all_metadata1_CC$cell_type, all_metadata1_CC$orig.ident)
  all_metadata1_table = t(all_metadata1_table) / colSums(all_metadata1_table)
  colnames(all_metadata1_table) = paste(x, colnames(all_metadata1_table) )
  all_metadata1_table
}

Brain_cell_pct = lapply(unique(all_metadata1$organ), Clust_brain_region)
Brain_cell_pct = do.call(cbind, Brain_cell_pct)
pheatmap::pheatmap(Brain_cell_pct, color = pals::parula(25), border_color = NA, cluster_cols = F)

pheatmap::pheatmap(Brain_cell_pct[,seq(1, 60, 10)], color = pals::parula(25), border_color = NA,
                   #clustering_method = "ward.D2",
                   cluster_cols = F)



### Clustering based on PC
ALL = readRDS("../AllCells-del-1.9.RDS")
ALL1 = ALL[, sample(colnames(ALL), 300000)]
ref_pca_data <- ALL1@reductions$ref.pca@cell.embeddings

Macaque = subset(ALL1, orig.ident == "Macaque")
Macaque[["ref.pca"]] = Macaque[["pca"]]
Macaque_ref_pca = Macaque@reductions$ref.pca@cell.embeddings

ref_pca_data = rbind(ref_pca_data, Macaque_ref_pca[,1:30])

metadata <- ALL1@meta.data
pca_metadata <- cbind(metadata, ref_pca_data[rownames(metadata), ])
pca_median <- pca_metadata %>%
  group_by(orig.ident, predicted.id) %>%
  dplyr::summarise(across(starts_with("refpca"), median))

pca_median <- pca_metadata[,c("orig.ident", "predicted.id", paste0("refpca_", 1:30))] %>%
  group_by(orig.ident, predicted.id) %>%
  dplyr::summarise_all(median)


ggplot(pca_median, aes(refpca_1, refpca_2, color = predicted.id, fill = predicted.id)) +
  geom_point(stroke =1 ) +
  facet_wrap(.~predicted.id)

library(dendextend)
library(circlize)
pca_values <- pca_median[, grep("^refpca_", colnames(pca_median))]
dist_matrix <- dist(pca_values, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "complete")
order_of_leaves <- hclust_result$order
labels <- paste(pca_median$orig.ident, pca_median$predicted.id, sep = "_")
ordered_labels <- labels[order_of_leaves]
dendrogram <- as.dendrogram(hclust_result)
dendrogram <- dendrogram %>% set("labels", ordered_labels)

par(mar = c(10, 4, 4, 2))  
plot(dendrogram, main = "Hierarchical Clustering of PCA Median Values", yaxt = "n")

pca_values = as.data.frame(pca_median)
rownames(pca_values) = paste(pca_values$predicted.id, pca_values$orig.ident)

cell_ann = pca_values[,1:2]
pca_values = pca_values[, -c(1:2)]

pdf("../plots/cell_type_speice_hclust.pdf", 30, 10)
pheatmap::pheatmap(t(pca_values), color = pals::parula(25), border_color = NA, cluster_rows = F,
                   annotation_col = cell_ann, clustering_method = "ward.D", fontsize_col = 6,
                   annotation_colors = ann_colors)
dev.off()



hc = as.dendrogram(hclust(dist(pca_values), method = "ward.D"))
hc1 = as.dendrogram(hclust(dist(pca_values), method = "ward.D2"))
pdf("../plots/cell_type_speice_hclust_circlize.pdf", 20, 20)
circlize_dendrogram(hc, labels_track_height = .3)
circlize_dendrogram(hc1, labels_track_height = .3)
dev.off()



#
pca_cell_cluster = reshape2::melt(as.data.frame(pca_median))
pca_cell_cluster$y = paste0(pca_cell_cluster$predicted.id, pca_cell_cluster$variable)

pca_cell_cluster = dcast(pca_cell_cluster, orig.ident ~ y, value.var = "value" )

rownames(pca_cell_cluster) = pca_cell_cluster[, 1]
pca_cell_cluster = pca_cell_cluster[, -1]

pdf("../plots/cell_type_speice_expr_hclust.pdf", 30, 10)
pheatmap::pheatmap(pca_cell_cluster, color = pals::parula(25), border_color = NA, show_colnames = F,
                   annotation_row = annotation_row[,1, drop = F],
                   annotation_colors = ann_colors)
dev.off()
plot(hclust(dist(pca_cell_cluster), method = "ward.D2"), hang = -1)



# pan-markers
ALL1 = JoinLayers(ALL1)
pan_marker = FindAllMarkers(ALL1, 
                            logfc.threshold = log2(1.5), 
                            min.diff.pct = 0.5, min.pct = 0.3,
                            only.pos = T)

library(dplyr)
pan_marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> top2

FeaturePlot(ALL1, "P2RY12")

DotPlot(ALL1, unique(top2$gene)) + RotatedAxis() + coord_flip()
DotPlot(ALL1, unique(top2$gene), split.by = "organ", cols = organ_colors) + RotatedAxis()

#
Idents(ALL1) = "cell_type"
pan_marker = FindAllMarkers(ALL1, 
                            logfc.threshold = log2(1.5), 
                            min.diff.pct = 0.5, min.pct = 0.3,
                            only.pos = T)

pan_marker$pct.diff = pan_marker$pct.1 - pan_marker$pct.2
pan_marker %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = pct.diff) -> top2

FeaturePlot(ALL1, "INA")

DotPlot(ALL1, unique(top2$gene)) + RotatedAxis() + coord_flip() +
  scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))

DotPlot(ALL1, unique(top2$gene), split.by = "orig.ident", cols = col1) + RotatedAxis() + coord_flip()


### evolution distance
# PC
spe = unique(ALL1$orig.ident)
ct = unique(ALL1$predicted.id)

spe_dist = lapply(ct, function(x){
  pca_median1 = pca_median[pca_median$predicted.id == x, ]
  as.numeric(cor(t(pca_median1[,-(1:2)])))
  #as.numeric(dist((pca_median1[,-(1:2)])))
})
names(spe_dist) = ct

spe_dist = melt(spe_dist)
ggplot(spe_dist, aes(L1, value))+
  geom_boxplot() + coord_flip() + cowplot::theme_cowplot() + ylim(0.7,1)

# 
pca_median_batch <- pca_metadata %>%
  group_by(orig.ident, batch, predicted.id) %>%
  dplyr::summarise(across(starts_with("refpca"), median))



spe_dist = function(x, y, z){
  tryCatch({
    pca_median1 = pca_median_batch[pca_median_batch$predicted.id == x, ]
    pca_median2 = pca_median1[pca_median1$orig.ident %in% c(y, z), ]
    n_batch1 = nrow(pca_median1[pca_median1$orig.ident == y, ])
    n_batch2 = nrow(pca_median1[pca_median1$orig.ident == z, ])
    cor_dist = 1 - cor(t(pca_median2[,-(1:3)]))
    batch_dist1 = (cor_dist[2:n_batch1, 1:(n_batch1-1)])
    batch_dist2 = (cor_dist[(n_batch1+2):(n_batch1+n_batch2), (n_batch1+1):(n_batch1+n_batch2-1)])
    spe_dist = as.numeric(cor_dist[(n_batch1+1):(n_batch1+n_batch2), 1:n_batch1])
    result = mean(spe_dist) / mean(c(batch_dist1, batch_dist2))
    return(result)
  },
  error = function(e) {return(NA)})
}
spe_dist("Eccentric medium spiny neuron", "Mouse", "Rat")

spe = spe[spe != "Lemur"] # only 1 individual
ct = ct[ct != "Miscellaneous"]
spe_dist_res = lapply(ct, function(x){
  res = c()
   for (i in spe){
     for (j in spe){
       if(i != j){
         result = spe_dist(x, i, j)
         res = c(res, result)
       }
     }
   }
  res
})
names(spe_dist_res) = ct

spe_dist_res = melt(spe_dist_res)
ggplot(spe_dist_res, aes(L1, value, fill = L1))+
  geom_boxplot(outliers=F) + coord_flip() + cowplot::theme_cowplot() +
  xlab("") + ylab("Distance")


