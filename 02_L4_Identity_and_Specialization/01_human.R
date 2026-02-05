# ==============================================================================
# Step 01: Human Cortex L4 Analysis
# ==============================================================================

Human_sobj <- readRDS("/home/zhangc/proj/Brain_scRNA/test/Human_sobj_downsampled.RDS")
Human_sobj_cortex <- subset(Human_sobj, supercluster_term %in% c("Deep-layer intratelencephalic", "Deep-layer corticothalamic and 6b", "Deep-layer near-projecting", "Upper-layer intratelencephalic") | cluster_id %in% c("113", "114", "115", "117", "118"))

Human_sobj_cortex <- subset(Human_sobj_cortex, tissue == "cerebral cortex")
set.seed(1)
Human_sobj_cortex <- Human_sobj_cortex[, sample(ncol(Human_sobj_cortex), 10000)]
Human_sobj_cortex <- Human_sobj_cortex[rownames(Human_sobj_cortex) %in% rownames(ALL1), ]


Human_sobj_cortex <- NormalizeData(Human_sobj_cortex)
Human_sobj_cortex <- FindVariableFeatures(Human_sobj_cortex, nfeatures = 3000)

Human_sobj_cortex <- ScaleData(Human_sobj_cortex)
Human_sobj_cortex <- RunPCA(Human_sobj_cortex)

Human_sobj_cortex <- RunUMAP(Human_sobj_cortex, dims = 1:10)


CellAnn <- read.table("human_cortex_cellAnno.txt", sep = "\t", header = T, row.names = 1)
Human_sobj_cortex$MTG.Label <- CellAnn[Human_sobj_cortex$cluster_id, "Transferred.MTG.Label"]

DimPlot(Human_sobj_cortex, group.by = c("MTG.Label", "cluster_id"), label = T)
DimPlot(Human_sobj_cortex, group.by = "supercluster_term")


Human_sobj_cortex <- FindNeighbors(Human_sobj_cortex, dims = 1:10)
Human_sobj_cortex <- FindClusters(Human_sobj_cortex, resolution = 0.4)

DimPlot(Human_sobj_cortex, label = T)

DimPlot(Human_sobj_cortex, group.by = c("MTG.Label", "cluster_id"), label = T)

new.cluster.ids <- c(
  "L2/3 IT", "L2/3 IT", "L5 IT", "L2/3 IT", "L6 IT", "L4 IT",
  "L6 CT", "L2/3 IT", "L6b", "L5 IT", "L6 IT Car3", "L2/3 IT", "L5/6 NP", "L5 ET"
)
names(new.cluster.ids) <- levels(Human_sobj_cortex)
Human_sobj_cortex <- RenameIdents(Human_sobj_cortex, new.cluster.ids)

DimPlot(Human_sobj_cortex, label = T)
Human_sobj_cortex$cell <- Idents(Human_sobj_cortex)

FeaturePlot(Human_sobj_cortex, c("FOXP2", "VIPR2", "COBLL1", "RORB"))

DotPlot(Human_sobj_cortex, c("CUX2", "VIPR2", "COBLL1", "RORB", "TRABD2A")) & coord_flip()

saveRDS(Human_sobj_cortex, "Human_sobj_cortex.RDS")

# Science

dir.path <- "../data/Science2023/"
read_Science2023 <- function(s) {
  counts <- readRDS(file.path(dir.path, paste0(s, "_mat.RDS")))

  cell.meta <- readRDS(file.path(dir.path, paste0(s, "_meta.RDS")))
  cell.meta <- as.data.frame(cell.meta)
  rownames(cell.meta) <- cell.meta$sample_id

  counts <- counts[, cell.meta$sample_id]

  scRNA <- CreateSeuratObject(counts = counts, meta.data = cell.meta)
  scRNA <- NormalizeData(scRNA)
  scRNA
}

Science2023_human <- read_Science2023("human")
set.seed(1)
Science2023_human <- Science2023_human[, sample(ncol(Science2023_human), 10000)]
Science2023_human_cortex <- subset(Science2023_human, subclass %in% c("L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 IT", "L6 CT", "L6 IT Car3", "L6b"))
Science2023_human.expr <- AverageExpression(Science2023_human_cortex, group.by = "subclass")$RNA


Human_sobj_cortex.expr <- AverageExpression(Human_sobj_cortex)$RNA

common_gene <- intersect(rownames(Science2023_human.expr), rownames(Human_sobj_cortex.expr))
common_gene <- common_gene[common_gene %in% VariableFeatures(Human_sobj_cortex)]

Human_avg_expr.scale <- t(scale(t(Human_sobj_cortex.expr[common_gene, ])))
Science2023_avg_expr.scale <- t(scale(t(Science2023_human.expr[common_gene, ])))

Human_avg_expr.scale <- na.omit(Human_avg_expr.scale)
Science2023_avg_expr.scale <- na.omit(Science2023_avg_expr.scale)

common_gene <- intersect(rownames(Human_avg_expr.scale), rownames(Science2023_avg_expr.scale))
Expr_cor <- cor(Human_avg_expr.scale[common_gene, ], Science2023_avg_expr.scale[common_gene, ])

pdf("../plots/compare_Science2023_L4.cluster.pdf", 5, 4)
pheatmap::pheatmap(Expr_cor[c(1, 4, 2, 3, 9, 8, 5, 7, 6), c(1, 2, 4, 7, 3, 5, 6, 8, 9)],
  cluster_rows = F, cluster_cols = F,
  color = colorRampPalette(c("white", "white", "gray95", "lightpink", "red3"))(30),
  # color =  colorRampPalette(pals::parula(20))(30),
  border_color = NA
)
dev.off()

#
merged <- readRDS("../data/Science2023/merged.RDS")
