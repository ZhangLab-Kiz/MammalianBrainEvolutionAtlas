# ==============================================================================
# Step 05: L4 Gene Expression Analysis
# ==============================================================================

library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

Idents(obj_filtered) <- "final.celltype"

DefaultAssay(obj_filtered) <- "RNA"

obj_filtered <- JoinLayers(obj_filtered)

# ==============================================================================
# DEG Analysis
# ==============================================================================

obj_pri <- subset(obj_filtered, subset = Order == "Primates")

markers_L4 <- FindMarkers(obj_pri,
  ident.1 = "L4 IT",
  ident.2 = NULL,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

markers_L4$gene <- rownames(markers_L4)

markers_L4$sig <- ifelse(markers_L4$p_val_adj < 0.05 & abs(markers_L4$avg_log2FC) > 0.5, "Significant", "No")
markers_L4$direction <- "No"
markers_L4$direction[markers_L4$sig == "Significant" & markers_L4$avg_log2FC > 0.5] <- "Up"
markers_L4$direction[markers_L4$sig == "Significant" & markers_L4$avg_log2FC < -0.5] <- "Down"

write.csv(markers_L4, file.path(DIR_FIG, "DEG_Primate_L4.csv"))

# ==============================================================================
# GO Enrichment
# ==============================================================================

genes_to_test <- markers_L4 %>%
  filter(direction == "Up") %>%
  pull(gene)

ego <- enrichGO(
  gene = genes_to_test,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.01
)

p_go <- dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment of Primate L4 Specific Genes") +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(DIR_FIG, "GO_Enrichment_L4.pdf"), p_go, width = 7, height = 8)


# ==============================================================================
# Volcano Plot
# ==============================================================================

markers_L4 <- read.csv("Query/Query_Figure/Paper/DEG_Primate_L4.csv")

plot_data <- markers_L4
plot_data$logP <- -log10(plot_data$p_val_adj)

max_p_val <- max(plot_data$logP[is.finite(plot_data$logP)])
plot_data$logP[is.infinite(plot_data$logP)] <- max_p_val + 20

logfc_cut <- 0.25
pval_cut <- 0.05

plot_data$is_sig <- ifelse(plot_data$p_val_adj < pval_cut & plot_data$avg_log2FC > logfc_cut, "Yes", "No")

manual_genes <- c("RORB", "FOXP2", "EPHA3", "RARB", "DCC")
top5_genes <- plot_data %>%
  filter(is_sig == "Yes") %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene)
label_genes <- unique(c(top5_genes, manual_genes))

p_vol <- ggplot() +
  geom_point(
    data = subset(plot_data, is_sig == "No"),
    aes(x = avg_log2FC, y = logP),
    color = "grey85", size = 1, alpha = 0.9
  ) +
  geom_point(
    data = subset(plot_data, is_sig == "Yes"),
    aes(x = avg_log2FC, y = logP),
    color = "#E41A1C", size = 2, alpha = 0.9
  ) +
  scale_color_gradientn(
    name = expression(-log[10](italic(P)[adj])),
    colors = c("#FF9999", "#E41A1C", "#800000")
  ) +
  geom_text_repel(
    data = subset(plot_data, gene %in% label_genes),
    aes(x = avg_log2FC, y = logP, label = gene),
    size = 3.5,
    fontface = "italic",
    color = "black",
    box.padding = 1.2,
    point.padding = 0.6,
    force = 5,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey30",
    segment.size = 0.3
  ) +
  geom_vline(xintercept = c(logfc_cut), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "grey60") +
  scale_x_continuous(name = expression(log[2] ~ Fold ~ Change)) +
  scale_y_continuous(name = expression(-log[10](italic(P)[adj]))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.position = "right"
  )

ggsave(file.path(DIR_FIG, "Volcano_L4.pdf"), p_vol, width = 6, height = 6)
