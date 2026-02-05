# Brain developmemt data
# https://github.com/yimingsun12138/BrainEvoDevo

Brain_dev = readRDS("Brain_integrated_RNA_Seurat.rds")

DimPlot(subset(Brain_dev, species %in% c("macaque", "mouse")) , group.by = "cell_type", label = T, raster=FALSE, pt.size = 0.01) & NoAxes() &
  scale_color_manual(values = pals::cols25(20)) 
ggsave("../plots/Deve_celltype.dimplot.pdf", width = 10, height = 8)
ggsave("../plots/Deve_celltype.dimplot.png", width = 10, height = 8)

DimPlot(subset(Brain_dev, species %in% c("macaque", "mouse")) , group.by = "species", raster=FALSE) & NoAxes() &
  scale_color_manual(values = as.character(pals::brewer.dark2(3))) 
ggsave("../plots/Deve_species.dimplot.png", width = 10, height = 8)

Brain_dev1 = subset(Brain_dev, cell_type %in% c("RG", "nIP", "ExN", "ExM", "ExUp", "ExDp") & species %in% c("macaque", "mouse")) 
Brain_dev1$cell_type = factor(Brain_dev1$cell_type, levels = c("RG", "nIP", "ExN", "ExM", "ExUp", "ExDp"))
Idents(Brain_dev1) = "cell_type"
DefaultAssay(Brain_dev1) = "RNA"

DimPlot(Brain_dev1, group.by = "cell_type", split.by = "species")

genes = c("MCPH1", "ALG13", "ZRANB3", "ANO5", "WDR11", "COG6", "MARK1", "PIK3CB", "MSH3", "DNAJB14", "GNPDA2", "DYNC1LI1", "DYNC1H1", "HOPX", "PAK1")

DoHeatmap(Brain_dev1, c("HOPX"), split.by = "species")

dev_DEG = FindMarkers(Brain_dev1, group.by = "species",
                      ident.1 = "macaque", ident.2 = "mouse")
dev_DEG$gouhui_DEG = "No"
dev_DEG$gouhui_DEG[rownames(dev_DEG) %in% c(Up_DEG_gene, Deep_DEG_gene)] = "Gouhui"
dev_DEG$gouhui_DEG[rownames(dev_DEG) %in% expr_cor_ct_p[expr_cor_ct_p$cell %in% c("Upper-layer excitatory", "Deep-layer excitatory"), "gene"]] = "EQ"

dev_DEG_p = dev_DEG[dev_DEG$avg_log2FC > 0 & dev_DEG$p_val_adj < 0.05, ] # 261/348 ~ 72.9%

genes = c("MCPH1", "ALG13", "ZRANB3", "ANO5", "WDR11", "COG6", "MARK1", "PIK3CB", "MSH3", "DNAJB14", "GNPDA2", "DYNC1LI1", "DYNC1H1", "HOPX", "PAK1")

genes = c( "ZRANB3","COG6", "MARK1","PIK3CB", "MSH3", "DNAJB14","PDGFRA", "INTU", "PRSS55", "CCPG1", "CYRIA", "GNPDA2","PAK5", "SNX14", "PPP2R2C", "CACNA1B", "OLFM3")
VlnPlot(Brain_dev1, features = genes, pt.size = 0, stack = T,
        split.by = "cell_type", group.by = "species", flip = T, fill.by = "ident") +
  stat_summary(fun = "mean", geom = "point", 
               shape = 16, size = 2, color = "yellow3",
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = rev(rainbow(10))[3:9]) 
ggsave("../plots/Deve_Maca_Mouse_DEgene.Vlnplot.pdf", width = 8, height = 8)



DefaultAssay(Brain_dev) = "RNA"
FeaturePlot( subset(Brain_dev, species %in% c("macaque", "mouse")) , pt.size = 0.1,
             features = c("COG6", "MARK1"),  split.by = "species", cols = c("gray90", "red3")) & NoAxes()
ggsave("../plots/Deve_gene.featureplot.png", width = 14, height = 12)

Brain_dev1 = ScaleData(Brain_dev1, features = rownames(dev_DEG_p[dev_DEG_p$gouhui_DEG != "No", ]))

Brain_dev1$species_celltype <- paste(Brain_dev1$species, levels(Brain_dev1$cell_type), sep = " ")
Brain_dev1$species_celltype = factor(Brain_dev1$species_celltype, levels = c(paste("macaque", levels(Brain_dev1$cell_type), sep = " "),
                                                                             paste("mouse", levels(Brain_dev1$cell_type), sep = " ") ))
DoHeatmap(subset(Brain_dev1, downsample = 1000), features =  rownames(dev_DEG_p[dev_DEG_p$gouhui_DEG != "No", ]), 
          disp.min = -2, disp.max = 2, size = 4,
          group.by = c("species_celltype") ) &
  scale_fill_viridis_c()
ggsave("../plots/EQ_Gouhui_deve.heatmap.pdf", width = 10, height = 10)
ggsave("../plots/EQ_Gouhui_deve.heatmap.png", width = 10, height = 10)



Brain_dev2 = subset(Brain_dev1, cell_type %in% c("RG")) 

Idents(Brain_dev2) = "species"
RG_DEG = FindMarkers(Brain_dev2, ident.1 = "macaque", ident.2 = "mouse")
RG_DEG$gouhui_DEG = "No"
RG_DEG$gouhui_DEG[rownames(RG_DEG) %in% c(Up_DEG_gene, Deep_DEG_gene)] = "Gouhui"
RG_DEG$gouhui_DEG[rownames(RG_DEG) %in% expr_cor_ct_p[expr_cor_ct_p$cell %in% c("Upper-layer excitatory", "Deep-layer excitatory"), "gene"]] = "EQ"


ggplot(RG_DEG[RG_DEG$gouhui_DEG != "No", ], aes(avg_log2FC, -log10(p_val),color = gouhui_DEG)) +
  geom_point(size = 0.2)

RG_DEG_p = RG_DEG[RG_DEG$avg_log2FC > 0 & RG_DEG$p_val_adj < 0.05, ] # 152/348 ~ 43.7%


DotPlot(Brain_dev1, features = rownames(RG_DEG_p[RG_DEG_p$gouhui_DEG != "No", ])[1:10],
        split.by = "cell_type", group.by = "species") & coord_flip()


