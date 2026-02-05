# ==============================================================================
# Step 04: Species Composition and L4 Characterization
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggtree)
library(tibble)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)

DIR_OBJ <- "./Query/Query_object"
DIR_FIG <- "/home/xiaotx/WholeBrain/Query/Query_Figure/Paper"
if (!dir.exists(DIR_FIG)) dir.create(DIR_FIG, recursive = TRUE)

# ==============================================================================
# Load Data and Preprocessing
# ==============================================================================

if (!exists("obj_filtered")) {
  cat("Loading integrated object...\n")
  obj_filtered <- readRDS(file.path(DIR_OBJ, "All_Query_Mapped_Final_0128.rds"))
}

order_levels <- c("Primates", "Rodentia", "Lagomorpha", "Carnivora", "Artiodactyla", "Eulipotyphla")
species_annotation <- data.frame(
  orig.ident = c(
    "Human", "Macaque", "Colobus", "Marmoset", "Lemur",
    "Mouse", "Rat", "Cangshu", "Ger", "Chin", "Cav", "Marmot", "Heter", "Rabbit",
    "Cat", "Dog", "Afox", "Rdog", "Ferret", "Meles",
    "Cattle", "Pig", "Goat", "Hippo", "Shrew"
  ),
  Order = c(
    rep("Primates", 5), rep("Rodentia", 8), "Lagomorpha",
    rep("Carnivora", 6), rep("Artiodactyla", 4), "Eulipotyphla"
  )
)

obj_filtered$Order <- species_annotation$Order[match(obj_filtered$orig.ident, species_annotation$orig.ident)]
obj_filtered$Order <- factor(obj_filtered$Order, levels = order_levels)

obj_filtered$Species_Group <- ifelse(obj_filtered$Order == "Primates", "Primates", "Non-Primates")
obj_filtered$Species_Group <- factor(obj_filtered$Species_Group, levels = c("Primates", "Non-Primates"))

species_order_list <- unique(species_annotation$orig.ident[order(match(species_annotation$Order, order_levels))])
species_order_list <- intersect(species_order_list, unique(obj_filtered$orig.ident))
obj_filtered$orig.ident <- factor(obj_filtered$orig.ident, levels = species_order_list)

# ==============================================================================
# Global UMAP Visualization
# ==============================================================================

my_colors <- c(
  "L2/3 IT"      = "#E41A1C",
  "L4 IT"        = "#377EB8",
  "L5 IT"        = "#4DAF4A",
  "L6 IT"        = "#0b290f",
  "L5 ET"        = "#984EA3",
  "L5/6 NP"      = "#FF7F00",
  "L6 CT"        = "#008080",
  "L6b"          = "#A65628",
  "L6 IT Car3"   = "#F781BF"
)

p_pure_cell <- DimPlot(obj_filtered,
  reduction = "ref.umap",
  group.by = "final.celltype",
  cols = my_colors,
  label = FALSE,
  raster = FALSE
) +
  NoAxes() +
  NoLegend() +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

output_filepath_cell_pure <- file.path(DIR_FIG, "Fig4_A_UMAP_CellType_Pure.png")
ggsave(
  filename = output_filepath_cell_pure,
  plot = p_pure_cell,
  device = "png",
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

p_for_legend_cell <- DimPlot(obj_filtered,
  reduction = "ref.umap",
  group.by = "final.celltype",
  cols = my_colors
) + theme(legend.position = "right")

g_cell <- ggplot2::ggplotGrob(p_for_legend_cell)
legend_index_cell <- which(sapply(g_cell$grobs, function(x) x$name) == "guide-box")

if (length(legend_index_cell) > 0) {
  legend_grob_cell <- g_cell$grobs[[legend_index_cell[1]]]

  output_filepath_legend_cell <- file.path(DIR_FIG, "Fig4_A_UMAP_CellType_Legend.pdf")
  pdf(output_filepath_legend_cell, width = 4, height = 6)
  grid::grid.newpage()
  grid::grid.draw(legend_grob_cell)
  dev.off()
  cat(sprintf("Saved CellType Legend to: %s\n", output_filepath_legend_cell))
}

# ==============================================================================
# L4 Proportion Analysis
# ==============================================================================

total_counts <- obj_filtered@meta.data %>%
  group_by(orig.ident, Order) %>%
  summarise(Total_Cells = n(), .groups = "drop")

l4_counts <- obj_filtered@meta.data %>%
  filter(final.celltype == "L4 IT") %>%
  group_by(orig.ident) %>%
  summarise(L4_Count = n(), .groups = "drop")

plot_data_l4 <- left_join(total_counts, l4_counts, by = "orig.ident") %>%
  mutate(L4_Count = ifelse(is.na(L4_Count), 0, L4_Count)) %>%
  mutate(Percentage = L4_Count / Total_Cells * 100)

plot_data_l4$Highlight <- ifelse(plot_data_l4$Order == "Primates", "Primates", "Non-Primates")
plot_data_l4$Highlight <- factor(plot_data_l4$Highlight, levels = c("Primates", "Non-Primates"))

table_export <- plot_data_l4 %>%
  rename(
    Species = orig.ident,
    Taxonomy_Order = Order,
    Total_Neurons = Total_Cells,
    L4_IT_Neurons = L4_Count,
    L4_Proportion_Percent = Percentage
  ) %>%
  mutate(L4_Proportion_Percent = round(L4_Proportion_Percent, 2)) %>%
  select(Species, Taxonomy_Order, Total_Neurons, L4_IT_Neurons, L4_Proportion_Percent)

write.csv(table_export, file.path(DIR_FIG, "Table_S1_L4_Proportion_Stats.csv"), row.names = FALSE)

p_bar_l4 <- ggplot(plot_data_l4, aes(x = orig.ident, y = Percentage)) +
  geom_bar(aes(fill = Highlight), stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("Primates" = "#D73027", "Non-Primates" = "grey80")) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_data_l4$Percentage) * 1.1)) + # 留点头部空间
  labs(
    title = "Proportion of L4 IT Neurons Across Species",
    y = "Proportion (% of Total Cells)",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(DIR_FIG, "Fig4_B_L4_Proportion_Only.pdf"), p_bar_l4, width = 8, height = 4)

# ==============================================================================
# FOXP2 Expression Analysis
# ==============================================================================

df_foxp2 <- FetchData(obj_filtered, vars = c("FOXP2", "final.celltype", "Order"))

df_foxp2$final.celltype <- as.character(df_foxp2$final.celltype)
df_foxp2$final.celltype[df_foxp2$final.celltype %in% c("L4 IT", "L5 IT")] <- "L4/5 IT"

celltype_levels <- c(
  "L2/3 IT",
  "L4/5 IT",
  "L5 ET",
  "L5/6 NP",
  "L6 IT",
  "L6 CT",
  "L6b",
  "L6 IT Car3"
)

df_foxp2 <- df_foxp2 %>% filter(final.celltype %in% celltype_levels)
df_foxp2$final.celltype <- factor(df_foxp2$final.celltype, levels = celltype_levels)

order_levels <- c("Primates", "Rodentia", "Lagomorpha", "Carnivora", "Artiodactyla", "Eulipotyphla")
df_foxp2$Order <- factor(df_foxp2$Order, levels = order_levels)

order_colors <- c(
  "Primates"     = "#E41A1C",
  "Rodentia"     = "#377EB8",
  "Lagomorpha"   = "#4DAF4A",
  "Carnivora"    = "#984EA3",
  "Artiodactyla" = "#FF7F00",
  "Eulipotyphla" = "#A65628"
)

p_stack_wide <- ggplot(df_foxp2, aes(x = final.celltype, y = FOXP2, fill = Order)) +
  geom_violin(scale = "width", adjust = 1.2, trim = TRUE, linewidth = 0.2) +
  facet_grid(Order ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = order_colors) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() +
  labs(
    title = "FOXP2 Expression Landscape Across Orders",
    y = "Expression Level",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12, color = "black"),
    strip.background = element_rect(fill = "transparent", color = NA),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "none"
  )

output_file <- file.path(DIR_FIG, "Fig4_D_Violin.pdf")
ggsave(output_file, p_stack_wide, width = 8, height = 5)

# ==============================================================================
# Split UMAP Visualization (Primates vs Non-Primates)
# ==============================================================================

target_celltype <- "L4 IT"
highlight_color <- "#D73027"
bg_color <- "grey90"

obj_filtered$plot_group <- "Background"

obj_filtered$plot_group[
  obj_filtered$Species_Group == "Primates" &
    obj_filtered$final.celltype == target_celltype
] <- "Pri_Target"

obj_filtered$plot_group[
  obj_filtered$Species_Group == "Non-Primates" &
    obj_filtered$final.celltype == target_celltype
] <- "NonPri_Target"

obj_pri <- subset(obj_filtered, subset = Species_Group == "Primates")

obj_pri$plot_group <- ifelse(
  obj_pri$final.celltype == target_celltype,
  "Target",
  "Background"
)

p_pure_pri <- DimPlot(obj_pri,
  reduction = "ref.umap",
  group.by = "plot_group",
  cols = c(
    "Target"     = highlight_color,
    "Background" = bg_color
  ),
  order = "Target",
  label = FALSE,
  raster = FALSE
) +
  NoAxes() +
  NoLegend() +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(
  filename = file.path(DIR_FIG, "Fig4_A_UMAP_L4_Highlight_Primates_Split.png"),
  plot = p_pure_pri,
  device = "png",
  width = 8, height = 8, units = "in", dpi = 600
)

obj_non <- subset(obj_filtered, subset = Species_Group == "Non-Primates")

obj_non$plot_group <- ifelse(
  obj_non$final.celltype == target_celltype,
  "Target",
  "Background"
)

p_pure_non <- DimPlot(obj_non,
  reduction = "ref.umap",
  group.by = "plot_group",
  cols = c(
    "Target"     = highlight_color,
    "Background" = bg_color
  ),
  order = "Target",
  label = FALSE,
  raster = FALSE
) +
  NoAxes() +
  NoLegend() +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(
  filename = file.path(DIR_FIG, "Fig4_A_UMAP_L4_Highlight_NonPrimates_Split.png"),
  plot = p_pure_non,
  device = "png",
  width = 8, height = 8, units = "in", dpi = 600
)

# ==============================================================================
# Additional Key Candidates Validation
# ==============================================================================

target_genes <- c("RORB", "FOXP2", "IL1RAPL2")

layer_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L5 ET", "L5/6 NP", "L6 CT", "L6b", "L6 IT Car3")

for (gene in target_genes) {
  message(paste("Processing gene:", gene))

  tryCatch(
    {
      avg_exp <- AverageExpression(
        obj_filtered,
        features = c(gene),
        group.by = c("final.celltype", "orig.ident")
      )$RNA

      avg_exp_long <- as.data.frame(as.table(as.matrix(avg_exp)))
      colnames(avg_exp_long) <- c("Gene", "Group", "Expression")

      avg_exp_long <- avg_exp_long %>%
        separate(Group, into = c("CellType", "Species"), sep = "_", extra = "merge")

      final_levels <- intersect(layer_order, unique(avg_exp_long$CellType))

      if (length(final_levels) > 0) {
        avg_exp_long$CellType <- factor(avg_exp_long$CellType, levels = rev(final_levels))

        if (exists("species_order_list")) {
          valid_species <- intersect(species_order_list, unique(avg_exp_long$Species))
          avg_exp_long$Species <- factor(avg_exp_long$Species, levels = valid_species)
        }

        p_dot <- ggplot(avg_exp_long, aes(x = Species, y = CellType)) +
          geom_point(aes(size = Expression, color = Expression)) +
          scale_color_gradientn(
            colours = c("white", "#FFF5F0", "#FC9272", "#CB181D", "#67000D"),
            values = scales::rescale(c(0, 0.1, 1, 2, 4)),
            name = "Avg Exp"
          ) +
          scale_size(range = c(0, 7), name = "Expression") +
          theme_bw() +
          labs(
            title = paste0(gene, ": Specificity in Primate L4"),
            x = "",
            y = ""
          ) +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", color = "black", size = 10),
            axis.text.y = element_text(color = "black", face = "bold", size = 10),
            panel.grid.major = element_line(colour = "grey92", linetype = "dashed"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1)
          )

        out_filename <- file.path(DIR_FIG, paste0("Fig1_C_", gene, "_DotPlot_WhiteRed.pdf"))
        ggsave(out_filename, p_dot, width = 12, height = 4.5)
      } else {
        warning(paste("No valid cell types found for gene:", gene))
      }
    },
    error = function(e) {
      message(paste("Error processing gene", gene, ":", e$message))
    }
  )
}
