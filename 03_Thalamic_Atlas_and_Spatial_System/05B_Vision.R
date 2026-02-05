library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(scales)

# ==============================================================================
# Setup and Data Loading
# ==============================================================================

mdi_file <- "/data/xiaotx/Th/Figure/03_MDI/03_MDI_Proportion_Data.csv"
prop_df <- read.csv(mdi_file, stringsAsFactors = FALSE)

db_file <- "/data/xiaotx/Th/Data/acuityblind.csv"
db_df <- read.csv(db_file, stringsAsFactors = FALSE, check.names = FALSE)

out_dir <- "/data/xiaotx/Th/Figure/08_Phenotype_Analysis_Final"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# Data Cleaning and Merging
# ==============================================================================

prop_df$Species_Clean <- trimws(gsub("_", " ", prop_df$Species))

name_corrections <- c(
    "Mustela putorius furo"       = "Mustela putorius",
    "Canis lupus familiaris"      = "Canis lupus",
    "Petaurus breviceps papuanus" = "Petaurus breviceps",
    "Rhinolophus blythi"          = "Rhinolophus pusillus"
)

prop_df$Species_Mapped <- ifelse(prop_df$Species_Clean %in% names(name_corrections),
    name_corrections[prop_df$Species_Clean],
    prop_df$Species_Clean
)

if ("MSW05_Binomial" %in% colnames(db_df)) {
    db_species_col <- "MSW05_Binomial"
} else {
    db_species_col <- grep("Binomial", colnames(db_df), ignore.case = TRUE, value = TRUE)[1]
}

db_df$Species_Mapped <- trimws(gsub("_", " ", db_df[[db_species_col]]))

cols_to_keep <- c("Species_Mapped", "VA")
db_subset <- db_df[, cols_to_keep, drop = FALSE]
db_subset[db_subset == -999] <- NA

final_df <- merge(prop_df, db_subset, by = "Species_Mapped", all.x = TRUE)

# ==============================================================================
# Add Taxonomic Order Information
# ==============================================================================

order_lookup <- data.frame(
    Species_Mapped = c(
        "Callithrix jacchus", "Colobus guereza", "Homo sapiens", "Lemur catta", "Macaca mulatta",
        "Canis lupus", "Felis catus", "Meles meles", "Mustela putorius",
        "Nyctereutes procyonoides", "Paguma larvata", "Procyon lotor", "Vulpes lagopus",
        "Bos taurus", "Capra hircus", "Cervus nippon", "Sus scrofa",
        "Callosciurus erythraeus", "Cavia porcellus", "Chinchilla lanigera", "Graphiurus kelleni",
        "Heterocephalus glaber", "Hystrix brachyura", "Marmota himalayana", "Meriones unguiculatus",
        "Mus musculus", "Phodopus sungorus", "Rattus norvegicus", "Rhizomys sinensis", "Tamias sibiricus",
        "Oryctolagus cuniculus", "Atelerix albiventris", "Scaptonyx fusicaudus", "Suncus murinus",
        "Hipposideros armiger", "Rhinolophus pusillus", "Rousettus leschenaultii", "Petaurus breviceps"
    ),
    Order = c(
        rep("Primates", 5), rep("Carnivora", 8), rep("Artiodactyla", 4), rep("Rodentia", 13),
        "Lagomorpha", rep("Eulipotyphla", 3), rep("Chiroptera", 3), "Diprotodontia"
    ),
    stringsAsFactors = FALSE
)

final_df <- merge(final_df, order_lookup, by = "Species_Mapped", all.x = TRUE)

plot_data <- final_df %>% filter(!is.na(VA))
cat(">>> 最终用于绘图的物种数量:", nrow(plot_data), "\n")

write.csv(plot_data, file.path(out_dir, "Data_With_VA_Cleaned.csv"), row.names = FALSE)

# ==============================================================================
# Boxplot: High vs Low Visual Acuity
# ==============================================================================

order_colors <- c(
    "Artiodactyla"  = "#E41A1C",
    "Carnivora"     = "#377EB8",
    "Chiroptera"    = "#4DAF4A",
    "Diprotodontia" = "#0b290f",
    "Eulipotyphla"  = "#984EA3",
    "Lagomorpha"    = "#FF7F00",
    "Primates"      = "#008080",
    "Rodentia"      = "#A65628"
)

target_labels <- c(
    "Homo sapiens",
    "Mus musculus",
    "Sus scrofa",
    "Canis lupus",
    "Oryctolagus cuniculus"
)

plot_data$Label_Show <- ifelse(plot_data$Species_Mapped %in% target_labels,
    plot_data$Species_Mapped,
    NA
)

threshold_va <- 3

plot_data$Visual_Group <- ifelse(plot_data$VA > threshold_va,
    "High Acuity (> 3 cpd)",
    "Low Acuity (<= 3 cpd)"
)

plot_data$Visual_Group <- factor(plot_data$Visual_Group,
    levels = c("Low Acuity (<= 3 cpd)", "High Acuity (> 3 cpd)")
)

my_comparisons <- list(c("Low Acuity (<= 3 cpd)", "High Acuity (> 3 cpd)"))

p_box <- ggplot(plot_data, aes(x = Visual_Group, y = Proportion)) +
    geom_boxplot(fill = "white", color = "black", width = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = Order), shape = 16, width = 0.2, size = 3, alpha = 0.8) +
    stat_compare_means(
        comparisons = my_comparisons, method = "wilcox.test",
        label = "p.signif", size = 8, vjust = 0.5, color = "black"
    ) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = "MDI Proportion", color = "Order") +
    theme(
        legend.position = "right",
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11),
        panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = order_colors)

ggsave(file.path(out_dir, "Boxplot_VA_Threshold_Group.pdf"), p_box, width = 5.5, height = 5)

# ==============================================================================
# Scatter Plot: Visual Acuity (Log Scale) vs MDI Proportion
# ==============================================================================

threshold_val <- 3

p_scatter_log <- ggplot(plot_data, aes(x = VA, y = Proportion)) +
    geom_vline(xintercept = threshold_val, linetype = "longdash", color = "grey40", size = 0.8, alpha = 0.8) +
    annotate("text",
        x = threshold_val,
        y = max(plot_data$Proportion, na.rm = TRUE) * 0.95,
        label = "Threshold (3 cpd)",
        angle = 90,
        vjust = -0.8,
        hjust = 1,
        size = 4,
        color = "grey40",
        fontface = "italic"
    ) +
    geom_smooth(method = "lm", color = "black", fill = "lightgray", alpha = 0.2) +
    geom_point(aes(color = Order), shape = 16, size = 3.5, alpha = 0.8) +
    stat_cor(
        method = "pearson", label.x.npc = "left", label.y.npc = "top",
        size = 5, p.accuracy = 0.001, r.accuracy = 0.01
    ) +
    geom_text_repel(
        aes(label = Label_Show),
        size = 4,
        fontface = "italic",
        min.segment.length = 0,
        box.padding = 0.8,
        max.overlaps = Inf
    ) +
    scale_x_log10(
        breaks = c(0.1, 0.3, 1, 3, 10, 30, 60),
        labels = c("0.1", "0.3", "1", "3", "10", "30", "60")
    ) +
    annotation_logticks(sides = "b") +
    scale_color_manual(values = order_colors) +
    theme_bw() +
    labs(
        x = "Visual Acuity (Cycles Per Degree, Log Scale)",
        y = "MDI Proportion",
        color = "Order"
    ) +
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right"
    )

out_dir <- "/data/xiaotx/Th/Figure/08_Phenotype_Analysis_Final"
ggsave(file.path(out_dir, "Scatter_VA_LogScale_WithThreshold.pdf"), p_scatter_log, width = 7.5, height = 6)