library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# ==============================================================================
# Setup and Data Loading
# ==============================================================================

mdi_file <- "/data/xiaotx/Th/Figure/03_MDI/03_MDI_Proportion_Data.csv"
prop_df <- read.csv(mdi_file, stringsAsFactors = FALSE)

mam_file <- "/data/xiaotx/Th/Data/MamFuncDat.csv"
mam_df <- read.csv(mam_file, stringsAsFactors = FALSE, check.names = FALSE)

out_dir <- "/data/xiaotx/Th/Figure/08_Phenotype_Analysis_Final"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# Annotation and Merging
# ==============================================================================

prop_df$Species_Clean <- trimws(gsub("_", " ", prop_df$Species))

col_sci <- grep("Scientific", colnames(mam_df), ignore.case = TRUE, value = TRUE)[1]
col_noct <- grep("Activity-Nocturnal", colnames(mam_df), ignore.case = TRUE, value = TRUE)[1]
col_crep <- grep("Activity-Crepuscular", colnames(mam_df), ignore.case = TRUE, value = TRUE)[1]
col_diur <- grep("Activity-Diurnal", colnames(mam_df), ignore.case = TRUE, value = TRUE)[1]

mam_subset <- mam_df[, c(col_sci, col_noct, col_crep, col_diur)]
colnames(mam_subset) <- c("Scientific", "Act_Nocturnal", "Act_Crepuscular", "Act_Diurnal")
mam_subset$Scientific <- trimws(mam_subset$Scientific)

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

merged_df <- merge(prop_df, mam_subset, by.x = "Species_Mapped", by.y = "Scientific", all.x = TRUE)

# ==============================================================================
# Phenotype Classification (Gyrification & Activity)
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
        rep("Primates", 5),
        rep("Carnivora", 8),
        rep("Artiodactyla", 4),
        rep("Rodentia", 13),
        "Lagomorpha", rep("Eulipotyphla", 3), rep("Chiroptera", 3), "Diprotodontia"
    ),
    stringsAsFactors = FALSE
)

merged_df <- merge(merged_df, order_lookup, by = "Species_Mapped", all.x = TRUE)

merged_df <- merged_df %>%
    mutate(
        Gyrification = case_when(
            Order %in% c("Primates", "Carnivora", "Artiodactyla", "Cetartiodactyla") ~ 1,
            TRUE ~ 0
        ),
        Activity_Strict = case_when(
            Act_Diurnal == 1 & Act_Nocturnal == 0 & Act_Crepuscular == 0 ~ "Diurnal",
            Act_Nocturnal == 1 & Act_Diurnal == 0 & Act_Crepuscular == 0 ~ "Nocturnal",
            Act_Crepuscular == 1 & Act_Diurnal == 0 & Act_Nocturnal == 0 ~ "Crepuscular",
            TRUE ~ "Mixed/Unknown"
        )
    )

final_data <- merged_df %>%
    dplyr::select(
        Species, Species_Mapped, Order,
        Gyrification, Activity_Strict,
        Proportion
    )

write.csv(final_data, "/data/xiaotx/Th/Data/Phenotype_MamFuncDat_Matched.csv", row.names = FALSE)


# ==============================================================================
# Plotting
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

plot_phenotype <- function(data, x_col, filename_suffix) {
    groups <- levels(data[[x_col]])
    my_comparisons <- list(c(groups[1], groups[2]))

    stats <- data %>%
        group_by(!!sym(x_col)) %>%
        summarise(N = n(), Mean = mean(Proportion), SD = sd(Proportion), SE = SD / sqrt(N))
    write.csv(stats, file.path(out_dir, paste0("Stats_", filename_suffix, ".csv")), row.names = FALSE)

    p <- ggplot(data, aes_string(x = x_col, y = "Proportion")) +
        geom_boxplot(fill = "white", color = "black", width = 0.5, outlier.shape = NA) +
        geom_jitter(aes(color = Order), shape = 16, width = 0.2, size = 3, alpha = 0.8) +
        stat_compare_means(
            comparisons = my_comparisons,
            method = "wilcox.test",
            label = "p.signif",
            size = 8,
            vjust = 0.5,
            color = "black"
        ) +
        theme_bw() +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        labs(x = NULL, y = "MDI Proportion", title = NULL, subtitle = NULL, color = "Order") +
        theme(
            legend.position = "right",
            legend.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold", color = "black"),
            axis.text.y = element_text(size = 11),
            panel.grid.minor = element_blank()
        ) +
        scale_color_manual(values = order_colors)

    ggsave(file.path(out_dir, paste0("Boxplot_CustomColor_", filename_suffix, ".pdf")), p, width = 5.5, height = 5)

    return(p)
}

# 1. Gyrification Analysis
df_brain <- final_data %>%
    filter(!is.na(Gyrification)) %>%
    mutate(Gyrification_Label = factor(Gyrification, levels = c(0, 1), labels = c("Lissencephalic", "Gyrencephalic")))

plot_phenotype(df_brain, "Gyrification_Label", "Gyrification")

# 2. Activity Pattern Analysis
df_activity <- final_data %>%
    filter(Activity_Strict %in% c("Diurnal", "Nocturnal")) %>%
    mutate(Activity_Strict = factor(Activity_Strict, levels = c("Nocturnal", "Diurnal")))

plot_phenotype(df_activity, "Activity_Strict", "Strict_Activity")
