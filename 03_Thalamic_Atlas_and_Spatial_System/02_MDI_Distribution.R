library(dplyr)
library(tidyr)
library(tibble)
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)
library(Seurat)

DIR_OBJ  <- "/data/xiaotx/Th/Object"
DIR_DATA <- "/data/xiaotx/Th/Data"
DIR_OUT  <- "/data/xiaotx/Th/Figure/03_MDI_Proportion"
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

rds_path <- file.path(DIR_OBJ, "Th_final_human.rds")
merged_plot <- readRDS(rds_path)

# ==============================================================================
# Calculate Neuron Proportions
# ==============================================================================

merged_plot$latin_name[merged_plot$latin_name == "Petaurus breviceps"] <- "Petaurus breviceps papuanus"
merged_plot$latin_name[merged_plot$latin_name == "Mustela putorius"] <- "Mustela putorius furo"
merged_plot$latin_name[merged_plot$orig.ident == "Human"] <- "Homo sapiens"

neuron_types <- c(
  "Thalamic excitatory",
  "Splatter",
  "Midbrain-derived inhibitory",
  "CGE interneuron",
  "MGE interneuron"
)

meta <- merged_plot@meta.data

total_neurons <- meta %>%
  filter(final.celltype %in% neuron_types) %>%
  group_by(latin_name) %>%
  summarise(total_neurons_in_species = n(), .groups = 'drop')

mdi_counts <- meta %>%
  filter(final.celltype == "Midbrain-derived inhibitory") %>%
  group_by(latin_name) %>%
  summarise(mdi_count = n(), .groups = 'drop')

tree_file_path <- file.path(DIR_DATA, "mammal.addhuman.final.tree")
phylo_tree <- read.tree(file = tree_file_path)
all_species_df <- data.frame(latin_name = gsub("_", " ", phylo_tree$tip.label))

all_species_proportions <- all_species_df %>%
  left_join(total_neurons, by = "latin_name") %>%
  left_join(mdi_counts, by = "latin_name") %>%
  mutate(
    total_neurons_in_species = ifelse(is.na(total_neurons_in_species), 0, total_neurons_in_species),
    mdi_count = ifelse(is.na(mdi_count), 0, mdi_count),
    proportion_of_neurons = ifelse(total_neurons_in_species == 0, 0, mdi_count / total_neurons_in_species)
  )

# ==============================================================================
# Data Export
# ==============================================================================

export_table <- all_species_proportions %>%
  dplyr::select(latin_name, mdi_count, total_neurons_in_species, proportion_of_neurons) %>%
  dplyr::rename(
    Species = latin_name,
    MDI_Count = mdi_count,
    Total_Neurons = total_neurons_in_species,
    Proportion = proportion_of_neurons
  ) %>%
  dplyr::arrange(desc(Proportion))

csv_output_path <- file.path(DIR_OUT, "03_MDI_Proportion_Data.csv")
write.csv(export_table, file = csv_output_path, row.names = FALSE, quote = FALSE)

# ==============================================================================
# Phylogenetic Tree Processing
# ==============================================================================

find_real_name <- function(target, all_tips) {
  if (target %in% all_tips) return(target)
  idx <- which(tolower(all_tips) == tolower(target))
  if (length(idx) > 0) return(all_tips[idx[1]])
  idx_grep <- grep(paste0("^", target), all_tips, ignore.case = TRUE)
  if (length(idx_grep) > 0) return(all_tips[idx_grep[1]])
  return(NA)
}

tips_to_drop <- c("Tupaia_chinensis", "Bos_grunniens")
real_drop <- tips_to_drop[tips_to_drop %in% phylo_tree$tip.label]

if (length(real_drop) > 0) {
  pruned_tree <- drop.tip(phylo_tree, tip = real_drop)
} else {
  pruned_tree <- phylo_tree
}

outgroup <- "Petaurus_breviceps_papuanus"
if (outgroup %in% pruned_tree$tip.label) {
  message("    正在以 ", outgroup, " 为外群定根...")
  pruned_tree <- root(pruned_tree, outgroup = outgroup, resolve.root = TRUE)
}

raw_manual_order <- c(
  "Homo_sapiens", "Macaca_mulatta", "Colobus_guereza", "Callithrix_jacchus", "Lemur_catta",
  "Mus_musculus", "Rattus_norvegicus", "Meriones_unguiculatus", "Phodopus_sungorus",
  "Graphiurus_kelleni", "Tamias_sibiricus", "Callosciurus_erythraeus", "Marmota_himalayana",
  "Cavia_porcellus", "Chinchilla_lanigera", "Heterocephalus_glaber", "Rhizomys_sinensis", "Hystrix_brachyura",
  "Oryctolagus_cuniculus",
  "Canis_lupus_familiaris", "Vulpes_lagopus", "Nyctereutes_procyonoides", "Mustela_putorius_furo",
  "Meles_meles", "Felis_catus", "Paguma_larvata", "Procyon_lotor",
  "Sus_scrofa", "Capra_hircus", "Bos_taurus", "Cervus_nippon",
  "Rousettus_leschenaultii", "Rhinolophus_blythi", "Hipposideros_armiger", "Suncus_murinus",
  "Scaptonyx_fusicaudus", "Atelerix_albiventris",
  "Petaurus_breviceps_papuanus"
)

current_tree_tips <- pruned_tree$tip.label
valid_constraints <- sapply(raw_manual_order, find_real_name, all_tips = current_tree_tips)
valid_constraints <- valid_constraints[!is.na(valid_constraints)]
valid_constraints <- unique(valid_constraints)

missing_tips <- setdiff(current_tree_tips, valid_constraints)
if (length(missing_tips) > 0) {
  valid_constraints <- c(valid_constraints, missing_tips)
}

final_constraints <- rev(valid_constraints)
pruned_tree_sorted <- rotateConstr(pruned_tree, constraint = final_constraints)

sorted_tips <- pruned_tree_sorted$tip.label

plot_data_final <- export_table %>%
  mutate(
    raw_name = gsub(" ", "_", Species),
    matched_name = sapply(raw_name, find_real_name, all_tips = phylo_tree$tip.label),
    label = matched_name,
    display_name = gsub("_", " ", label),
    display_name = gsub(" Gray", "", display_name),
    log_mdi_proportion = log1p(Proportion),
    final_label = paste0(display_name, " (", sprintf("%.2f%%", Proportion * 100), ")")
  ) %>%
  filter(!is.na(label)) %>%
  dplyr::select(label, log_mdi_proportion, final_label)

# ==============================================================================
# Tree Visualization
# ==============================================================================

p_tree <- ggtree(pruned_tree_sorted, ladderize = FALSE) %<+% plot_data_final +
 
  geom_tree(size = 0.8, color = "#4d4d4d") +

  geom_tippoint(
    aes(fill = log_mdi_proportion, size = log_mdi_proportion),
    shape = 21, color = "white", stroke = 0.6, alpha = 0.9
  ) +

  geom_tiplab(
    aes(label = final_label),
    size = 4, fontface = "italic", color = "black",
    offset = 0.2, align = TRUE, linesize = 0.2, linetype = "dotted"
  ) +

  scale_fill_viridis_c(
    option = "plasma", name = "log1p(Proportion)",
    guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 6)
  ) +
  
  scale_size_continuous(
    range = c(3, 10), name = "log1p(Proportion)",
    guide = guide_legend(title.position = "top", override.aes = list(fill = "black", color = "black"))
  ) +

  theme_tree2() +
  ggtree::hexpand(0.5) +
  theme(
    legend.position = c(0.1, 0.85),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key = element_blank()
  )

pdf_output_path <- file.path(DIR_OUT, "03_MDI_Proportion_Final.pdf")

ggsave(
  filename = pdf_output_path,
  plot = p_tree,
  device = "pdf",
  width = 14, height = 11, units = "in", dpi = 600
)
