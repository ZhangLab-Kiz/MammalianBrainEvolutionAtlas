library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(ggrepel)
library(ggraph)
library(tidygraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)
library(ggplot2)
library(dplyr)

# ==============================================================================
# Setup and Data Loading
# ==============================================================================

DIR_OUT <- "/data/xiaotx/Th/Figure/06_MDI_Gene_Strict"
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

rds_path <- "/data/xiaotx/Th/Object/th_joinlayers_neuron.rds"
seu <- readRDS(rds_path)

DefaultAssay(seu) <- "RNA"
Idents(seu) <- "final.celltype"

# ==============================================================================
# Data Preparation: Subset and Balance
# ==============================================================================

target_group <- "Midbrain-derived inhibitory"
ref_groups <- c("CGE interneuron", "MGE interneuron")

cell_stats <- table(seu$species, seu$final.celltype)

if (target_group %in% colnames(cell_stats)) {
    valid_species <- rownames(cell_stats)[cell_stats[, target_group] > 100]
} else {
    stop("check final.celltype")
}

target_types <- c(target_group, ref_groups)

seu_inhib <- subset(seu,
    idents = target_types,
    subset = species %in% valid_species
)

print(table(seu_inhib$final.celltype))

cells_mdi <- WhichCells(seu_inhib, idents = target_group)
cells_ref <- WhichCells(seu_inhib, idents = ref_groups)

n_mdi <- length(cells_mdi)
n_ref <- length(cells_ref)


if (n_mdi > n_ref) {

    set.seed(42)
    keep_mdi <- sample(cells_mdi, size = n_ref)

    final_cells <- c(keep_mdi, cells_ref)

    seu_balanced <- subset(seu_inhib, cells = final_cells)
} else {
    seu_balanced <- seu_inhib
}

print(table(seu_balanced$final.celltype))

# ==============================================================================
# MDI vs CGE+MGE
# ==============================================================================

markers_balanced <- FindMarkers(
    object = seu_balanced,
    ident.1 = target_group,
    ident.2 = ref_groups,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = FALSE
)

markers_balanced <- markers_balanced %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
        raw_logP = -log10(p_val_adj + 1e-300),
        logP_capped = pmin(raw_logP, 100),

        significance = case_when(
            p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up (MDI Specific)",
            p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down (Ref High)",
            TRUE ~ "Not Sig"
        )
    )

out_csv <- file.path(DIR_OUT, "DEG_MDI_vs_Inhibitory_Balanced.csv")
write.csv(markers_balanced, out_csv, row.names = FALSE)

saveRDS(markers_balanced, file.path(DIR_OUT, "DEG_MDI_vs_Inhibitory_Balanced.rds"))

rds_path <- file.path(DIR_OUT, "DEG_MDI_vs_Inhibitory_Balanced.rds")

markers_MDI <- readRDS(rds_path) %>%
    mutate(
        direction = case_when(
            "significance" %in% colnames(.) & significance == "Up (MDI Specific)" ~ "Up (MDI Specific)",
            "significance" %in% colnames(.) & significance == "Down (Ref High)" ~ "Down",
            avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up (MDI Specific)",
            avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
            TRUE ~ "No"
        )
    )

# ==============================================================================
# GO Enrichment Analysis
# ==============================================================================

genes_to_test <- markers_MDI %>%
    filter(direction == "Up (MDI Specific)") %>%
    pull(gene)

if (length(genes_to_test) > 0) {
    ego <- enrichGO(
        gene = genes_to_test,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05
    )

    if (!is.null(ego) && nrow(ego) > 0) {
        p_go <- dotplot(ego, showCategory = 15) +
            ggtitle("GO Enrichment: MDI Specific Genes (Strict)") +
            theme(axis.text.y = element_text(size = 10))

        ggsave(file.path(DIR_OUT, "GO_Enrichment_MDI_Strict.pdf"), p_go, width = 6, height = 8)
    } else {
        message("Can not.")
    }
}

ggsave(file.path(DIR_OUT, "GO_Enrichment.pdf"), p_go, width = 7, height = 8)

# ==============================================================================
# Volcano Plot
# ==============================================================================

plot_data <- markers_MDI

logfc_cut <- 0.25
pval_cut <- 0.05

plot_data$logP <- -log10(plot_data$p_val_adj + 1e-300)
max_p_val <- 100
plot_data$logP <- pmin(plot_data$logP, max_p_val)

plot_data$is_sig <- ifelse(plot_data$p_val_adj < pval_cut & plot_data$avg_log2FC > logfc_cut, "Yes", "No")

manual_genes <- c("OTX2")
top_genes <- plot_data %>%
    filter(direction == "Up (MDI Specific)") %>%
    top_n(n = 8, wt = avg_log2FC) %>%
    pull(gene)

label_genes <- unique(c(top_genes, manual_genes))

label_genes <- label_genes[label_genes %in% plot_data$gene]

p_vol <- ggplot(plot_data, aes(x = avg_log2FC, y = logP)) +
    geom_point(
        data = subset(plot_data, is_sig == "No"),
        color = "grey85", size = 1, alpha = 0.9
    ) +
    geom_point(
        data = subset(plot_data, is_sig == "Yes"),
        color = "#E41A1C", size = 2, alpha = 0.9
    ) +
    geom_text_repel(
        data = subset(plot_data, gene %in% label_genes),
        aes(label = gene),
        size = 3.5, fontface = "italic", color = "black",
        box.padding = 1.0, max.overlaps = 30, point.padding = 0.5
    ) +
    geom_vline(xintercept = c(logfc_cut), linetype = "dashed", color = "grey60") +
    theme_bw() +
    labs(
        title = "Volcano: MDI vs Inhibitory (Balanced)",
        x = expression(log[2] ~ Fold ~ Change),
        y = expression(-log[10](italic(P)[adj]))
    )

ggsave(file.path(DIR_OUT, "Volcano_MDI.pdf"), p_vol, width = 7, height = 6)


# ==============================================================================
# Cross-Species DotPlots
# ==============================================================================

target_genes <- label_genes

full_rds_path <- "/data/xiaotx/Th/Object/th_joinlayers_neuron.rds"
sub_seu <- readRDS(full_rds_path)

order_levels <- c(
    "Primates", "Rodentia", "Lagomorpha", "Carnivora",
    "Artiodactyla", "Chiroptera", "Eulipotyphla", "Marsupialia"
)
species_annotation <- data.frame(
    species = c(
        "Human", "Macaque", "Colobus", "Marmoset", "Lemur",
        "Mouse", "Rat", "Ger", "Cangshu",
        "Glis", "Tam", "Cal",
        "Marmot",
        "Cav", "Chin", "Heter", "Rhiz", "Hys",
        "Rabbit",
        "Dog", "Afox", "Rdog", "Ferret",
        "Meles", "Cat", "Pag", "Pro",
        "Pig", "Goat", "Cattle", "Deer",
        "Rous", "Rhin", "Hippo",
        "Shrew", "Scap", "Eri",
        "Pbv"
    ),
    Order = c(
        rep("Primates", 5),
        rep("Rodentia", 13),
        "Lagomorpha",
        rep("Carnivora", 8),
        rep("Artiodactyla", 4),
        rep("Chiroptera", 3),
        rep("Eulipotyphla", 3),
        "Marsupialia"
    )
)

sub_seu$Order <- species_annotation$Order[match(sub_seu$species, species_annotation$species)]
sub_seu$Order <- factor(sub_seu$Order, levels = order_levels)

sub_seu$Species_Group <- ifelse(sub_seu$Order == "Primates", "Primates", "Non-Primates")
sub_seu$Species_Group <- factor(sub_seu$Species_Group, levels = c("Primates", "Non-Primates"))

species_order_list <- unique(species_annotation$species[order(match(species_annotation$Order, order_levels))])
species_order_list <- intersect(species_order_list, unique(sub_seu$species))
sub_seu$species <- factor(sub_seu$species, levels = species_order_list)


neuron_types_order <- c(
    "Midbrain-derived inhibitory",
    "Thalamic excitatory",
    "Splatter",
    "CGE interneuron",
    "MGE interneuron"
)

for (gene in target_genes) {
    if (!gene %in% rownames(sub_seu)) next

    message(paste("Processing gene:", gene))

    tryCatch(
        {
            avg_exp <- AverageExpression(
                sub_seu,
                features = c(gene),
                group.by = c("final.celltype", "species")
            )$RNA

            avg_exp_long <- as.data.frame(as.table(as.matrix(avg_exp)))
            colnames(avg_exp_long) <- c("Gene", "Group", "Expression")

            avg_exp_long <- avg_exp_long %>%
                separate(Group, into = c("CellType", "Species"), sep = "_", extra = "merge")

            final_levels <- intersect(neuron_types_order, unique(avg_exp_long$CellType))

            if (length(final_levels) > 0) {
                avg_exp_long$CellType <- factor(avg_exp_long$CellType, levels = rev(final_levels))

                if (exists("species_order_list")) {
                    valid_sp <- intersect(species_order_list, unique(avg_exp_long$Species))
                    avg_exp_long$Species <- factor(avg_exp_long$Species, levels = valid_sp)
                }

                p_dot <- ggplot(avg_exp_long, aes(x = Species, y = CellType)) +
                    geom_point(aes(size = Expression, color = Expression)) +
                    scale_color_gradientn(
                        colours = c("white", "#FFF5F0", "#FC9272", "#CB181D", "#67000D"),
                        values = scales::rescale(c(0, 0.1, 1, 2, max(avg_exp_long$Expression, 4))),
                        name = "Avg Exp"
                    ) +
                    scale_size(range = c(0, 7), name = "Expression") +
                    theme_bw() +
                    labs(
                        title = paste0(gene, ": Specificity in MDI"),
                        x = "",
                        y = ""
                    ) +
                    theme(
                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", color = "black", size = 10),
                        axis.text.y = element_text(color = "black", face = "bold", size = 10),
                        panel.grid.major = element_line(colour = "grey92", linetype = "dashed"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1)
                    )

                out_filename <- file.path(DIR_OUT, paste0("DotPlot_", gene, ".pdf"))
                ggsave(out_filename, p_dot, width = 12, height = 4.5)
            } else {
                message(paste("Skipping", gene, "- no matching cell types found."))
            }
        },
        error = function(e) {
            message(paste("Error processing", gene, ":", e$message))
        }
    )
}

if (!exists("markers_balanced")) {
    markers_balanced <- readRDS(file.path(DIR_OUT, "DEG_MDI_vs_Inhibitory_Balanced.rds"))
}

up_genes <- markers_balanced %>%
    filter(significance == "Up (MDI Specific)") %>%
    pull(gene) %>%
    unique()

message("Gene number: ", length(up_genes))

dbs <- c(
    "GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "KEGG_2016", "KEGG_2019_Human", "Disease_Perturbations_from_GEO_up",
    "WikiPathways_2024_Human", "Panther_2016", "SynGO_2024", "Reactome_Pathways_2024",
    "GWAS_Catalog_2025", "UK_Biobank_GWAS_v1", "ClinVar_2025", "OMIM_Disease", "OMIM_Expanded"
)

setEnrichrSite("Enrichr")

enriched_results <- enrichr(genes = up_genes, databases = dbs)

enrich_out_dir <- file.path(DIR_OUT, "EnrichR_Results")
if (!dir.exists(enrich_out_dir)) dir.create(enrich_out_dir)

enrich_out_dir <- file.path(DIR_OUT, "EnrichR_Results_BossStyle")
if (!dir.exists(enrich_out_dir)) dir.create(enrich_out_dir)

for (db_name in names(enriched_results)) {
    res_df <- enriched_results[[db_name]]

    if (is.null(res_df) || nrow(res_df) == 0) {
        next
    }

    sig_check <- res_df %>% filter(P.value < 0.05)
    if (nrow(sig_check) == 0) {
        message("Skip: ", db_name)
        next
    }

    p <- plotEnrich(
        enriched_results[[db_name]],
        showTerms = 20,
        numChar = 50,
        y = "Ratio",
        orderBy = "P.value",
        title = paste0("Enrichment: ", db_name)
    )


    ggsave(file.path(enrich_out_dir, paste0(db_name, ".pdf")), p, width = 7, height = 8)

    write.csv(res_df, file.path(enrich_out_dir, paste0(db_name, ".csv")), row.names = FALSE)

    message("Saved: ", db_name)
}
