library(dplyr)
library(ggplot2)
library(tidyr)
library(pals)
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(purrr)
library(ggrepel)
library(ggpmisc)


metadata <- readRDS("/project/luyz/work/singlecell/Allmetadata_nohuman_1.8.RDS")
setwd("/project/luyz/work/singlecell/Paper_plots/Ratio/Oligo_ratio")

sub <- subset(metadata,Organ != "CE")

tab <- table(sub$Species,sub$cell_type)


cell_proportions <- prop.table(tab, margin = 1)
Oligo_Neuron <- cell_proportions[,"Oligodendrocyte"]/cell_proportions[,"Neuron"]
Oligo_all_Neuron <- (cell_proportions[,"Oligodendrocyte"]+cell_proportions[,"Oligodendrocyte precursor"])/cell_proportions[,"Neuron"]
Oligo <- cell_proportions[,"Oligodendrocyte"]+cell_proportions[,"Oligodendrocyte precursor"]
Glia_Neuron <- (cell_proportions[,"Oligodendrocyte"] + 
                cell_proportions[,"Oligodendrocyte precursor"] + 
                cell_proportions[,"Microglia"]) / 
               cell_proportions[,"Neuron"]
Glia <- cell_proportions[,"Oligodendrocyte"] + 
                cell_proportions[,"Oligodendrocyte precursor"] + 
                cell_proportions[,"Microglia"] +
                cell_proportions[,"Astrocyte"]
                
cell_proportions <- cbind(cell_proportions,Oligo_Neuron,Glia_Neuron,Oligo_all_Neuron,Oligo,Glia)


sample_info <- read.csv("/project/luyz/work/singlecell/Paper_plots/phenofile/phenofile.csv", stringsAsFactors = FALSE)
sample_info <- sample_info[, c("Species", "BW")]

proportions_df <- as.data.frame.matrix(cell_proportions)  
proportions_df$Species <- rownames(proportions_df)

merged_data <- proportions_df %>%
  left_join(sample_info, by = "Species") %>%
  mutate(
    log_BW = log10(BW),  
    Species = as.character(Species)
  )

summary(lm(Oligodendrocyte ~ log_BW,data = merged_data))

merged_data_long <- merged_data %>%
  pivot_longer(
    cols = -c(Species, BW, log_BW), 
    names_to = "CellType",
    values_to = "Proportion") %>%
  mutate(
    Species = as.character(Species),
    CellType = as.character(CellType)
  )


plot_data <- merged_data_long %>%
  filter(CellType == "Oligodendrocyte") %>%
  filter(!is.na(BW) & !is.na(Proportion))

p1 <- ggplot(plot_data, aes(x = log_BW, y = Proportion)) +
  geom_point(size = 8, alpha = 0.7, color = "#084b08") +
  geom_smooth(method = "lm", se = TRUE, color = "#f72092", linewidth = 0.5,fullrange = TRUE) +
    geom_text_repel(
    aes(label = Species),
    size = 10,                    
    color = "black",            
    point.padding = 0.5,         
    segment.color = "grey40",    
    segment.alpha = 0.8,         
    force = 10,                 
    max.iter = 10000,            
    max.overlaps = 50,           
    nudge_x = 0.05,              
    nudge_y = 0.05) +
  labs(
    x = expression(log10("Brain Weight")),
    y = "Oligo Proportion",
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.ticks = element_line(color = "black", linewidth = 0.5), 
    axis.title = element_text(face = "bold",color = "black",size = 8),
    axis.text.x = element_text(size = 8,color = "black"),  
    axis.text.y = element_text(size = 8,color = "black"),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8,color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 8,color = "black"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8, ,color = "black")
  )

ggsave("/project/luyz/work/singlecell/Paper_plots/Ratio/Oligo_ratio/Oligo_BW_dotplot.pdf",plot = p1 ,dpi = 600,width = 2,height = 2)


