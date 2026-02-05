# plot gene-term netwrok
library(tidygraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)
library(ggforce) 

###
# tidy_gene_term is a data.frame including enriched Terms and Genes.
# Term                     gene
# INTELLECTUAL DISABILITY CYFIP2
# INTELLECTUAL DISABILITY TRIO  
# INTELLECTUAL DISABILITY MTR  
#

colnames(tidy_gene_term) = c("fun", "gene")
edges = tidy_gene_term
g  <- as_tbl_graph(edges, directed = FALSE)
ig <- as.igraph(g)

V(ig)$type <- ifelse(V(ig)$name %in% edges$fun, "function", "gene")

fun_id  <- V(ig)[type == "function"]
gene_id <- V(ig)[type == "gene"]

lay_inner <- layout_in_circle(induced_subgraph(ig, fun_id))      
lay_outer <- layout_in_circle(induced_subgraph(ig, gene_id)) * 2.2  

lay_df <- tibble(
  name = V(ig)$name,
  type = V(ig)$type,
  x    = c(lay_inner[,1], lay_outer[,1]),
  y    = c(lay_inner[,2], lay_outer[,2])
)

fun_long <- edges %>% distinct(gene, fun)   
lay_df <- lay_df %>%
  left_join(fun_long %>% rename(name = gene), by = "name") 

pie_data <- edges %>%
  distinct(gene, fun) %>%                                  
  left_join(lay_df %>% filter(type == "gene") %>% select(name, x, y), by = c("gene" = "name")) %>%
  group_by(gene) %>%
  mutate(n   = n(),
         ang = 360 / n,
         id  = row_number(),
         start = lag(cumsum(ang), default = 0)) %>%
  ungroup() %>%
  mutate(r = 0.08)   # 外环扇形半径

label_data <- pie_data %>%
  mutate(
    mid_angle_rad = (start + ang/2) * pi / 180,
    label_x = x + (r + 0.05) * cos(mid_angle_rad),
    label_y = y + (r + 0.05) * sin(mid_angle_rad),
    angle_deg = (start + ang/2) * pi / 180 * 180/pi,
    angle_adj = ifelse(angle_deg > 90 & angle_deg < 270, 
                       angle_deg + 180,
                       angle_deg)
  )

segment_data <- edges %>%
  left_join(
    lay_df %>% 
      filter(type == "function") %>% 
      select(fun = name, x_fun = x, y_fun = y),
    by = "fun"
  ) %>%
  left_join(
    lay_df %>% 
      filter(type == "gene") %>% 
      select(gene = name, x_gene = x, y_gene = y),
    by = "gene"
  ) %>%
  distinct()

all_funs <- unique(edges$fun)

library(RColorBrewer)
color_palette <- brewer.pal(length(all_funs), "Dark2")

fun_colors <- setNames(color_palette[1:length(all_funs)], all_funs)

pie_data <- pie_data %>%
  mutate(fun_color = fun_colors[fun])

lay_df_fun <- lay_df %>%
  filter(type == "function") %>%
  mutate(fun_color = fun_colors[name])

p <- ggplot() +
  geom_curve(
    data = segment_data,
    aes(x = x_fun, y = y_fun, xend = x_gene, yend = y_gene),
    curvature = 0.25,
    angle = 90,
    ncp = 20,
    colour = "grey80", 
    size = 0.5, 
    alpha = 1
  ) +
  
  geom_point(
    data = lay_df_fun,
    aes(x, y, colour = name),
    size = 15, stroke = 1.5
  ) +
  
  geom_text(
    data = lay_df_fun,
    aes(x, y, label = name),
    size = 4, fontface = "bold", vjust = 0.5, hjust = 0.5
  ) +
  
  ggforce::geom_arc_bar(
    data = pie_data,
    aes(x0 = x, y0 = y, r = r, r0 = 0,
        start = start * pi / 180, end = (start + ang) * pi / 180,
        fill = fun),
    colour = NA, alpha = .9
  ) +
  
  geom_text(
    data = lay_df %>% filter(type == "gene"),
    aes(x = x * 1.1, y = y * 1.1, label = name),
    size = 3.5, vjust = 0.5, hjust = 0.5
  ) +
  
  scale_colour_manual(
    name = "功能模块",
    values = fun_colors
  ) +
  
  scale_fill_manual(
    name = "功能模块",
    values = fun_colors
  ) +
  
  coord_fixed(ratio = 1, xlim = c(-2.4, 2.4), ylim = c(-2.4, 2.4)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

print(p)

ggsave("../plots/EQ_Neuron_functions.netplot.pdf", width = 7, height = 7)



