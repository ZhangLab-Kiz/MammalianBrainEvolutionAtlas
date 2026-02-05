# species-wise correlation network

library(tidyverse)
library(igraph)

avg_expr = lapply(names(avg_expr_list)[-c(7,22)], function(x){
  res = avg_expr_list[[x]]
  rownames(res) = paste(x, rownames(res), sep = "|")
  res
})
avg_expr = do.call(rbind, avg_expr)
avg_expr = t(avg_expr)

avg_expr.cor = cor(avg_expr)

nodes_df <- tibble(name = colnames(avg_expr.cor)) %>%
  separate(name, into = c("cell", "species"), sep = "\\|", remove = FALSE) %>%
  distinct()

edge_df <- nodes_df %>%
  dplyr::inner_join(nodes_df, by = "species", suffix = c("_i", "_j")) %>%
  dplyr::filter(cell_i != cell_j) %>%
  mutate(cor = avg_expr.cor[cbind(match(name_i, rownames(avg_expr.cor)),
                                  match(name_j, colnames(avg_expr.cor)))],
         length = 1 - cor) %>%
  transmute(from  = name_i,
            to    = name_j,
            cor,
            length,
            species,
            cell)

edge_df = reshape2::melt(avg_expr.cor)
edge_df$length = 1-edge_df$value
colnames(edge_df) = c("from", "to", "cor", "length")
edge_df$from_spe = str_remove_all(edge_df$from,".*\\|")
edge_df$from_cell = str_remove_all( edge_df$from, "\\|.*")
edge_df$to_spe = str_remove_all(edge_df$to,".*\\|")
edge_df$to_cell = str_remove_all( edge_df$to, "\\|.*")
edge_df = edge_df[edge_df$from_spe != edge_df$to_spe,]
edge_df = edge_df[ (edge_df$from_cell != edge_df$to_cell & edge_df$cor > 0.62) | (edge_df$from_cell == edge_df$to_cell & edge_df$cor > 0.45),]


edge_df <- edge_df %>%
  mutate(
    pair_id = pmap_chr(list(from, to), ~ paste(sort(c(..1, ..2)), collapse = "|"))
  ) %>%
  group_by(pair_id) %>%
  slice(1) %>%  
  ungroup() %>%
  dplyr::select(-pair_id)

g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = nodes_df)

set.seed(123)
E(g)$length <- edge_df$length
E(g)$weight <- edge_df$cor

layout_coords <- layout_with_kk(g, weights = E(g)$distance)

g_tbl <- as_tbl_graph(g) %>% 
  activate(nodes) %>% 
  mutate(x = layout_coords[, 1], y = layout_coords[, 2])

ggraph(g_tbl, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(alpha = cor, width = cor), colour = "grey50") +
  scale_edge_alpha_continuous(range = c(0.1, 1)) +
  scale_edge_width_continuous(range = c(0.1, 1.5)) +
  geom_node_point(aes(color = cell), size = 2) +
  scale_color_manual(values = ct_col_new) +
  coord_equal() +
  theme_graph() +
  labs(title = "Intra-species cell-type similarity network")


