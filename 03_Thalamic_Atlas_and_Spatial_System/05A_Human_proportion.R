# ==============================================================================
# Calculate Proportion of IN1 in Human Thalamus
# ==============================================================================

seu <- readRDS("Human_Develop_Th/second-tri-thalamus_norm.rds")

target_clusters <- c("EN1", "EN2", "EN3", "IN1", "IN2", "IN3", "IN4", "IN5")

seu_neurons <- subset(seu, clusters %in% target_clusters)

prop_IN1 <- sum(seu_neurons$clusters == "IN1") / length(seu_neurons$clusters)

message(paste0("IN1 Proportion: ", round(prop_IN1 * 100, 2), "%"))
