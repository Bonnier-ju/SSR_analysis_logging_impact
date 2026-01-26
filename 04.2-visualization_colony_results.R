# ============================================================
# Script: 04.2-visualisation_colony_results.R
# Purpose: Visualize COLONY parentage outputs (BestConfig + BestCluster), map parent-offspring links, compute distances,
#          and relate DBH to reproductive success.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: <plot>.BestConfig.csv ; <plot>.BestCluster.csv ; <plot>_tree_info.csv (ID, lat, long, dbh, etc.)
# Outputs: figures/04_parentage_colony/<plot>/*.png ; results/04_parentage_colony/<plot>/*.csv
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
plot_name <- "PAI74"
base_path <- "."

# Inputs
bestconfig_file <- file.path(base_path, "results", "colony", plot_name, paste0(plot_name, "_colony_project.BestConfig.csv"))
bestcluster_file <- file.path(base_path, "results", "colony", plot_name, paste0(plot_name, "_colony_project.BestCluster.csv"))
treeinfo_file <- file.path(base_path, "data", plot_name, paste0(plot_name, "_tree_info.csv"))

# Outputs
fig_dir <- file.path("figures", "04_parentage_colony", plot_name)
res_dir <- file.path("results", "04_parentage_colony", plot_name)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

out_map_links <- file.path(fig_dir, paste0("colony_parent_offspring_map_", plot_name, ".png"))
out_map_clusters <- file.path(fig_dir, paste0("colony_family_clusters_map_", plot_name, ".png"))
out_hist_dist <- file.path(fig_dir, paste0("colony_parent_offspring_distance_hist_", plot_name, ".png"))
out_parent_offspring_csv <- file.path(res_dir, paste0("parent_offspring_counts_", plot_name, ".csv"))
out_distance_summary_csv <- file.path(res_dir, paste0("distance_summary_", plot_name, ".csv"))

# Cluster filter
min_cluster_size <- 10

# Map styling
offspring_fill <- "green3"
parent_fill <- "royalblue3"
link_color <- "azure4"
## ======================================================

## ===================== LIBRARIES ======================
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(geosphere)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(bestconfig_file))
stopifnot(file.exists(bestcluster_file))
stopifnot(file.exists(treeinfo_file))

best_config <- read.csv(bestconfig_file, header = TRUE, stringsAsFactors = FALSE)
best_cluster <- read.csv(bestcluster_file, header = TRUE, stringsAsFactors = FALSE)
tree_info <- read.csv(treeinfo_file, header = TRUE, stringsAsFactors = FALSE)

# Ensure unique IDs in metadata
tree_info <- tree_info %>% distinct(ID, .keep_all = TRUE)

# Basic checks
stopifnot(all(c("ID", "lat", "long") %in% names(tree_info)))
stopifnot(all(c("OffspringID", "FatherID", "MotherID") %in% names(best_config)))
stopifnot(all(c("OffspringID", "FatherID", "MotherID", "ClusterIndex") %in% names(best_cluster)))
## ======================================================

################################################################################
############################# 1) BESTCONFIG MAP ################################
################################################################################

# Offspring coords
offspring_positions <- best_config %>%
  left_join(tree_info, by = c("OffspringID" = "ID")) %>%
  rename(OffspringLat = lat, OffspringLong = long)

# Father coords
father_positions <- best_config %>%
  left_join(tree_info, by = c("FatherID" = "ID")) %>%
  rename(FatherLat = lat, FatherLong = long)

# Mother coords
mother_positions <- best_config %>%
  left_join(tree_info, by = c("MotherID" = "ID")) %>%
  rename(MotherLat = lat, MotherLong = long)

# Merge into one table
map_data <- offspring_positions %>%
  select(OffspringID, OffspringLat, OffspringLong, FatherID, MotherID) %>%
  bind_cols(father_positions %>% select(FatherLat, FatherLong)) %>%
  bind_cols(mother_positions %>% select(MotherLat, MotherLong))

# Identify sampled vs inferred parents (#...)
map_data <- map_data %>%
  mutate(
    has_real_father = !str_starts(FatherID, "#") & !is.na(FatherID) & FatherID != "",
    has_real_mother = !str_starts(MotherID, "#") & !is.na(MotherID) & MotherID != ""
  )

n_with_only_one_real_parent <- sum(xor(map_data$has_real_father, map_data$has_real_mother), na.rm = TRUE)
n_with_both_real_parents <- sum(map_data$has_real_father & map_data$has_real_mother, na.rm = TRUE)

cat("Offspring with ONE sampled parent:", n_with_only_one_real_parent, "\n")
cat("Offspring with TWO sampled parents:", n_with_both_real_parents, "\n")

# Plot links + points
p_links <- ggplot() +
  geom_segment(
    data = map_data %>% filter(has_real_father, !is.na(FatherLat)),
    aes(x = OffspringLong, y = OffspringLat, xend = FatherLong, yend = FatherLat),
    color = link_color, alpha = 0.9, linewidth = 0.8
  ) +
  geom_segment(
    data = map_data %>% filter(has_real_mother, !is.na(MotherLat)),
    aes(x = OffspringLong, y = OffspringLat, xend = MotherLong, yend = MotherLat),
    color = link_color, alpha = 0.9, linewidth = 0.8
  ) +
  geom_point(
    data = map_data,
    aes(x = OffspringLong, y = OffspringLat, fill = "Offspring"),
    shape = 21, size = 3, color = "black"
  ) +
  geom_point(
    data = map_data %>% filter(has_real_father, !is.na(FatherLat)),
    aes(x = FatherLong, y = FatherLat, fill = "Parent"),
    shape = 24, size = 3, color = "black"
  ) +
  geom_point(
    data = map_data %>% filter(has_real_mother, !is.na(MotherLat)),
    aes(x = MotherLong, y = MotherLat, fill = "Parent"),
    shape = 24, size = 3, color = "black"
  ) +
  scale_fill_manual(values = c("Offspring" = offspring_fill, "Parent" = parent_fill), name = "Type") +
  labs(title = plot_name, x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right")

print(p_links)
ggsave(out_map_links, p_links, width = 7, height = 7, dpi = 300, bg = "white")
cat("Saved:", out_map_links, "\n")

################################################################################
############################# 2) BESTCLUSTER MAP ###############################
################################################################################

# Replace inferred parent IDs by NA (only when the string starts with '#')
best_cluster <- best_cluster %>%
  mutate(
    FatherID = ifelse(str_detect(FatherID, "^#"), NA_character_, FatherID),
    MotherID = ifelse(str_detect(MotherID, "^#"), NA_character_, MotherID)
  )

# Filter clusters by size
valid_clusters <- best_cluster %>%
  count(ClusterIndex, name = "total_count") %>%
  filter(total_count >= min_cluster_size) %>%
  pull(ClusterIndex)

filtered_best_cluster <- best_cluster %>% filter(ClusterIndex %in% valid_clusters)

# Merge coords
offspring_positions_c <- filtered_best_cluster %>%
  left_join(tree_info, by = c("OffspringID" = "ID")) %>%
  rename(OffspringLat = lat, OffspringLong = long)

father_positions_c <- filtered_best_cluster %>%
  left_join(tree_info, by = c("FatherID" = "ID")) %>%
  rename(FatherLat = lat, FatherLong = long)

mother_positions_c <- filtered_best_cluster %>%
  left_join(tree_info, by = c("MotherID" = "ID")) %>%
  rename(MotherLat = lat, MotherLong = long)

cluster_data <- offspring_positions_c %>%
  bind_cols(father_positions_c %>% select(FatherLat, FatherLong)) %>%
  bind_cols(mother_positions_c %>% select(MotherLat, MotherLong))

# Cluster plot (same color for links + individuals)
n_clust <- length(unique(cluster_data$ClusterIndex))
cluster_cols <- scales::hue_pal()(n_clust)

p_clusters <- ggplot() +
  geom_segment(
    data = cluster_data %>% filter(!is.na(FatherLat)),
    aes(x = OffspringLong, y = OffspringLat, xend = FatherLong, yend = FatherLat, color = factor(ClusterIndex)),
    alpha = 0.6, linewidth = 0.8, show.legend = FALSE
  ) +
  geom_segment(
    data = cluster_data %>% filter(!is.na(MotherLat)),
    aes(x = OffspringLong, y = OffspringLat, xend = MotherLong, yend = MotherLat, color = factor(ClusterIndex)),
    alpha = 0.6, linewidth = 0.8, show.legend = FALSE
  ) +
  geom_point(
    data = cluster_data,
    aes(x = OffspringLong, y = OffspringLat, fill = factor(ClusterIndex)),
    shape = 23, size = 3, color = "black"
  ) +
  geom_point(
    data = cluster_data %>% filter(!is.na(FatherLat)),
    aes(x = FatherLong, y = FatherLat, fill = factor(ClusterIndex)),
    shape = 24, size = 3, color = "black"
  ) +
  geom_point(
    data = cluster_data %>% filter(!is.na(MotherLat)),
    aes(x = MotherLong, y = MotherLat, fill = factor(ClusterIndex)),
    shape = 24, size = 3, color = "black"
  ) +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  scale_color_manual(values = cluster_cols) +
  labs(
    title = paste("Clusters of Parent-Offspring relations (", plot_name, ")", sep = ""),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p_clusters)
ggsave(out_map_clusters, p_clusters, width = 8, height = 7, dpi = 300, bg = "white")
cat("Saved:", out_map_clusters, "\n")

################################################################################
############################# 3) DISTANCE STATS ################################
################################################################################

# Distances (use geosphere, vectorized)
map_data <- map_data %>%
  mutate(
    FatherDistance_m = ifelse(
      has_real_father & !is.na(FatherLat),
      distHaversine(cbind(OffspringLong, OffspringLat), cbind(FatherLong, FatherLat)),
      NA_real_
    ),
    MotherDistance_m = ifelse(
      has_real_mother & !is.na(MotherLat),
      distHaversine(cbind(OffspringLong, OffspringLat), cbind(MotherLong, MotherLat)),
      NA_real_
    ),
    ParentDistance_m = pmin(FatherDistance_m, MotherDistance_m, na.rm = TRUE)
  )

total_links <- sum(!is.na(map_data$FatherDistance_m)) + sum(!is.na(map_data$MotherDistance_m))

distance_summary <- map_data %>%
  summarise(
    MeanDistance_m = mean(ParentDistance_m, na.rm = TRUE),
    MinDistance_m  = min(ParentDistance_m, na.rm = TRUE),
    MaxDistance_m  = max(ParentDistance_m, na.rm = TRUE)
  ) %>%
  mutate(TotalLinks = total_links)

print(distance_summary)
write.csv(distance_summary, out_distance_summary_csv, row.names = FALSE)

# Optional histogram
p_hist <- ggplot(map_data, aes(x = ParentDistance_m)) +
  geom_histogram(binwidth = 10, fill = "grey50", color = "white") +
  labs(
    title = paste("Parent-offspring distances (", plot_name, ")", sep = ""),
    x = "Distance (m)",
    y = "Number of links"
  ) +
  theme_minimal()

print(p_hist)
ggsave(out_hist_dist, p_hist, width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved:", out_hist_dist, "\n")

################################################################################
############## 4) DBH vs OFFSPRING NUMBER (REPRODUCTIVE SUCCESS) ################
################################################################################

# Build parent -> offspring counts (father + mother pooled)
parent_offspring <- map_data %>%
  select(FatherID, MotherID) %>%
  pivot_longer(cols = c(FatherID, MotherID), names_to = "ParentType", values_to = "ParentID") %>%
  filter(!is.na(ParentID), ParentID != "", !str_starts(ParentID, "#")) %>%   # keep sampled parents only
  group_by(ParentID) %>%
  summarise(NumOffspring = n(), .groups = "drop") %>%
  left_join(tree_info, by = c("ParentID" = "ID"))

write.csv(parent_offspring, out_parent_offspring_csv, row.names = FALSE)
cat("Saved:", out_parent_offspring_csv, "\n")

# Linear model
model <- lm(NumOffspring ~ dbh, data = parent_offspring)
print(summary(model))

p_dbh <- ggplot(parent_offspring, aes(x = dbh, y = NumOffspring)) +
  geom_point(size = 3, alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  labs(
    title = paste(plot_name, "- DBH vs offspring number"),
    x = "DBH",
    y = "Number of offspring"
  ) +
  theme_minimal()

print(p_dbh)
ggsave(file.path(fig_dir, paste0("dbh_vs_offspring_", plot_name, ".png")),
       p_dbh, width = 7, height = 5, dpi = 300, bg = "white")
