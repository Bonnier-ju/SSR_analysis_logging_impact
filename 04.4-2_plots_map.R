########################################################################
######################### 2 plots map ##################################
########################################################################

library(readr)
library(dplyr)
library(stringr)


path_filtered <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo/filtered_data_PAI74_HKO50.csv"
path_cluster  <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.2-parentage_Colony/PAI74_HKO50/results_colony_HKO50_PAI74/Colony_HKO50_PAI74.BestCluster.csv"
output_path   <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis/df_map.csv"

filtered <- read_csv(path_filtered, show_col_types = FALSE)
clusters <- read_csv(path_cluster,  show_col_types = FALSE)


# Derive Plot/DBH for Mother and Father 
filtered_aug <- filtered %>%
  mutate(
    Mother_ID = as.character(Mother_ID),
    Father_ID = as.character(Father_ID),
    Parent1_Cervus = as.character(Parent1_Cervus),
    Parent2_Cervus = as.character(Parent2_Cervus),
        Plot_Mother = dplyr::if_else(Parent1_Cervus == Mother_ID, Plot_Parent1,
                                 dplyr::if_else(Parent2_Cervus == Mother_ID, Plot_Parent2, NA_character_)),
    Plot_Father = dplyr::if_else(Parent1_Cervus == Father_ID, Plot_Parent1,
                                 dplyr::if_else(Parent2_Cervus == Father_ID, Plot_Parent2, NA_character_)),
        DBH_Mother = dplyr::if_else(Parent1_Cervus == Mother_ID, DBH_Parent1,
                                dplyr::if_else(Parent2_Cervus == Mother_ID, DBH_Parent2, NA_real_)),
    DBH_Father = dplyr::if_else(Parent1_Cervus == Father_ID, DBH_Parent1,
                                dplyr::if_else(Parent2_Cervus == Father_ID, DBH_Parent2, NA_real_))
  )

# Prepare cluster table
clusters_min <- clusters %>%
  select(OffspringID, ClusterIndex) %>%
  mutate(
    OffspringID  = as.character(str_trim(OffspringID)),
    ClusterIndex = as.character(ClusterIndex)
  ) %>%
  distinct(OffspringID, .keep_all = TRUE)

# Build final df_map 
df_map <- filtered_aug %>%
  mutate(OffspringID = as.character(str_trim(OffspringID))) %>%
  left_join(clusters_min, by = "OffspringID") %>%
  rename(Offspring_Cluster = ClusterIndex) %>%
  mutate(
    Mother_Cluster = Offspring_Cluster,
    Father_Cluster = Offspring_Cluster
  ) %>%
  select(
    OffspringID, Plot_Offspring,
    Latitude_Offspring, Longitude_Offspring, DBH_Offspring,
    Offspring_Cluster,
    Mother_ID, Plot_Mother, Latitude_Mother, Longitude_Mother, DBH_Mother, Distance_To_Mother, Mother_Cluster,
    Father_ID, Plot_Father, Latitude_Father, Longitude_Father, DBH_Father, Distance_To_Father, Father_Cluster
  )

# Export 
write_csv(df_map, output_path)
cat("Final dataframe exported to:", output_path, "\n")



###########################################################################################



########################################################################
######################### 2 plots map - visualization ##################
########################################################################

########################################################################
######################### 2 plots map - visualization ##################
########################################################################

library(readr)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(ggspatial)

# ------------------ Paths ------------------
path_dfmap <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis/df_map.csv"
path_isolines <- "C:/Users/bonni/Desktop/Fichiers_cartes_Qgis/Isolignes/Isolignes_Regina_5m/fr_662043116_lidar_regstgmult2013_02.shp"
output_plot <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis/parentage_map.png"

# ------------------ Read data ------------------
df_map <- read_csv(path_dfmap, show_col_types = FALSE)
isolines <- st_read(path_isolines, quiet = TRUE)

# ------------------ Create sf objects (initially WGS84) ------------------
offspring_sf <- st_as_sf(df_map,
                         coords = c("Longitude_Offspring","Latitude_Offspring"),
                         crs = 4326)

mother_sf <- st_as_sf(df_map,
                      coords = c("Longitude_Mother","Latitude_Mother"),
                      crs = 4326)

father_sf <- st_as_sf(df_map,
                      coords = c("Longitude_Father","Latitude_Father"),
                      crs = 4326)

# ------------------ Create line features ------------------
# Seed dispersal (offspring → mother, red)
seed_lines <- lapply(1:nrow(df_map), function(i) {
  st_linestring(matrix(c(df_map$Longitude_Offspring[i], df_map$Latitude_Offspring[i],
                         df_map$Longitude_Mother[i],    df_map$Latitude_Mother[i]),
                       ncol=2, byrow=TRUE))
})
seed_sf <- st_sfc(seed_lines, crs = 4326) %>%
  st_sf(OffspringID = df_map$OffspringID)

# Pollen dispersal (offspring → father, blue)
pollen_lines <- lapply(1:nrow(df_map), function(i) {
  st_linestring(matrix(c(df_map$Longitude_Offspring[i], df_map$Latitude_Offspring[i],
                         df_map$Longitude_Father[i],    df_map$Latitude_Father[i]),
                       ncol=2, byrow=TRUE))
})
pollen_sf <- st_sfc(pollen_lines, crs = 4326) %>%
  st_sf(OffspringID = df_map$OffspringID)

# ------------------ Reproject all layers to isolines CRS ------------------
target_crs <- st_crs(isolines)

offspring_sf <- st_transform(offspring_sf, target_crs)
mother_sf    <- st_transform(mother_sf, target_crs)
father_sf    <- st_transform(father_sf, target_crs)
seed_sf      <- st_transform(seed_sf, target_crs)
pollen_sf    <- st_transform(pollen_sf, target_crs)

# ------------------ Crop isolines to sampling extent ------------------
# Collect all coordinates into one multipoint
sampling_coords <- rbind(
  st_coordinates(offspring_sf),
  st_coordinates(mother_sf),
  st_coordinates(father_sf)
)
sampling_points <- st_multipoint(sampling_coords)
sampling_points_sfc <- st_sfc(sampling_points, crs = target_crs)

# Buffer and bbox around sampling points (e.g., 200 m buffer)
bbox <- st_as_sfc(st_bbox(st_buffer(sampling_points_sfc, dist = 200)))

# Crop isolines
isolines_crop <- st_intersection(isolines, bbox)

# Convert cluster to factor for proper discrete colors
offspring_sf$Offspring_Cluster <- as.factor(df_map$Offspring_Cluster)
mother_sf$Offspring_Cluster    <- as.factor(df_map$Offspring_Cluster)
father_sf$Offspring_Cluster    <- as.factor(df_map$Offspring_Cluster)
seed_sf$Offspring_Cluster      <- as.factor(df_map$Offspring_Cluster)
pollen_sf$Offspring_Cluster    <- as.factor(df_map$Offspring_Cluster)

# Define custom colors for clusters
my_colors <- c(
  "1" = "#CD5555",  
  "3" = "#00BFFF",  
  "7" = "#4DAF4A",  
  "8" = "#984EA3",  
  "9" = "#FF7F00",  
  "10" = "#FFFF33", 
  "11" = "#A65628",  
  "12" = "#F781BF",  
  "13" = "#999999",  
  "14" = "#66C2A5"  
)


# ------------------ Plot ------------------

p <- ggplot() +
  # Topographic isolines
  geom_sf(data = isolines_crop, color = "grey70", size = 0.3) +
  # Seed dispersal lines (red)
  geom_sf(data = seed_sf, aes(linetype = "Seed dispersal (to mother)"), 
          color = "red", size = 0.4) +
  # Parents (squares, grouped)
  geom_sf(data = mother_sf, aes(color = Offspring_Cluster, shape = "Parent"), size = 2) +
  geom_sf(data = father_sf, aes(color = Offspring_Cluster, shape = "Parent"), size = 2) +
  # Offspring (circles)
  geom_sf(data = offspring_sf, aes(color = Offspring_Cluster, shape = "Offspring"), size = 2) +
  # Scales
  scale_shape_manual(name = "Individuals",
                     values = c("Offspring" = 16, "Parent" = 15)) +
  scale_linetype_manual(name = "Dispersal",
                        values = c("Seed dispersal (to mother)" = "solid")) +
  scale_color_manual(name = "Family cluster", values = my_colors) +
  theme_minimal() +
  labs(title = "Parent–offspring seed dispersal with topographic isolines",
       x = "Longitude", y = "Latitude") +
  # Add scale bar in bottom-left
  annotation_scale(location = "bl", 
                   width_hint = 0.3,      # proportion of map width
                   line_width = 0.5, 
                   text_cex = 0.8)

print(p)






########################################################################
################## Dispersal counts within vs between plots ############
########################################################################

library(readr)
library(dplyr)
library(stringr)


# ---- 1) Basic sanity: keep rows with available plot info for parents ----
df_seed  <- df_map %>% filter(!is.na(Plot_Mother), !is.na(Plot_Offspring))
df_pollen <- df_map %>% filter(!is.na(Plot_Father), !is.na(Plot_Offspring))

# ---- 2) SEED dispersal: mother (source) -> offspring (destination) ----
seed_counts <- df_seed %>%
  mutate(
    flow_type = if_else(Plot_Mother == Plot_Offspring, "within_plot", "between_plots"),
    direction = paste0(Plot_Mother, " → ", Plot_Offspring)
  ) %>%
  count(flow_type, direction, name = "n") %>%
  arrange(flow_type, desc(n))

# Overall totals for seed
seed_totals <- df_seed %>%
  mutate(between = Plot_Mother != Plot_Offspring) %>%
  summarise(
    n_total = n(),
    n_within = sum(!between),
    n_between = sum(between),
    pct_between = round(100 * n_between / n_total, 1)
  )

# Optional: mean distances by flow_type (useful QA)
seed_distance_summary <- df_seed %>%
  mutate(flow_type = if_else(Plot_Mother == Plot_Offspring, "within_plot", "between_plots")) %>%
  group_by(flow_type) %>%
  summarise(
    n = n(),
    mean_distance_m = mean(Distance_To_Mother, na.rm = TRUE),
    median_distance_m = median(Distance_To_Mother, na.rm = TRUE),
    max_distance_m = max(Distance_To_Mother, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 3) POLLEN dispersal: father (source) -> offspring (destination) ----
pollen_counts <- df_pollen %>%
  mutate(
    flow_type = if_else(Plot_Father == Plot_Offspring, "within_plot", "between_plots"),
    direction = paste0(Plot_Father, " → ", Plot_Offspring)
  ) %>%
  count(flow_type, direction, name = "n") %>%
  arrange(flow_type, desc(n))

# Overall totals for pollen
pollen_totals <- df_pollen %>%
  mutate(between = Plot_Father != Plot_Offspring) %>%
  summarise(
    n_total = n(),
    n_within = sum(!between),
    n_between = sum(between),
    pct_between = round(100 * n_between / n_total, 1)
  )

# Optional: mean distances by flow_type (useful QA)
pollen_distance_summary <- df_pollen %>%
  mutate(flow_type = if_else(Plot_Father == Plot_Offspring, "within_plot", "between_plots")) %>%
  group_by(flow_type) %>%
  summarise(
    n = n(),
    mean_distance_m = mean(Distance_To_Father, na.rm = TRUE),
    median_distance_m = median(Distance_To_Father, na.rm = TRUE),
    max_distance_m = max(Distance_To_Father, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 4) Nice 2x2 directional matrices (HKO50 vs PAI74) ----
# Seed matrix: rows = source (Mother plot), cols = destination (Offspring plot)
seed_matrix <- df_seed %>%
  count(Plot_Mother, Plot_Offspring, name = "n") %>%
  tidyr::pivot_wider(names_from = Plot_Offspring, values_from = n, values_fill = 0) %>%
  arrange(Plot_Mother)

# Pollen matrix: rows = source (Father plot), cols = destination (Offspring plot)
pollen_matrix <- df_pollen %>%
  count(Plot_Father, Plot_Offspring, name = "n") %>%
  tidyr::pivot_wider(names_from = Plot_Offspring, values_from = n, values_fill = 0) %>%
  arrange(Plot_Father)

# ---- 5) Print quick summaries to console ----
cat("\n=== SEED dispersal (mother → offspring) ===\n")
print(seed_totals)
print(seed_counts)
print(seed_distance_summary)

cat("\n=== POLLEN dispersal (father → offspring) ===\n")
print(pollen_totals)
print(pollen_counts)
print(pollen_distance_summary)

cat("\n=== Directional matrices ===\n")
cat("\n[SEED] rows = Mother plot, cols = Offspring plot\n")
print(seed_matrix)
cat("\n[POLLEN] rows = Father plot, cols = Offspring plot\n")
print(pollen_matrix)

# ---- 6) (Optional) Export summaries to CSV next to your map outputs ----
out_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis"
readr::write_csv(seed_counts, file.path(out_dir, "seed_counts_direction.csv"))
readr::write_csv(seed_totals, file.path(out_dir, "seed_totals.csv"))
readr::write_csv(seed_distance_summary, file.path(out_dir, "seed_distance_summary.csv"))

readr::write_csv(pollen_counts, file.path(out_dir, "pollen_counts_direction.csv"))
readr::write_csv(pollen_totals, file.path(out_dir, "pollen_totals.csv"))
readr::write_csv(pollen_distance_summary, file.path(out_dir, "pollen_distance_summary.csv"))

readr::write_csv(seed_matrix, file.path(out_dir, "seed_direction_matrix.csv"))
readr::write_csv(pollen_matrix, file.path(out_dir, "pollen_direction_matrix.csv"))



########################################################################
############ Parental reproductive success & dispersal distances #######
########################################################################

library(readr)
library(dplyr)
library(ggplot2)

# ------------------ Paths ------------------
path_dfmap <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis/df_map.csv"
output_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis"

# ------------------ Read ------------------
df_map <- read_csv(path_dfmap, show_col_types = FALSE)

# ------------------ SEED (mothers) ------------------
mother_summary <- df_map %>%
  filter(!is.na(Mother_ID), !is.na(Plot_Mother)) %>%
  group_by(Plot_Mother, Mother_ID) %>%
  summarise(
    n_offspring = n(),                                   # number of offspring
    mean_distance = mean(Distance_To_Mother, na.rm = TRUE), # mean dispersal distance
    median_distance = median(Distance_To_Mother, na.rm = TRUE),
    max_distance = max(Distance_To_Mother, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Plot_Mother, desc(n_offspring))

# ------------------ POLLEN (fathers) ------------------
father_summary <- df_map %>%
  filter(!is.na(Father_ID), !is.na(Plot_Father)) %>%
  group_by(Plot_Father, Father_ID) %>%
  summarise(
    n_offspring = n(),
    mean_distance = mean(Distance_To_Father, na.rm = TRUE),
    median_distance = median(Distance_To_Father, na.rm = TRUE),
    max_distance = max(Distance_To_Father, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Plot_Father, desc(n_offspring))

# ------------------ EXPORT ------------------
write_csv(mother_summary, file.path(output_dir, "mother_reproductive_success.csv"))
write_csv(father_summary, file.path(output_dir, "father_reproductive_success.csv"))

# ------------------ QUICK OVERVIEW ------------------
cat("\n=== Mothers ===\n")
print(head(mother_summary, 10))
cat("\n=== Fathers ===\n")
print(head(father_summary, 10))


# Define custom colors for plots
plot_colors <- c("HKO50" = "#CDAD00", "PAI74" = "mediumorchid4")

# ------------------ Visualisation: number of offspring ------------------
p1 <- ggplot(mother_summary, aes(x = Plot_Mother, y = n_offspring, fill = Plot_Mother)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = plot_colors) +
  theme_minimal() +
  labs(title = "Number of offspring per mother", x = "Plot", y = "Offspring count")

p2 <- ggplot(father_summary, aes(x = Plot_Father, y = n_offspring, fill = Plot_Father)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = plot_colors) +
  theme_minimal() +
  labs(title = "Number of offspring per father", x = "Plot", y = "Offspring count")

# ------------------ Visualisation: mean dispersal distances ------------------
p3 <- ggplot(mother_summary, aes(x = Plot_Mother, y = mean_distance, fill = Plot_Mother)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = plot_colors) +
  theme_minimal() +
  labs(title = "Mean seed dispersal distance per mother", x = "Plot", y = "Distance (m)")

p4 <- ggplot(father_summary, aes(x = Plot_Father, y = mean_distance, fill = Plot_Father)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = plot_colors) +
  theme_minimal() +
  labs(title = "Mean pollen dispersal distance per father", x = "Plot", y = "Distance (m)")

# ------------------ Visualisation: density plots (dispersal distances) ------------------
p5 <- ggplot(df_map, aes(x = Distance_To_Mother, fill = Plot_Mother, color = Plot_Mother)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(limits = c(0, 2000)) +
  theme_minimal() +
  labs(title = "Distribution of seed dispersal distances", x = "Distance (m)", y = "Density")

p6 <- ggplot(df_map, aes(x = Distance_To_Father, fill = Plot_Father, color = Plot_Father)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(limits = c(0, 2000)) +
  theme_minimal() +
  labs(title = "Distribution of pollen dispersal distances", x = "Distance (m)", y = "Density")



########################################################################
############# Statistical comparison between plots #####################
########################################################################

# ---- Mothers: number of offspring ----
wilcox_mother_offspring <- wilcox.test(
  n_offspring ~ Plot_Mother,
  data = mother_summary
)
cat("\n[Wilcoxon] Number of offspring per mother (HKO50 vs PAI74):\n")
print(wilcox_mother_offspring)

# ---- Fathers: number of offspring ----
wilcox_father_offspring <- wilcox.test(
  n_offspring ~ Plot_Father,
  data = father_summary
)
cat("\n[Wilcoxon] Number of offspring per father (HKO50 vs PAI74):\n")
print(wilcox_father_offspring)

# ---- Mothers: mean dispersal distance ----
wilcox_mother_distance <- wilcox.test(
  mean_distance ~ Plot_Mother,
  data = mother_summary
)
cat("\n[Wilcoxon] Mean seed dispersal distance per mother (HKO50 vs PAI74):\n")
print(wilcox_mother_distance)

# ---- Fathers: mean dispersal distance ----
wilcox_father_distance <- wilcox.test(
  mean_distance ~ Plot_Father,
  data = father_summary
)
cat("\n[Wilcoxon] Mean pollen dispersal distance per father (HKO50 vs PAI74):\n")
print(wilcox_father_distance)


########################################################################
########## Combined plot: seed vs pollen dispersal by plot #############
########################################################################

library(dplyr)
library(ggplot2)
library(tidyr)

# ---- Reshape data into long format ----
df_long <- df_map %>%
  select(Plot_Offspring, Plot_Mother, Distance_To_Mother,
         Plot_Father, Distance_To_Father) %>%
  # Gather seed and pollen distances into one column
  pivot_longer(
    cols = c(Distance_To_Mother, Distance_To_Father),
    names_to = "type",
    values_to = "distance"
  ) %>%
  mutate(
    type = ifelse(type == "Distance_To_Mother", "Seed", "Pollen"),
    plot = ifelse(type == "Seed", Plot_Mother, Plot_Father)
  ) %>%
  filter(!is.na(distance), !is.na(plot))


# Combined boxplot 
p_combined_box <- ggplot(df_long, aes(x = type, y = distance, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ plot, ncol = 2) +
  scale_fill_manual(values = c("Seed" = "tan1", "Pollen" = "#F0E68C")) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_minimal() +
  labs(title = "Seed vs Pollen dispersal distances per plot",
       x = "Dispersal type", y = "Distance (m)")

print(p_combined_box)


########################################################################
########## Barplots of reproductive success per parent #################
########################################################################

library(readr)
library(dplyr)
library(ggplot2)

# ------------------ Paths ------------------
path_dfmap <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis/df_map.csv"
output_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.4-2_plots_analysis"

# ------------------ Read ------------------
df_map <- read_csv(path_dfmap, show_col_types = FALSE)

# ------------------ Custom plot colors ------------------
plot_colors <- c("HKO50" = "#CDAD00", "PAI74" = "mediumorchid4")

# ------------------ Summaries ------------------
mother_summary <- df_map %>%
  filter(!is.na(Mother_ID), !is.na(Plot_Mother)) %>%
  group_by(Plot_Mother, Mother_ID) %>%
  summarise(n_offspring = n(), .groups = "drop")

father_summary <- df_map %>%
  filter(!is.na(Father_ID), !is.na(Plot_Father)) %>%
  group_by(Plot_Father, Father_ID) %>%
  summarise(n_offspring = n(), .groups = "drop")

# ------------------ Function for plotting ------------------
plot_bar <- function(data, id_col, plot_col, title) {
  data <- data %>%
    arrange(desc(n_offspring)) %>%
    slice_head(n = 25)
  
  id_var <- sym(id_col)
  plot_var <- sym(plot_col)
  
  id_levels <- data %>%
    arrange(desc(n_offspring)) %>%
    pull(!!id_var) %>%
    unique()
  
  data <- data %>%
    mutate(!!id_var := factor(!!id_var, levels = id_levels))
  
  ggplot(data, aes(x = !!id_var, y = n_offspring, fill = !!plot_var)) +
    geom_col() +
    scale_fill_manual(values = plot_colors, name = "Plot") +
    theme_minimal(base_size = 14) +  # increase global base font size
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title   = element_text(size = 15, face = "bold")
    ) +
    labs(title = title, x = "Parent ID", y = "Number of offspring")
}


# Mothers HKO50
p_mother_hko50 <- plot_bar(
  mother_summary %>% filter(Plot_Mother == "HKO50"),
  id_col = "Mother_ID", plot_col = "Plot_Mother",
  title = "Top 25 mothers - HKO50"
)

# Mothers PAI74
p_mother_pai74 <- plot_bar(
  mother_summary %>% filter(Plot_Mother == "PAI74"),
  id_col = "Mother_ID", plot_col = "Plot_Mother",
  title = "Top 25 mothers - PAI74"
)

# Fathers HKO50
p_father_hko50 <- plot_bar(
  father_summary %>% filter(Plot_Father == "HKO50"),
  id_col = "Father_ID", plot_col = "Plot_Father",
  title = "Top 25 fathers - HKO50"
)

# Fathers PAI74
p_father_pai74 <- plot_bar(
  father_summary %>% filter(Plot_Father == "PAI74"),
  id_col = "Father_ID", plot_col = "Plot_Father",
  title = "Top 25 fathers - PAI74"
)

# ------------------ Save plots ------------------
ggsave(file.path(output_dir, "barplot_mothers_HKO50.png"), p_mother_hko50, width = 11, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "barplot_mothers_PAI74.png"), p_mother_pai74, width = 11, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "barplot_fathers_HKO50.png"), p_father_hko50, width = 11, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "barplot_fathers_PAI74.png"), p_father_pai74, width = 11, height = 7, dpi = 300, bg = "white")






