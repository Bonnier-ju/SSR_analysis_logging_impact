############### Visualisation of Cervus results ###############
###############################################################

library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(geosphere)
library(readr)

# Load CERVUS result file
results <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/PAI74/summary_PAI74.csv", sep = ",")
geo_data <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/PAI74/coord_PAI74.csv", sep = ",")

# === Define plot name and color ===
plot_name <- "PAI74"
plot_color <- "mediumorchid4"

# Merge coordinates with offspring and candidate parents
results_geo <- results %>%
  left_join(geo_data, by = c("Offspring_ID" = "ID")) %>%
  rename(lat_offspring = lat, long_offspring = long) %>%
  left_join(geo_data, by = c("First_candidate_ID" = "ID")) %>%
  rename(lat_first_parent = lat, long_first_parent = long) %>%
  left_join(geo_data, by = c("Second_candidate_ID" = "ID")) %>%
  rename(lat_second_parent = lat, long_second_parent = long)

# Extract high-confidence pair links
links_first <- results_geo %>%
  filter(Pair_confidence1 %in% c("*", "+")) %>%
  select(Offspring_ID, lat_offspring, long_offspring,
         First_candidate_ID, lat_first_parent, long_first_parent, Pair_confidence1) %>%
  rename(parent_id = First_candidate_ID,
         lat_parent = lat_first_parent,
         long_parent = long_first_parent,
         confidence = Pair_confidence1)

links_second <- results_geo %>%
  filter(Pair_confidence2 %in% c("*", "+")) %>%
  select(Offspring_ID, lat_offspring, long_offspring,
         Second_candidate_ID, lat_second_parent, long_second_parent, Pair_confidence2) %>%
  rename(parent_id = Second_candidate_ID,
         lat_parent = lat_second_parent,
         long_parent = long_second_parent,
         confidence = Pair_confidence2)

links <- bind_rows(links_first, links_second)
head(links)

# Plot results
ggplot() +
  # Draw parent-offspring links (single color)
  geom_segment(data = links, aes(x = long_offspring, y = lat_offspring,
                                 xend = long_parent, yend = lat_parent),
               color = "gray40", size = 0.5) +
  
  # Draw offspring points
  geom_point(data = links, aes(x = long_offspring, y = lat_offspring, shape = "Offspring", color = "Offspring"), size = 3) +
  
  # Draw parent points
  geom_point(data = links, aes(x = long_parent, y = lat_parent, shape = "Parent", color = "Parent"), size = 3) +
  
  # Manual color and shape scales
  scale_color_manual(values = c("Offspring" = "forestgreen", "Parent" = "royalblue")) +
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17)) +
  
  theme_minimal() +
  labs(title = paste("Parent-Offspring Links –", plot_name),
       x = "Longitude", y = "Latitude", color = "Type", shape = "Type")




#### Plot results with topographic line ####
# Load isolines and transform CRS
file_path_isolines <- "C:/Users/bonni/Desktop/Fichiers_cartes_Qgis/Isolignes/Isolignes_Regina_5m/fr_662043116_lidar_regstgmult2013_02.shp"
isolines <- st_read(file_path_isolines, quiet = TRUE)
isolines_wgs84 <- st_transform(isolines, crs = 4326)

# Create sf points from parent-offspring coordinates
coords_all <- data.frame(
  ID = c(links$Offspring.ID, links$parent_id),
  long = c(links$long_offspring, links$long_parent),
  lat = c(links$lat_offspring, links$lat_parent)
)
coords_sf <- st_as_sf(coords_all, coords = c("long", "lat"), crs = 4326)

# Compute bounding box and expand
bbox <- st_bbox(coords_sf)
bbox_exp <- structure(
  c(
    xmin = bbox["xmin"] - 0.001,
    ymin = bbox["ymin"] - 0.001,
    xmax = bbox["xmax"] + 0.001,
    ymax = bbox["ymax"] + 0.001
  ),
  class = "bbox",
  crs = st_crs(isolines_wgs84)
)

# Crop isolines
isolines_crop <- st_crop(isolines_wgs84, bbox_exp)

# Create plot object
parentage_plot <- ggplot() +
  geom_sf(data = isolines_crop, color = "grey50", size = 0.3, alpha = 0.7) +
  geom_segment(data = links, aes(x = long_offspring, y = lat_offspring,
                                 xend = long_parent, yend = lat_parent),
               color = "darkorange", size = 0.7, alpha = 0.8) +
  geom_point(data = links, aes(x = long_offspring, y = lat_offspring, shape = "Offspring"),
             color = "forestgreen", size = 3) +
  geom_point(data = links, aes(x = long_parent, y = lat_parent, shape = "Parent"),
             color = "royalblue", size = 3) +
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
  theme_minimal() +
  theme(
    legend.position = c(0.99, 0.9),  # Top right inside the plot
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white")
  ) +
  labs(
    title = paste("Parent-Offspring Links –", plot_name),
    x = "Longitude", y = "Latitude"
  ) +
  coord_sf(
    xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
    ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
    expand = FALSE
  )

# Display the plot
print(parentage_plot)

# Save as PNG
ggsave(
  filename = paste0("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/Results_graphs/parent_offspring_map_", plot_name, ".png"),
  plot = parentage_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

#### Analysis on distances parents_offspring

# Compute Haversine distance in meters between each offspring and its assigned parent
links <- links %>%
  mutate(
    distance = distHaversine(
      cbind(long_offspring, lat_offspring),
      cbind(long_parent, lat_parent)
    )
  )

# Summary statistics of parent-offspring distances
distance_stats <- links %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    min_distance = min(distance, na.rm = TRUE),
    max_distance = max(distance, na.rm = TRUE)
  )

# Print the results
print(distance_stats)

# Count total number of links
total_links <- nrow(links)
cat("Total number of parent-offspring links:", total_links, "\n")


# histogram of parents-offspring distances
hist_distance <- ggplot(links, aes(x = distance)) +
  geom_histogram(binwidth = 10, fill = plot_color, color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(distance, na.rm = TRUE)),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(
    title = paste("Histogram of Parent-Offspring Distances –", plot_name),
    x = "Distance (meters)",
    y = "Number of Pairs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
print(hist_distance)

ggsave(
  filename = paste0("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/Results_graphs/histogram_distances_", plot_name, ".png"),
  plot = hist_distance,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
