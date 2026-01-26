# ============================================================
# Script: 04.1-visualisation_cervus_results.R
# Purpose: Visualize high-confidence parent-offspring links from CERVUS, map them (optional isolines), and summarize dispersal distances.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: summary_<PLOT>.csv ; coord_<PLOT>.csv ; optional isolines shapefile
# Outputs: figures/04_parentage/parent_offspring_map_<PLOT>.png ; figures/04_parentage/histogram_distances_<PLOT>.png
# Figures/Tables: Parent-offspring map + distance histogram
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
plot_name <- "PAI74"
plot_color <- "mediumorchid4"

base_path <- "."

# Input files (relative paths)
results_file <- file.path(base_path, "results", "cervus", plot_name, paste0("summary_", plot_name, ".csv"))
coords_file <- file.path(base_path, "data", plot_name, paste0("coord_", plot_name, ".csv"))

# Optional isolines (set to NA to disable). Must be in WGS84 or will be transformed.
isolines_file <- NA_character_
# Example (if you keep isolines in repo):
# isolines_file <- file.path(base_path, "data", "isolines", "Isolignes_Regina_5m.shp")

# Outputs
fig_dir <- file.path("figures", "04_parentage", plot_name)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_map <- file.path(fig_dir, paste0("parent_offspring_map_", plot_name, ".png"))
out_hist <- file.path(fig_dir, paste0("histogram_distances_", plot_name, ".png"))

# Crop buffer around study area (degrees in WGS84)
bbox_buffer <- 0.001
## ======================================================

## ===================== LIBRARIES ======================
library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(geosphere)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(results_file))
stopifnot(file.exists(coords_file))

results <- read.csv(results_file, sep = ",", stringsAsFactors = FALSE)
geo_data <- read.csv(coords_file, sep = ",", stringsAsFactors = FALSE)

# Basic checks
stopifnot(all(c("ID", "lat", "long") %in% names(geo_data)))
stopifnot(all(c("Offspring_ID", "First_candidate_ID", "Second_candidate_ID",
                "Pair_confidence1", "Pair_confidence2") %in% names(results)))
## ======================================================

## ===================== MERGE COORDS ===================
results_geo <- results %>%
  left_join(geo_data, by = c("Offspring_ID" = "ID")) %>%
  rename(lat_offspring = lat, long_offspring = long) %>%
  left_join(geo_data, by = c("First_candidate_ID" = "ID")) %>%
  rename(lat_first_parent = lat, long_first_parent = long) %>%
  left_join(geo_data, by = c("Second_candidate_ID" = "ID")) %>%
  rename(lat_second_parent = lat, long_second_parent = long)
## ======================================================

## ===================== LINKS (HIGH CONF) ==============
links_first <- results_geo %>%
  filter(Pair_confidence1 %in% c("*", "+")) %>%
  select(
    Offspring_ID, lat_offspring, long_offspring,
    First_candidate_ID, lat_first_parent, long_first_parent, Pair_confidence1
  ) %>%
  rename(
    parent_id = First_candidate_ID,
    lat_parent = lat_first_parent,
    long_parent = long_first_parent,
    confidence = Pair_confidence1
  )

links_second <- results_geo %>%
  filter(Pair_confidence2 %in% c("*", "+")) %>%
  select(
    Offspring_ID, lat_offspring, long_offspring,
    Second_candidate_ID, lat_second_parent, long_second_parent, Pair_confidence2
  ) %>%
  rename(
    parent_id = Second_candidate_ID,
    lat_parent = lat_second_parent,
    long_parent = long_second_parent,
    confidence = Pair_confidence2
  )

links <- bind_rows(links_first, links_second) %>%
  filter(
    !is.na(lat_offspring), !is.na(long_offspring),
    !is.na(lat_parent), !is.na(long_parent)
  )

cat("Total number of high-confidence parent-offspring links:", nrow(links), "\n")
## ======================================================

## ===================== PLOT (BASIC) ===================
basic_plot <- ggplot() +
  geom_segment(
    data = links,
    aes(x = long_offspring, y = lat_offspring, xend = long_parent, yend = lat_parent),
    color = "gray40", linewidth = 0.5
  ) +
  geom_point(
    data = links,
    aes(x = long_offspring, y = lat_offspring, shape = "Offspring", color = "Offspring"),
    size = 3
  ) +
  geom_point(
    data = links,
    aes(x = long_parent, y = lat_parent, shape = "Parent", color = "Parent"),
    size = 3
  ) +
  scale_color_manual(values = c("Offspring" = "forestgreen", "Parent" = "royalblue")) +
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17)) +
  theme_minimal() +
  labs(
    title = paste("Parent-Offspring Links –", plot_name),
    x = "Longitude", y = "Latitude", color = "Type", shape = "Type"
  )

print(basic_plot)
## ======================================================

## ============ PLOT WITH OPTIONAL TOPOGRAPHIC ISOLINES ============
# Build bbox from points
coords_all <- data.frame(
  long = c(links$long_offspring, links$long_parent),
  lat  = c(links$lat_offspring, links$lat_parent)
)
coords_sf <- st_as_sf(coords_all, coords = c("long", "lat"), crs = 4326)

bbox <- st_bbox(coords_sf)
bbox_exp <- bbox
bbox_exp["xmin"] <- bbox["xmin"] - bbox_buffer
bbox_exp["ymin"] <- bbox["ymin"] - bbox_buffer
bbox_exp["xmax"] <- bbox["xmax"] + bbox_buffer
bbox_exp["ymax"] <- bbox["ymax"] + bbox_buffer

isolines_crop <- NULL
if (!is.na(isolines_file) && file.exists(isolines_file)) {
  isolines <- st_read(isolines_file, quiet = TRUE)
  isolines_wgs84 <- st_transform(isolines, crs = 4326)
  isolines_crop <- st_crop(isolines_wgs84, bbox_exp)
}

parentage_plot <- ggplot() +
  { if (!is.null(isolines_crop)) geom_sf(data = isolines_crop, color = "grey50", linewidth = 0.3, alpha = 0.7) } +
  geom_segment(
    data = links,
    aes(x = long_offspring, y = lat_offspring, xend = long_parent, yend = lat_parent),
    color = "darkorange", linewidth = 0.7, alpha = 0.8
  ) +
  geom_point(
    data = links,
    aes(x = long_offspring, y = lat_offspring, shape = "Offspring"),
    color = "forestgreen", size = 3
  ) +
  geom_point(
    data = links,
    aes(x = long_parent, y = lat_parent, shape = "Parent"),
    color = "royalblue", size = 3
  ) +
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
  theme_minimal() +
  theme(
    legend.position = c(0.99, 0.9),
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

print(parentage_plot)

ggsave(
  filename = out_map,
  plot = parentage_plot,
  width = 10, height = 8, dpi = 300, bg = "white"
)

cat("Saved map:", out_map, "\n")
## ======================================================

## ===================== DISTANCES ======================
links <- links %>%
  mutate(
    distance_m = distHaversine(
      cbind(long_offspring, lat_offspring),
      cbind(long_parent, lat_parent)
    )
  )

distance_stats <- links %>%
  summarise(
    mean_distance_m = mean(distance_m, na.rm = TRUE),
    min_distance_m  = min(distance_m, na.rm = TRUE),
    max_distance_m  = max(distance_m, na.rm = TRUE)
  )

print(distance_stats)
## ======================================================

## ===================== HISTOGRAM ======================
hist_distance <- ggplot(links, aes(x = distance_m)) +
  geom_histogram(binwidth = 10, fill = plot_color, color = "white", alpha = 0.8) +
  geom_vline(
    aes(xintercept = mean(distance_m, na.rm = TRUE)),
    linetype = "dashed", color = "black", linewidth = 0.8
  ) +
  labs(
    title = paste("Histogram of Parent-Offspring Distances –", plot_name),
    x = "Distance (m)",
    y = "Number of pairs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(hist_distance)

ggsave(
  filename = out_hist,
  plot = hist_distance,
  width = 8, height = 6, dpi = 300, bg = "white"
)

cat("Saved histogram:", out_hist, "\n")
## ======================================================
