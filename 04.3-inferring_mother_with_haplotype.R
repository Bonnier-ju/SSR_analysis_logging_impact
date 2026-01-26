# ============================================================
# Script: 04.3_compare_CERVUS_COLONY_infer_mother_haplotype.R
# Purpose:
#   (1) Compare parentage assignments between CERVUS and COLONY
#   (2) Use chloroplast haplotypes to infer the mother among two CERVUS candidates
#   (3) Compute seed (mother-offspring) and pollen (father-offspring) distances
#   (4) Map links (optionally with topographic isolines)
#   (5) Plot seed/pollen distance distributions (single plot + all plots combined)
#
# Notes:
#   - Maternal inference uses chloroplast haplotype match (offspring vs candidate parents).
#   - Keeps only cases where EXACTLY ONE of the two candidates matches offspring haplotype.
#   - Make sure haplotype file has: ID, plot, lat, long, dbh, alt, Haplotypes
#   - Make sure CERVUS file has: Offspring_ID, First_candidate_ID, Second_candidate_ID,
#       Pair_confidence1, Pair_confidence2
#   - Make sure COLONY BestConfig has: OffspringID, FatherID, MotherID
#
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
plot_name <- "PAI74"

# --- Input files (edit paths) ---
file_cervus <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/PAI74/summary_PAI74.csv"
file_colony <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.2-parentage_Colony/PAI74/results_colony_PAI74/PAI74_colony_projet.BestConfig.csv"
file_haplo  <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/Chloroplatics_markers/Haplotype_full.csv"

# Optional: isolines (set to NA if you don't want isolines)
file_isolines <- "C:/Users/bonni/Desktop/Fichiers_cartes_Qgis/Isolignes/Isolignes_Regina_5m/fr_662043116_lidar_regstgmult2013_02.shp"
bbox_expand_deg <- 0.001  # expansion in degrees for map zoom

# --- Outputs ---
out_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_comparison_csv <- file.path(out_dir, paste0("comparison_data_", plot_name, ".csv"))
out_filtered_csv   <- file.path(out_dir, paste0("filtered_data_", plot_name, ".csv"))

out_map_basic_png   <- file.path(out_dir, paste0("map_mother_father_", plot_name, ".png"))
out_map_isoline_png <- file.path(out_dir, paste0("map_mother_father_isolines_", plot_name, ".png"))
out_map_mother_png  <- file.path(out_dir, paste0("map_mother_only_isolines_", plot_name, ".png"))

out_seed_hist_png   <- file.path(out_dir, paste0("seed_distance_hist_", plot_name, ".png"))
out_pollen_hist_png <- file.path(out_dir, paste0("pollen_distance_hist_", plot_name, ".png"))

# Multi-plot combined distributions (optional)
make_all_plots_combined <- TRUE
file_filtered_SPR <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_SPR.csv"
file_filtered_REG <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_REG.csv"
file_filtered_NOU <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_NOU.csv"
file_filtered_PAR <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_PAR.csv"

out_seed_global_png  <- file.path(out_dir, "seed_distance_global_all_plots.png")
out_pollen_global_png <- file.path(out_dir, "pollen_distance_global_all_plots.png")

# Plot style
plot_color <- "mediumorchid4"
## ======================================================

## ===================== LIBRARIES ======================
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(geosphere)

# isolines map
suppressWarnings({
  ok_sf <- requireNamespace("sf", quietly = TRUE)
})
## ======================================================

## ===================== HELPERS =========================
stop_if_missing <- function(x) {
  if (!file.exists(x)) stop("File not found: ", x)
}

# Robust distance (meters)
dist_m <- function(lon1, lat1, lon2, lat2) {
  ifelse(is.na(lon1) | is.na(lat1) | is.na(lon2) | is.na(lat2),
         NA_real_,
         geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)))
}
## ======================================================

################################################################################
######################## 1) LOAD + FILTER CERVUS ###############################
################################################################################

stop_if_missing(file_cervus)
stop_if_missing(file_colony)
stop_if_missing(file_haplo)

# If your CERVUS CSV uses ";" delimiter, switch to: read_delim(file_cervus, delim=";")
data_cervus_raw <- read_csv(file_cervus, show_col_types = FALSE)

# Keep high-confidence links in either candidate column (* or +)
data_cervus <- data_cervus_raw %>%
  filter(
    Pair_confidence1 %in% c("*", "+") |
      Pair_confidence2 %in% c("*", "+")
  ) %>%
  rename(
    OffspringID = Offspring_ID,
    Parent1_Cervus = First_candidate_ID,
    Parent2_Cervus = Second_candidate_ID
  ) %>%
  select(OffspringID, Parent1_Cervus, Parent2_Cervus)

################################################################################
######################## 2) LOAD + FILTER COLONY ###############################
################################################################################

data_colony_raw <- read_csv(file_colony, show_col_types = FALSE)

# Keep records where at least one parent is not inferred (#...)
data_colony <- data_colony_raw %>%
  filter(!(str_starts(FatherID, "#") & str_starts(MotherID, "#"))) %>%
  transmute(
    OffspringID = OffspringID,
    Parent1_Colony = FatherID,
    Parent2_Colony = MotherID
  )

################################################################################
######################## 3) COMPARE CERVUS vs COLONY ###########################
################################################################################

comparison_data <- left_join(data_cervus, data_colony, by = "OffspringID")

# Concordance stats (simple exact match)
concordance <- comparison_data %>%
  mutate(
    Parent1_Match = Parent1_Cervus == Parent1_Colony,
    Parent2_Match = Parent2_Cervus == Parent2_Colony
  ) %>%
  summarise(
    Concordance_Parent1 = sum(Parent1_Match, na.rm = TRUE),
    Concordance_Parent2 = sum(Parent2_Match, na.rm = TRUE),
    Total = n()
  )

print(concordance)

# Remove rows where COLONY provides no usable parent at all
comparison_data <- comparison_data %>%
  filter(!(is.na(Parent1_Colony) & is.na(Parent2_Colony)))

################################################################################
######################## 4) LOAD HAPLOTYPES + MERGE ############################
################################################################################

# If your haplotype file is ";" separated, set delim=";".
haplo <- read_delim(
  file_haplo,
  delim = ",",
  col_types = cols(
    ID = col_character(),
    plot = col_character(),
    lat = col_double(),
    long = col_double(),
    dbh = col_double(),
    alt = col_double(),
    Haplotypes = col_character()
  )
)

# Merge offspring haplotype + coordinates, then parent1 and parent2 info from haplo
comparison_data <- comparison_data %>%
  left_join(
    haplo %>%
      transmute(
        OffspringID = ID,
        Plot_Offspring = plot,
        Latitude_Offspring = lat,
        Longitude_Offspring = long,
        DBH_Offspring = dbh,
        Haplotypes_Offspring = Haplotypes
      ),
    by = "OffspringID"
  ) %>%
  left_join(
    haplo %>%
      transmute(
        Parent1_Cervus = ID,
        Plot_Parent1 = plot,
        Latitude_Parent1 = lat,
        Longitude_Parent1 = long,
        DBH_Parent1 = dbh,
        Haplotypes_Parent1 = Haplotypes
      ),
    by = "Parent1_Cervus"
  ) %>%
  left_join(
    haplo %>%
      transmute(
        Parent2_Cervus = ID,
        Plot_Parent2 = plot,
        Latitude_Parent2 = lat,
        Longitude_Parent2 = long,
        DBH_Parent2 = dbh,
        Haplotypes_Parent2 = Haplotypes
      ),
    by = "Parent2_Cervus"
  )

# Haplotype match flags (careful with NA)
comparison_data <- comparison_data %>%
  mutate(
    Match_Parent1 = ifelse(!is.na(Haplotypes_Offspring) & !is.na(Haplotypes_Parent1) &
                             Haplotypes_Offspring == Haplotypes_Parent1, "Yes", "No"),
    Match_Parent2 = ifelse(!is.na(Haplotypes_Offspring) & !is.na(Haplotypes_Parent2) &
                             Haplotypes_Offspring == Haplotypes_Parent2, "Yes", "No")
  )

# Save full comparison table
write.csv(comparison_data, out_comparison_csv, row.names = FALSE)
cat("Saved:", out_comparison_csv, "\n")

################################################################################
######################## 5) INFER MOTHER / FATHER ##############################
################################################################################

# Keep only cases where EXACTLY one parent matches the offspring haplotype
filtered_data <- comparison_data %>%
  filter(
    (Match_Parent1 == "Yes" & Match_Parent2 == "No") |
      (Match_Parent1 == "No" & Match_Parent2 == "Yes")
  ) %>%
  mutate(
    Mother_ID = ifelse(Match_Parent1 == "Yes", Parent1_Cervus, Parent2_Cervus),
    Father_ID = ifelse(Match_Parent1 == "Yes", Parent2_Cervus, Parent1_Cervus),
    
    Latitude_Mother  = ifelse(Mother_ID == Parent1_Cervus, Latitude_Parent1, Latitude_Parent2),
    Longitude_Mother = ifelse(Mother_ID == Parent1_Cervus, Longitude_Parent1, Longitude_Parent2),
    Latitude_Father  = ifelse(Father_ID == Parent1_Cervus, Latitude_Parent1, Latitude_Parent2),
    Longitude_Father = ifelse(Father_ID == Parent1_Cervus, Longitude_Parent1, Longitude_Parent2)
  ) %>%
  mutate(
    Distance_To_Mother = dist_m(Longitude_Offspring, Latitude_Offspring, Longitude_Mother, Latitude_Mother),
    Distance_To_Father = dist_m(Longitude_Offspring, Latitude_Offspring, Longitude_Father, Latitude_Father)
  )

write.csv(filtered_data, out_filtered_csv, row.names = FALSE)
cat("Saved:", out_filtered_csv, "\n")

cat("Mean distance to mothers (m):", round(mean(filtered_data$Distance_To_Mother, na.rm = TRUE), 2), "\n")
cat("Mean distance to fathers (m):", round(mean(filtered_data$Distance_To_Father, na.rm = TRUE), 2), "\n")

################################################################################
######################## 6) MAPS (BASIC) #######################################
################################################################################

# Build long table for links
map_links <- filtered_data %>%
  select(
    OffspringID, Latitude_Offspring, Longitude_Offspring,
    Mother_ID, Latitude_Mother, Longitude_Mother,
    Father_ID, Latitude_Father, Longitude_Father
  ) %>%
  pivot_longer(
    cols = c(Latitude_Mother, Longitude_Mother, Latitude_Father, Longitude_Father),
    names_to = c(".value", "Parent_Type"),
    names_pattern = "(Latitude|Longitude)_(.*)"
  ) %>%
  mutate(
    Parent_ID = ifelse(Parent_Type == "Father", Father_ID, Mother_ID),
    Link_Type = ifelse(Parent_Type == "Mother", "Mother-Offspring", "Father-Offspring")
  )

p_map_basic <- ggplot() +
  geom_segment(
    data = map_links,
    aes(
      x = Longitude_Offspring, y = Latitude_Offspring,
      xend = Longitude, yend = Latitude, color = Link_Type
    ),
    linewidth = 0.9, alpha = 0.7
  ) +
  geom_point(
    data = filtered_data,
    aes(x = Longitude_Offspring, y = Latitude_Offspring),
    color = "forestgreen", size = 3, shape = 16, alpha = 0.85
  ) +
  geom_point(
    data = filtered_data,
    aes(x = Longitude_Mother, y = Latitude_Mother),
    color = "royalblue", size = 3, shape = 17, alpha = 0.85
  ) +
  geom_point(
    data = filtered_data,
    aes(x = Longitude_Father, y = Latitude_Father),
    color = "royalblue", size = 3, shape = 17, alpha = 0.85
  ) +
  scale_color_manual(values = c("Mother-Offspring" = "red", "Father-Offspring" = "blue")) +
  theme_minimal() +
  labs(
    title = paste("Parentage map with haplotype-based mother inference -", plot_name),
    x = "Longitude", y = "Latitude", color = "Link"
  )

print(p_map_basic)
ggsave(out_map_basic_png, p_map_basic, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved:", out_map_basic_png, "\n")

################################################################################
######################## 7) MAPS (WITH ISOLINES) ###############################
################################################################################

if (!is.na(file_isolines) && ok_sf && file.exists(file_isolines)) {
  library(sf)
  
  isolines <- st_read(file_isolines, quiet = TRUE)
  isolines_wgs84 <- st_transform(isolines, crs = 4326)
  
  # bbox from offspring points
  pts_sf <- st_as_sf(
    filtered_data,
    coords = c("Longitude_Offspring", "Latitude_Offspring"),
    crs = 4326
  )
  bbox <- st_bbox(pts_sf)
  bbox_exp <- bbox
  bbox_exp["xmin"] <- bbox["xmin"] - bbox_expand_deg
  bbox_exp["ymin"] <- bbox["ymin"] - bbox_expand_deg
  bbox_exp["xmax"] <- bbox["xmax"] + bbox_expand_deg
  bbox_exp["ymax"] <- bbox["ymax"] + bbox_expand_deg
  
  isolines_crop <- st_crop(isolines_wgs84, bbox_exp)
  
  p_map_iso <- ggplot() +
    geom_sf(data = isolines_crop, color = "gray50", linewidth = 0.3, alpha = 0.7) +
    geom_segment(
      data = map_links,
      aes(
        x = Longitude_Offspring, y = Latitude_Offspring,
        xend = Longitude, yend = Latitude, color = Link_Type
      ),
      linewidth = 0.9, alpha = 0.7
    ) +
    geom_point(
      data = filtered_data,
      aes(x = Longitude_Offspring, y = Latitude_Offspring, shape = "Offspring"),
      color = "forestgreen", size = 3, alpha = 0.85
    ) +
    geom_point(
      data = filtered_data,
      aes(x = Longitude_Mother, y = Latitude_Mother, shape = "Parent"),
      color = "royalblue", size = 3, alpha = 0.85
    ) +
    geom_point(
      data = filtered_data,
      aes(x = Longitude_Father, y = Latitude_Father, shape = "Parent"),
      color = "royalblue", size = 3, alpha = 0.85
    ) +
    scale_color_manual(values = c("Mother-Offspring" = "red", "Father-Offspring" = "blue"),
                       name = "Parent-Offspring link") +
    scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
    theme_minimal() +
    theme(
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4)
    ) +
    labs(
      title = paste("Seed & pollen dispersal inferred from haplotypes -", plot_name),
      x = "Longitude", y = "Latitude"
    ) +
    coord_sf(
      xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
      ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
      expand = FALSE
    )
  
  print(p_map_iso)
  ggsave(out_map_isoline_png, p_map_iso, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved:", out_map_isoline_png, "\n")
  
  # Mother-only map
  mother_links <- map_links %>% filter(Link_Type == "Mother-Offspring")
  
  p_mother_iso <- ggplot() +
    geom_sf(data = isolines_crop, color = "gray50", linewidth = 0.3, alpha = 0.7) +
    geom_segment(
      data = mother_links,
      aes(
        x = Longitude_Offspring, y = Latitude_Offspring,
        xend = Longitude, yend = Latitude
      ),
      color = "red", linewidth = 0.9, alpha = 0.7
    ) +
    geom_point(
      data = filtered_data,
      aes(x = Longitude_Offspring, y = Latitude_Offspring, shape = "Offspring"),
      color = "forestgreen", size = 3, alpha = 0.85
    ) +
    geom_point(
      data = filtered_data,
      aes(x = Longitude_Mother, y = Latitude_Mother, shape = "Parent"),
      color = "royalblue", size = 3, alpha = 0.85
    ) +
    scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
    theme_minimal() +
    theme(
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4)
    ) +
    labs(
      title = paste("Seed dispersal (mother-offspring only) -", plot_name),
      x = "Longitude", y = "Latitude"
    ) +
    coord_sf(
      xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
      ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
      expand = FALSE
    )
  
  print(p_mother_iso)
  ggsave(out_map_mother_png, p_mother_iso, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved:", out_map_mother_png, "\n")
}

################################################################################
######################## 8) DISTANCE DISTRIBUTIONS (PLOT) ######################
################################################################################

# --- Pollen distance histogram (10 m bins) ---
pollen_summary <- filtered_data %>%
  mutate(Pollen_Distance_Bin = floor(Distance_To_Father / 10) * 10) %>%
  filter(!is.na(Pollen_Distance_Bin)) %>%
  count(Pollen_Distance_Bin, name = "Events")

p_pollen <- ggplot(pollen_summary, aes(x = Pollen_Distance_Bin, y = Events)) +
  geom_col(fill = plot_color, color = "black") +
  theme_minimal() +
  labs(
    title = paste("Pollen dispersal distance distribution -", plot_name),
    x = "Pollen dispersal distance (m)",
    y = "Number of events"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_pollen)
ggsave(out_pollen_hist_png, p_pollen, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved:", out_pollen_hist_png, "\n")

# --- Seed distance histogram (10 m bins) ---
seed_summary <- filtered_data %>%
  mutate(Seed_Distance_Bin = floor(Distance_To_Mother / 10) * 10) %>%
  filter(!is.na(Seed_Distance_Bin)) %>%
  count(Seed_Distance_Bin, name = "Events")

p_seed <- ggplot(seed_summary, aes(x = Seed_Distance_Bin, y = Events)) +
  geom_col(fill = plot_color, color = "black") +
  theme_minimal() +
  labs(
    title = paste("Seed dispersal distance distribution -", plot_name),
    x = "Seed dispersal distance (m)",
    y = "Number of events"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_seed)
ggsave(out_seed_hist_png, p_seed, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved:", out_seed_hist_png, "\n")

################################################################################
######################## 9) OPTIONAL: ALL PLOTS COMBINED #######################
################################################################################

if (make_all_plots_combined) {
  stop_if_missing(file_filtered_SPR)
  stop_if_missing(file_filtered_REG)
  stop_if_missing(file_filtered_NOU)
  stop_if_missing(file_filtered_PAR)
  
  # If your filtered_data_* files are ";" separated, switch sep=";"
  spr <- read.csv(file_filtered_SPR, sep = ";") %>% mutate(Plot = "Sparouine")
  reg <- read.csv(file_filtered_REG, sep = ";") %>% mutate(Plot = "Regina")
  nou <- read.csv(file_filtered_NOU, sep = ";") %>% mutate(Plot = "Nouragues")
  par <- read.csv(file_filtered_PAR, sep = ";") %>% mutate(Plot = "Paracou")
  
  all_data <- bind_rows(spr, reg, nou, par)
  
  # ---- Seed global ----
  mean_seed <- mean(all_data$Distance_To_Mother, na.rm = TRUE)
  
  seed_global <- all_data %>%
    mutate(Seed_Distance_Bin = floor(Distance_To_Mother / 10) * 10) %>%
    filter(!is.na(Seed_Distance_Bin)) %>%
    count(Seed_Distance_Bin, name = "Events")
  
  p_seed_global <- ggplot(seed_global, aes(x = Seed_Distance_Bin, y = Events)) +
    geom_col(fill = "#B5B5B5", color = "black") +
    geom_vline(xintercept = mean_seed, color = "red", linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    labs(
      title = "Global seed dispersal distance distribution (all plots combined)",
      x = "Seed dispersal distance (m)",
      y = "Number of events"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_seed_global)
  ggsave(out_seed_global_png, p_seed_global, width = 10, height = 5, dpi = 300, bg = "white")
  cat("Saved:", out_seed_global_png, "\n")
  
  # ---- Pollen global ----
  mean_pollen <- mean(all_data$Distance_To_Father, na.rm = TRUE)
  
  pollen_global <- all_data %>%
    mutate(Pollen_Distance_Bin = floor(Distance_To_Father / 10) * 10) %>%
    filter(!is.na(Pollen_Distance_Bin)) %>%
    count(Pollen_Distance_Bin, name = "Events")
  
  p_pollen_global <- ggplot(pollen_global, aes(x = Pollen_Distance_Bin, y = Events)) +
    geom_col(fill = "#B5B5B5", color = "black") +
    geom_vline(xintercept = mean_pollen, color = "red", linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    labs(
      title = "Global pollen dispersal distance distribution (all plots combined)",
      x = "Pollen dispersal distance (m)",
      y = "Number of events"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_pollen_global)
  ggsave(out_pollen_global_png, p_pollen_global, width = 10, height = 5, dpi = 300, bg = "white")
  cat("Saved:", out_pollen_global_png, "\n")
}
