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



