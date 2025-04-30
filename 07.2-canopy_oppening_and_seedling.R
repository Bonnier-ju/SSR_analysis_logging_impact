########################################################################
######### Link with seedling establisment and canopy openning ##########
########################################################################


# Load required packages
library(sf)
library(dplyr)
library(ggplot2)
library(logistf)

# File paths
base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Echantillonnage/Regina/Données initiales/"
tree_data_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/PAI74/PAI74_tree_info.csv"

gaps_path <- paste0(base_path, "REG74_trouées_abattage.shp")
roads_path <- paste0(base_path, "REG74_pistes_total.shp")
landings_path <- paste0(base_path, "REG74_places_dépot.shp")

# Read shapefiles
gaps <- st_read(gaps_path)
roads <- st_read(roads_path)
landings <- st_read(landings_path)

# Tree data
trees <- read.csv(tree_data_path)

# Convert to sf object and reproject if needed
trees_sf <- st_as_sf(trees, coords = c("long", "lat"), crs = 4326)  # adapt column names if necessary
target_crs <- st_crs(gaps)
trees_sf <- st_transform(trees_sf, crs = target_crs)

# Filter seedlings 
seedlings <- trees_sf %>% filter(dbh < 5)

# Spatial intersection with each logging feature
seedlings$in_gap <- lengths(st_intersects(seedlings, gaps)) > 0
seedlings$in_road <- lengths(st_intersects(seedlings, roads)) > 0
seedlings$in_landing <- lengths(st_intersects(seedlings, landings)) > 0

#Combined variable: inside any canopy opening
seedlings$in_any_opening <- seedlings$in_gap | seedlings$in_road | seedlings$in_landing

#Summary tables
table_gap <- seedlings %>% group_by(in_gap) %>% summarise(count = n())
table_road <- seedlings %>% group_by(in_road) %>% summarise(count = n())
table_landing <- seedlings %>% group_by(in_landing) %>% summarise(count = n())
table_any <- seedlings %>% group_by(in_any_opening) %>% summarise(count = n())

# Print summaries
print("Seedlings inside felling gaps:"); print(table_gap)
print("Seedlings inside roads:"); print(table_road)
print("Seedlings inside landings:"); print(table_landing)
print("Seedlings inside any type of opening:"); print(table_any)

# Fit Firth's logistic regression
firth_model <- logistf(in_any_opening ~ in_gap + in_road + in_landing, data = seedlings)
summary(firth_model)

### Ploting results ###
# Get the bounding box of the seedlings
bbox <- st_bbox(seedlings)

# Plot: cropped to sampling zone
ggplot() +
  geom_sf(data = gaps, fill = "lightgreen", color = NA, alpha = 0.5) +
  geom_sf(data = roads, fill = "orange", color = NA, alpha = 0.5) +
  geom_sf(data = landings, fill = "sienna", color = NA, alpha = 0.6) +
  
  geom_sf(data = seedlings, aes(color = in_any_opening), size = 1.8, shape = 21, stroke = 0.3, fill = "white") +
  
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),
    labels = c("Outside opening", "Inside opening"),
    name = "Seedling location"
  ) +
  
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  
  labs(
    title = "Seedling establishment relative to logging-induced canopy openings",
    subtitle = "Visualization restricted to seedling sampling zone",
    x = "Easting (m)", y = "Northing (m)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )



##################### Link canopy opening and seed distances ###################
library(dplyr)
library(sf)
library(readr)


parentage_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo/filtered_data_PAI74.csv"
df <- read_csv(parentage_path)

# Convert offspring coordinates to sf object
offspring_sf <- st_as_sf(df, coords = c("Longitude_Offspring", "Latitude_Offspring"), crs = 4326)
offspring_sf <- st_transform(offspring_sf, crs = st_crs(gaps))  # same CRS as the shapefiles

# Intersect with logging features
offspring_sf$in_gap <- lengths(st_intersects(offspring_sf, gaps)) > 0
offspring_sf$in_road <- lengths(st_intersects(offspring_sf, roads)) > 0
offspring_sf$in_landing <- lengths(st_intersects(offspring_sf, landings)) > 0
offspring_sf$in_any_opening <- offspring_sf$in_gap | offspring_sf$in_road | offspring_sf$in_landing

# Wilcoxon rank-sum test
wilcox.test(Distance_To_Mother ~ in_any_opening, data = offspring_sf)

# Convert logical to factor for labeling
offspring_sf$in_any_opening <- factor(offspring_sf$in_any_opening, levels = c(FALSE, TRUE),
                                      labels = c("Outside opening", "Inside opening"))

# Enhanced boxplot
ggplot(offspring_sf, aes(x = in_any_opening, y = Distance_To_Mother, fill = in_any_opening)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, color = "black", size = 1) +
  scale_fill_manual(values = c("Outside opening" = "#1f78b4", "Inside opening" = "#e31a1c")) +
  labs(
    title = "Mother–offspring distances by canopy opening",
    x = NULL,
    y = "Distance to mother (m)",
    fill = "Seedling location"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

########### Comparing genetic diversity inside and outside canopy opening #############


