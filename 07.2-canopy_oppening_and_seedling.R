########################################################################
######### Link with seedling establisment and canopy openning ##########
########################################################################


# 1. Load required packages
library(sf)
library(dplyr)
library(ggplot2)

# 2. Read shapefiles
gaps <- st_read("REG74_trouées_abattage.shp")
roads <- st_read("REG74_pistes_total.shp")
landings <- st_read("REG74_places_dépot.shp")

# 3. Read tree dataset
trees <- read.csv("PAI74_tree_info.csv")

# 4. Convert to sf object
# Adapt column names if needed
trees_sf <- st_as_sf(trees, coords = c("long", "lat"), crs = 4326)
# Reproject to match the shapefiles' CRS
target_crs <- st_crs(gaps)
trees_sf <- st_transform(trees_sf, crs = target_crs)

# 5. Filter seedlings (e.g., DBH < 2 cm)
seedlings <- trees_sf %>% filter(dbh < 2)

# 6. Check intersection with each opening type
seedlings$in_gap <- lengths(st_intersects(seedlings, gaps)) > 0
seedlings$in_road <- lengths(st_intersects(seedlings, roads)) > 0
seedlings$in_landing <- lengths(st_intersects(seedlings, landings)) > 0

# 7. Combined variable: inside any opening
seedlings$in_any_opening <- seedlings$in_gap | seedlings$in_road | seedlings$in_landing

# 8. Summary tables
table_gap <- seedlings %>% group_by(in_gap) %>% summarise(count = n())
table_road <- seedlings %>% group_by(in_road) %>% summarise(count = n())
table_landing <- seedlings %>% group_by(in_landing) %>% summarise(count = n())
table_any <- seedlings %>% group_by(in_any_opening) %>% summarise(count = n())

# Display summaries
print("Seedlings inside felling gaps:"); print(table_gap)
print("Seedlings inside roads:"); print(table_road)
print("Seedlings inside landings:"); print(table_landing)
print("Seedlings inside any opening:"); print(table_any)

# 9. Proportion test (inside vs outside any opening)
prop.test(x = table_any$count, n = sum(table_any$count))

# 10. Visualization
ggplot() +
  geom_sf(data = gaps, fill = "lightgreen", alpha = 0.4) +
  geom_sf(data = roads, fill = "orange", alpha = 0.4) +
  geom_sf(data = landings, fill = "brown", alpha = 0.5) +
  geom_sf(data = seedlings, aes(color = in_any_opening), size = 1.5) +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  labs(title = "Seedling establishment inside logging-related openings",
       color = "Inside opening") +
  theme_minimal()
