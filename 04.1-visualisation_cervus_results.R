############### Visualisation of Cervus results ###############
###############################################################

library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(geosphere)
library(readr)


# Load CERVUS result file
results <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/HKO50/results/summary_cervus_HKO50.csv", sep = ";")
geo_data <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/metadata_HKO50.csv", sep = ";")


# Merge coordinates with offspring and candidate parents
results_geo <- results %>%
  left_join(geo_data, by = c("Offspring.ID" = "ID")) %>%
  rename(lat_offspring = lat, long_offspring = long) %>%
  left_join(geo_data, by = c("First.candidate.ID" = "ID")) %>%
  rename(lat_first_parent = lat, long_first_parent = long) %>%
  left_join(geo_data, by = c("Second.candidate.ID" = "ID")) %>%
  rename(lat_second_parent = lat, long_second_parent = long)

# Extract high-confidence pair links
links_first <- results_geo %>%
  filter(Pair.confidence %in% c("*", "+")) %>%
  select(Offspring.ID, lat_offspring, long_offspring,
         First.candidate.ID, lat_first_parent, long_first_parent, Pair.confidence) %>%
  rename(parent_id = First.candidate.ID,
         lat_parent = lat_first_parent,
         long_parent = long_first_parent,
         confidence = Pair.confidence)

links_second <- results_geo %>%
  filter(Pair.confidence.1 %in% c("*", "+")) %>%
  select(Offspring.ID, lat_offspring, long_offspring,
         Second.candidate.ID, lat_second_parent, long_second_parent, Pair.confidence.1) %>%
  rename(parent_id = Second.candidate.ID,
         lat_parent = lat_second_parent,
         long_parent = long_second_parent,
         confidence = Pair.confidence.1)

links <- bind_rows(links_first, links_second)


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
  labs(title = "Parent-offspring links – HKO50",
       x = "Longitude", y = "Latitude", color = "Type", shape = "Type")

# Plot results with topographic line 
