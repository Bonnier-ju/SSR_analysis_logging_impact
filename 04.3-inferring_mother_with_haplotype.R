##############################################################################################
############ Comparing Cervus and Colony results + inferring mother with haplotype ###########
##############################################################################################


library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(utils)
library(geosphere)

# Load and prepare Cervus data
file_path_cervus <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.1-parentage_cervus/PAI74/summary_PAI74.csv"
data_cervus <- read_csv(file_path_cervus)  # Assuming ';' as delimiter
data_cervus <- data_cervus %>%
  filter(grepl("\\+|\\*", Pair_confidence1) | grepl("\\+|\\*", Pair_confidence2)) %>%
  rename(OffspringID = Offspring_ID) %>%  # Rename to match Colony data
  select(OffspringID, Parent1_Cervus = First_candidate_ID, Parent2_Cervus = Second_candidate_ID)

# Load and prepare Colony data
file_path_colony <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.2-parentage_Colony/PAI74/results_colony_PAI74/PAI74_colony_projet.BestConfig.csv"
data_colony <- read_csv(file_path_colony) %>%
  filter(!(str_starts(FatherID, "#") & str_starts(MotherID, "#"))) %>%  # Exclude only if both parents are supposed
  select(OffspringID, Parent1_Colony = FatherID, Parent2_Colony = MotherID)

# Merge Cervus and Colony data
comparison_data <- left_join(data_cervus, data_colony, by = "OffspringID")

# Analyze the concordance
comparison_data %>%
  mutate(Parent1_Match = Parent1_Cervus == Parent1_Colony,
         Parent2_Match = Parent2_Cervus == Parent2_Colony) %>%
  summarise(Concordance_Parent1 = sum(Parent1_Match, na.rm = TRUE),
            Concordance_Parent2 = sum(Parent2_Match, na.rm = TRUE),
            Total = n())


# Remove rows where both Parent1_Colony and Parent2_Colony are NA
comparison_data <- comparison_data %>%
  filter(!(is.na(Parent1_Colony) & is.na(Parent2_Colony)))

# Print the cleaned comparison data
print(comparison_data)



# Load haplotype data
file_path_haplotypes <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/Chloroplatics_markers/Haplotype_full.csv"

# Load haplotype data using read.table
data_haplotypes <- read.table(file_path_haplotypes, 
                              sep = ",",               # Set semicolon as the field separator
                              dec = ".",               # Set period as the decimal mark
                              header = TRUE,           # Assumes the first line in your file has column headers
                              stringsAsFactors = FALSE, # Ensures text data isn't converted to factors
                              colClasses = c(ID="character", 
                                             plot="character",
                                             lat="numeric",
                                             long="numeric",
                                             dbh="numeric",
                                             alt="numeric",
                                             Haplotypes="character") # Explicitly define column classes
)

# Merge haplotype data into comparison_data for offspring and both parents
comparison_data <- comparison_data %>%
  left_join(data_haplotypes, by = c("OffspringID" = "ID")) %>%
  left_join(data_haplotypes %>% rename(Haplotypes_Parent1 = Haplotypes), by = c("Parent1_Cervus" = "ID")) %>%
  left_join(data_haplotypes %>% rename(Haplotypes_Parent2 = Haplotypes), by = c("Parent2_Cervus" = "ID"))

# Add columns to indicate haplotype match with each parent
comparison_data <- comparison_data %>%
  mutate(
    Match_Parent1 = ifelse(Haplotypes == Haplotypes_Parent1, "Yes", "No"),
    Match_Parent2 = ifelse(Haplotypes == Haplotypes_Parent2, "Yes", "No")
  )

# Display the results with haplotype match information
print(comparison_data %>% select(OffspringID, Parent1_Cervus, Parent2_Cervus, Match_Parent1, Match_Parent2, Haplotypes, Haplotypes_Parent1, Haplotypes_Parent2))

print(comparison_data)

# Clean up the comparison_data by properly structuring the information for offspring, Parent1, and Parent2
comparison_data <- comparison_data %>%
  select(
    # Offspring Information
    OffspringID,
    Plot_Offspring = plot.x,
    Latitude_Offspring = lat.x,
    Longitude_Offspring = long.x,
    DBH_Offspring = dbh.x,
    Haplotypes_Offspring = Haplotypes,
    
    # Parent 1 Information
    Parent1_Cervus,
    Plot_Parent1 = plot.y,
    Latitude_Parent1 = lat.y,
    Longitude_Parent1 = long.y,
    DBH_Parent1 = dbh.y,
    Haplotypes_Parent1,
    
    # Parent 2 Information
    Parent2_Cervus,
    Plot_Parent2 = plot,
    Latitude_Parent2 = lat,
    Longitude_Parent2 = long,
    DBH_Parent2 = dbh,
    Haplotypes_Parent2,
    
    # Haplotype Matching Results
    Match_Parent1,
    Match_Parent2
  ) %>%
  distinct()  # Remove any duplicate rows if they exist

# Print the cleaned and structured data
print(comparison_data)

#Save csv file
output_file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo/comparison_data_PAI74.csv"
write.csv(comparison_data, file = output_file_path, row.names = FALSE)


# Create a new dataset where exactly one parent matches the offspring's haplotype
filtered_data <- comparison_data %>%
  filter((Match_Parent1 == "Yes" & Match_Parent2 == "No") | 
           (Match_Parent1 == "No" & Match_Parent2 == "Yes")) %>%
  
  # Assign mother and father based on haplotype match
  mutate(
    Mother_ID = ifelse(Match_Parent1 == "Yes", Parent1_Cervus, Parent2_Cervus),
    Father_ID = ifelse(Match_Parent1 == "No", Parent1_Cervus, Parent2_Cervus)
  )


# Print the filtered dataset with maternal and paternal assignments
print(filtered_data)



# Function to calculate Haversine distance (in meters)
calculate_distance <- function(lat1, lon1, lat2, lon2) {
  if (is.na(lat1) | is.na(lon1) | is.na(lat2) | is.na(lon2)) {
    return(NA)  # Return NA if any coordinate is missing
  }
  return(distHaversine(c(lon1, lat1), c(lon2, lat2)))  # Returns distance in meters
}

# Compute the correct latitude and longitude for the mother and father
filtered_data <- filtered_data %>%
  mutate(
    # Assign correct coordinates for Mother
    Latitude_Mother = ifelse(Mother_ID == Parent1_Cervus, Latitude_Parent1, Latitude_Parent2),
    Longitude_Mother = ifelse(Mother_ID == Parent1_Cervus, Longitude_Parent1, Longitude_Parent2),
    
    # Assign correct coordinates for Father
    Latitude_Father = ifelse(Father_ID == Parent1_Cervus, Latitude_Parent1, Latitude_Parent2),
    Longitude_Father = ifelse(Father_ID == Parent1_Cervus, Longitude_Parent1, Longitude_Parent2)
  ) %>%
  
  # Calculate distances
  rowwise() %>%
  mutate(
    Distance_To_Mother = calculate_distance(Latitude_Offspring, Longitude_Offspring, Latitude_Mother, Longitude_Mother),
    Distance_To_Father = calculate_distance(Latitude_Offspring, Longitude_Offspring, Latitude_Father, Longitude_Father)
  ) %>%
  ungroup()

# Compute the mean distance separately for mothers and fathers
mean_distance_mother <- mean(filtered_data$Distance_To_Mother, na.rm = TRUE)
mean_distance_father <- mean(filtered_data$Distance_To_Father, na.rm = TRUE)

# Print results
print(paste("Mean Distance to Mothers:", round(mean_distance_mother, 2), "meters"))
print(paste("Mean Distance to Fathers:", round(mean_distance_father, 2), "meters"))

# Define the file path for saving filtered data
output_filtered_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo/filtered_data_HKO50.csv"
write.csv(filtered_data, file = output_filtered_path, row.names = FALSE)






########### Plotting results on a map ##########
################################################


library(ggplot2)
library(ggrepel)
library(tidyr)

# Prepare the dataset using pivot_longer()
map_data <- filtered_data %>%
  select(OffspringID, Latitude_Offspring, Longitude_Offspring, 
         Mother_ID, Latitude_Mother, Longitude_Mother,
         Father_ID, Latitude_Father, Longitude_Father) %>%
  pivot_longer(cols = c(Latitude_Mother, Longitude_Mother, Latitude_Father, Longitude_Father),
               names_to = c(".value", "Parent_Type"),
               names_pattern = "(Latitude|Longitude)_(.*)") %>%
  mutate(Parent_ID = ifelse(Parent_Type == "Father", Father_ID, Mother_ID),
         Link_Type = ifelse(Parent_Type == "Mother", "Mother-Offspring", "Father-Offspring"))

# Create the map
ggplot() +
  # Plot parent-offspring connections
  geom_segment(data = map_data, aes(x = Longitude_Offspring, y = Latitude_Offspring, 
                                    xend = Longitude, yend = Latitude, 
                                    color = Link_Type), size = 1, alpha = 0.7) +
  
  # Plot offspring points (Green Circles)
  geom_point(data = filtered_data, aes(x = Longitude_Offspring, y = Latitude_Offspring), 
             color = "green", size = 3, shape = 16, alpha = 0.8) +  # Shape 16 = Circle
  
  # Plot parent points (Blue Triangles)
  geom_point(data = filtered_data, aes(x = Longitude_Mother, y = Latitude_Mother), 
             color = "blue", size = 3, shape = 17, alpha = 0.8) +  # Shape 17 = Triangle
  geom_point(data = filtered_data, aes(x = Longitude_Father, y = Latitude_Father), 
             color = "blue", size = 3, shape = 17, alpha = 0.8) +  # Shape 17 = Triangle
  
  # Customize legend and theme
  scale_color_manual(values = c("Mother-Offspring" = "red", "Father-Offspring" = "blue")) +
  labs(title = "Parentage Map of Offspring, Mothers, and Fathers",
       x = "Longitude", y = "Latitude", color = "Link Type") +
  theme_minimal()




########### Map of results with isolignes ############
######################################################

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(ggrepel)

# Load the isoline shapefile
file_path_isolines <- "C:/Users/bonni/Desktop/Fichiers_cartes_Qgis/Isolignes/Isolignes_Regina_5m/fr_662043116_lidar_regstgmult2013_02.shp"
isolines <- st_read(file_path_isolines)

# Transform the isolines to WGS84 (EPSG:4326) to match the offspring-parent data
isolines_wgs84 <- st_transform(isolines, crs = 4326)

# Convert offspring-parent data into sf object
filtered_data_sf <- st_as_sf(filtered_data, coords = c("Longitude_Offspring", "Latitude_Offspring"), crs = 4326)

# Compute bounding box for the study area
bbox_coords <- st_bbox(filtered_data_sf)

# Expand bounding box slightly
xmin <- bbox_coords$xmin - 0.001
ymin <- bbox_coords$ymin - 0.001
xmax <- bbox_coords$xmax + 0.001
ymax <- bbox_coords$ymax + 0.001

# Prepare the dataset for visualization
map_data <- filtered_data %>%
  select(OffspringID, Latitude_Offspring, Longitude_Offspring, 
         Mother_ID, Latitude_Mother, Longitude_Mother,
         Father_ID, Latitude_Father, Longitude_Father) %>%
  pivot_longer(cols = c(Latitude_Mother, Longitude_Mother, Latitude_Father, Longitude_Father),
               names_to = c(".value", "Parent_Type"),
               names_pattern = "(Latitude|Longitude)_(.*)") %>%
  mutate(Parent_ID = ifelse(Parent_Type == "Father", Father_ID, Mother_ID),
         Link_Type = ifelse(Parent_Type == "Mother", "Mother-Offspring", "Father-Offspring"))

# Create the final map with legend inside the plot
# 1000X800

ggplot() +
  # Add isolines
  geom_sf(data = isolines_wgs84, color = "gray50", size = 0.3, alpha = 0.7) +
  
  # Plot parent-offspring connections
  geom_segment(data = map_data, aes(x = Longitude_Offspring, y = Latitude_Offspring, 
                                    xend = Longitude, yend = Latitude, 
                                    color = Link_Type), size = 1, alpha = 0.7) +
  
  # Plot offspring points (Green Circles)
  geom_point(data = filtered_data, aes(x = Longitude_Offspring, y = Latitude_Offspring, shape = "Offspring"), 
             color = "green", size = 3, alpha = 0.8) +  # Shape 16 = Circle
  
  # Plot parent points (Blue Triangles)
  geom_point(data = filtered_data, aes(x = Longitude_Mother, y = Latitude_Mother, shape = "Parent"), 
             color = "blue", size = 3, alpha = 0.8) +  # Shape 17 = Triangle
  geom_point(data = filtered_data, aes(x = Longitude_Father, y = Latitude_Father, shape = "Parent"), 
             color = "blue", size = 3, alpha = 0.8) +  # Shape 17 = Triangle
  
  # Customize legend (inside the plot)
  scale_color_manual(values = c("Mother-Offspring" = "red", "Father-Offspring" = "blue"), name = "Parent-Offspring Link") +
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
  
  # Position the legend inside the plot
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85),  # Bottom right inside the plot
        legend.background = element_rect(fill = "white", color = "black", size = 0.5),
        legend.key = element_rect(fill = "white")) +
  
  labs(title = "Seed and Pollen dispersal in HKO50",
       x = "Longitude", y = "Latitude") +
  
  # Set zoom limits
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  )


# === Filter only Mother-Offspring links ===
mother_links <- map_data %>%
  filter(Link_Type == "Mother-Offspring")

# Create the map with only mother-offspring connections
ggplot() +
  # Add isolines
  geom_sf(data = isolines_wgs84, color = "gray50", size = 0.3, alpha = 0.7) +
  
  # Plot mother-offspring connections
  geom_segment(data = mother_links, aes(x = Longitude_Offspring, y = Latitude_Offspring, 
                                        xend = Longitude, yend = Latitude),
               color = "red", size = 1, alpha = 0.7) +
  
  # Plot offspring points
  geom_point(data = filtered_data, aes(x = Longitude_Offspring, y = Latitude_Offspring, shape = "Offspring"), 
             color = "green", size = 3, alpha = 0.8) +
  
  # Plot mother points
  geom_point(data = filtered_data, aes(x = Longitude_Mother, y = Latitude_Mother, shape = "Parent"), 
             color = "blue", size = 3, alpha = 0.8) +
  
  # Customize legend
  scale_shape_manual(values = c("Offspring" = 16, "Parent" = 17), name = "Individuals") +
  
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5),
        legend.key = element_rect(fill = "white")) +
  
  labs(title = "Seed Dispersal (Mother-Offspring Links Only) – PAI74",
       x = "Longitude", y = "Latitude") +
  
  # Zoom on the study area
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  )



###################### New visualisation of the results #######################
###############################################################################

library(ggplot2)
library(dplyr)
library(readr)


# Loading files
filtered_data <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo/filtered_data_HKO50.csv")  


########################### Pollen dispersal bar chart #########################

# Create 10-meter bins for Distance_To_Father
filtered_data <- filtered_data %>%
  mutate(Pollen_Distance_Bin = floor(Distance_To_Father / 10) * 10)

# Summarize the number of events per bin
pollen_summary <- filtered_data %>%
  group_by(Pollen_Distance_Bin) %>%
  summarise(Events = n())

# plot 800x500
ggplot(pollen_summary, aes(x = Pollen_Distance_Bin, y = Events)) +
  geom_bar(stat = "identity", fill = "#CDAD00", color = "black") +
  scale_x_continuous(breaks = seq(0, 500, by = 50), limits = c(0, 500)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), limits = c(0, 8)) +
  theme_minimal() +
  labs(title = "Pollen Dispersal Distance Distribution - HKO50",
       x = "Pollen Dispersal Distance (m)",
       y = "Number of Events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################### Seed dispersal bar chart #########################


# Create 10-meter bins for Distance_To_Mother (seed dispersal)
filtered_data <- filtered_data %>%
  mutate(Seed_Distance_Bin = floor(Distance_To_Mother / 10) * 10)

# Summarize the number of events per bin
seed_summary <- filtered_data %>%
  group_by(Seed_Distance_Bin) %>%
  summarise(Events = n())

# Plot the histogram : 800x500
ggplot(seed_summary, aes(x = Seed_Distance_Bin, y = Events)) +
  geom_bar(stat = "identity", fill = "#CDAD00", color = "black") +
  scale_x_continuous(breaks = seq(0, 440, by = 20), limits = c(0, 440)) +
  scale_y_continuous(breaks = seq(0, 20, by = 2), limits = c(0, 20)) +
  theme_minimal() +
  labs(title = "Seed Dispersal Distance Distribution - Regina",
       x = "Seed Dispersal Distance (m)",
       y = "Number of Events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#################### Plots with four plot data ###########################

###### Seed dispersal 

# Load the 4 filtered datasets
spr <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_SPR.csv", sep=";") %>% mutate(Plot = "Sparouine")
reg <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_REG.csv", sep=";") %>% mutate(Plot = "Regina")
nou <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_NOU.csv", sep=";") %>% mutate(Plot = "Nouragues")
par <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype/filtered_data_PAR.csv", sep=";") %>% mutate(Plot = "Paracou")


# Combine all datasets into one
all_data <- bind_rows(spr, reg, nou, par)

# Create 10m bins for seed dispersal distances
all_data <- all_data %>%
  mutate(Seed_Distance_Bin = floor(Distance_To_Mother / 10) * 10)

mean_seed_distance <- mean(all_data$Distance_To_Mother, na.rm = TRUE)


# Summarize the number of events per bin (no plot distinction)
seed_summary <- all_data %>%
  group_by(Seed_Distance_Bin) %>%
  summarise(Events = n(), .groups = "drop")

# Plot 1000x500 
ggplot(seed_summary, aes(x = Seed_Distance_Bin, y = Events)) +
  geom_bar(stat = "identity", fill = "#B5B5B5", color = "black") +
  geom_vline(xintercept = mean_seed_distance, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(0, 500, by = 50), limits = c(0, 500)) +
  scale_y_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 40)) +
  theme_minimal() +
  labs(title = "Global Seed Dispersal Distance Distribution",
       x = "Seed Dispersal Distance (m)",
       y = "Number of Events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




###### Pollen dispersal 

# Create 10m bins for pollen dispersal distances
all_data <- all_data %>%
  mutate(Pollen_Distance_Bin = floor(Distance_To_Father / 10) * 10)

mean_pollen_distance <- mean(all_data$Distance_To_Father, na.rm = TRUE)


# Summarize the number of events per bin (no plot distinction)
pollen_summary <- all_data %>%
  group_by(Pollen_Distance_Bin) %>%
  summarise(Events = n(), .groups = "drop")

# Plot
ggplot(pollen_summary, aes(x = Pollen_Distance_Bin, y = Events)) +
  geom_bar(stat = "identity", fill = "#B5B5B5", color = "black") +
  geom_vline(xintercept = mean_pollen_distance, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(0, 600, by = 50), limits = c(0, 700)) +
  theme_minimal() +
  labs(title = "Global Pollen Dispersal Distance Distribution (All Plots Combined)",
       x = "Pollen Dispersal Distance (m)",
       y = "Number of Events") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



