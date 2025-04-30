##############################################################
############ Visualization of Colony results #################
##############################################################

# Charger les biblioth?ques
library(ggplot2)
library(dplyr)
library(stringr)

# Charger les fichiers
file_path_1 <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.2-parentage_Colony/PAI74/results_colony_PAI74/PAI74_colony_projet.BestConfig.csv"
file_path_2 <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/PAI74/PAI74_tree_info.csv"


# Lire les fichiers CSV
best_config <- read.csv(file_path_1, header = TRUE, stringsAsFactors = FALSE, sep = ",")
tree_info <- read.csv(file_path_2, header = TRUE, stringsAsFactors = FALSE)

# Merges data (OffspringID)
offspring_positions <- best_config %>%
  left_join(tree_info, by = c("OffspringID" = "ID"))

tree_info <- tree_info %>% distinct(ID, .keep_all = TRUE)


# Merge father ID 
father_positions <- best_config %>%
  left_join(tree_info, by = c("FatherID" = "ID")) %>%
  rename(FatherLat = lat, FatherLong = long)

# Merge MotherID
mother_positions <- best_config %>%
  left_join(tree_info, by = c("MotherID" = "ID")) %>%
  rename(MotherLat = lat, MotherLong = long)


# Merge in one table
map_data <- offspring_positions %>%
  select(OffspringID, lat, long, FatherID, MotherID) %>% # Ajouter FatherID et MotherID
  rename(OffspringLat = lat, OffspringLong = long) %>%
  bind_cols(father_positions %>% select(FatherLat, FatherLong)) %>%
  bind_cols(mother_positions %>% select(MotherLat, MotherLong))

head(map_data)

#Identifier les parents "réels" (qui ne commencent pas par "#")
map_data <- map_data %>%
  mutate(
    has_real_father = !str_starts(FatherID, "#") & FatherID != "",
    has_real_mother = !str_starts(MotherID, "#") & MotherID != ""
  )

# Offspring avec un seul parent identifié
n_with_only_one_real_parent <- sum(xor(map_data$has_real_father, map_data$has_real_mother))

# Offspring avec les deux parents identifiés
n_with_both_real_parents <- sum(map_data$has_real_father & map_data$has_real_mother)

# Affichage
cat("✅ Nombre d'offspring avec UN SEUL parent identifié :", n_with_only_one_real_parent, "\n")
cat("✅ Nombre d'offspring avec les DEUX parents identifiés :", n_with_both_real_parents, "\n")




# Ploting results on a map 
# 700x700
ggplot() +
  # Ajouter les liens p?re-descendant (gris, moins transparent)
  geom_segment(data = map_data %>% filter(!is.na(FatherLat)),
               aes(x = OffspringLong, y = OffspringLat, xend = FatherLong, yend = FatherLat),
               color = "azure4", alpha = 0.9, size = 0.8) +
  
  # Ajouter les liens m?re-descendant (gris, moins transparent)
  geom_segment(data = map_data %>% filter(!is.na(MotherLat)),
               aes(x = OffspringLong, y = OffspringLat, xend = MotherLong, yend = MotherLat),
               color = "azure4", alpha = 0.9, size = 0.8) +
  
  # Ajouter les points pour descendants (ronds verts, avec l?gende)
  geom_point(data = map_data,
             aes(x = OffspringLong, y = OffspringLat, fill = "Offspring"),
             shape = 21, size = 3) +
  
  # Ajouter les points pour parents (triangles bleus, avec l?gende)
  geom_point(data = map_data %>% filter(!is.na(FatherLat)),
             aes(x = FatherLong, y = FatherLat, fill = "Parent"),
             shape = 24, size = 3) +
  geom_point(data = map_data %>% filter(!is.na(MotherLat)),
             aes(x = MotherLong, y = MotherLat, fill = "Parent"),
             shape = 24, size = 3) +
  
  # Personnalisation
  scale_fill_manual(values = c("Offspring" = "green", "Parent" = "blue"),
                    name = "Type") +
  labs(title = "PAI74",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right")






################# Map of family clusters #####################
##############################################################

# Charger les biblioth?ques
library(ggplot2)
library(dplyr)

file_path_1 <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/04-parentage_analysis/04.2-parentage_Colony/PAI74/results_colony_PAI74/PAI74_colony_projet.BestCluster.csv"
best_cluster <- read.csv(file_path_1, header = TRUE, stringsAsFactors = FALSE)
head(best_cluster)

file_path_2 <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/PAI74/PAI74_tree_info.csv"
tree_info <- read.csv(file_path_2, header = TRUE, stringsAsFactors = FALSE)

# Remplacer les "#" par NA dans FatherID et MotherID
best_cluster$FatherID <- gsub("^#", NA, best_cluster$FatherID)
best_cluster$MotherID <- gsub("^#", NA, best_cluster$MotherID)

# Filtrer les clusters ayant au moins 10 individus
cluster_sizes <- best_cluster %>%
  group_by(ClusterIndex) %>%
  summarize(total_count = n())

valid_clusters <- cluster_sizes %>%
  filter(total_count >= 10) %>%
  pull(ClusterIndex)

filtered_best_cluster <- best_cluster %>%
  filter(ClusterIndex %in% valid_clusters)

# Fusionner les donn?es pour ajouter les coordonn?es des OffspringID
offspring_positions <- filtered_best_cluster %>%
  left_join(tree_info, by = c("OffspringID" = "ID"))

# Ajouter les coordonn?es des FatherID
father_positions <- filtered_best_cluster %>%
  left_join(tree_info, by = c("FatherID" = "ID")) %>%
  rename(FatherLat = lat, FatherLong = long)

# Ajouter les coordonn?es des MotherID
mother_positions <- filtered_best_cluster %>%
  left_join(tree_info, by = c("MotherID" = "ID")) %>%
  rename(MotherLat = lat, MotherLong = long)

# Combiner les coordonn?es dans un seul tableau
cluster_data <- offspring_positions %>%
  bind_cols(father_positions %>% select(FatherLat, FatherLong)) %>%
  bind_cols(mother_positions %>% select(MotherLat, MotherLong))

# Cr?er une carte avec ggplot2
ggplot() +
  # Ajouter les liens p?re-descendant (sans l?gende)
  geom_segment(data = cluster_data %>% filter(!is.na(FatherLat)),
               aes(x = long, y = lat, xend = FatherLong, yend = FatherLat, color = factor(ClusterIndex)),
               alpha = 0.6, size = 0.8, show.legend = FALSE) +
  
  # Ajouter les liens m?re-descendant (sans l?gende)
  geom_segment(data = cluster_data %>% filter(!is.na(MotherLat)),
               aes(x = long, y = lat, xend = MotherLong, yend = MotherLat, color = factor(ClusterIndex)),
               alpha = 0.6, size = 0.8, show.legend = FALSE) +
  
  # Ajouter les descendants (losanges remplis de leur couleur de cluster)
  geom_point(data = cluster_data,
             aes(x = long, y = lat, fill = factor(ClusterIndex)),
             shape = 23, size = 3, color = "black") +
  
  # Ajouter les p?res (triangles remplis de leur couleur de cluster)
  geom_point(data = cluster_data %>% filter(!is.na(FatherLat)),
             aes(x = FatherLong, y = FatherLat, fill = factor(ClusterIndex)),
             shape = 24, size = 3, color = "black") +
  
  # Ajouter les m?res (triangles remplis de leur couleur de cluster)
  geom_point(data = cluster_data %>% filter(!is.na(MotherLat)),
             aes(x = MotherLong, y = MotherLat, fill = factor(ClusterIndex)),
             shape = 24, size = 3, color = "black") +
  
  # Personnalisation des couleurs et axes
  scale_fill_manual(values = scales::hue_pal()(length(unique(cluster_data$ClusterIndex)))) +
  scale_color_manual(values = scales::hue_pal()(length(unique(cluster_data$ClusterIndex)))) +
  labs(title = "Clusters of relations Parent-Offspring (PAI74)",
       x = "Longitude", y = "Latitude", fill = "Cluster") +  # L?gende uniquement pour les points
  theme_minimal() +
  theme(legend.position = "right")





################# Distances links #####################
#######################################################

# Fonction pour calculer la distance haversine en m?tres
haversine_distance <- function(lat1, lon1, lat2, lon2) {
  R <- 6371 * 1000 # Rayon de la Terre en m?tres
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c # Distance en m?tres
}

# Calculer les distances p?re-descendant
map_data <- map_data %>%
  mutate(FatherDistance = ifelse(!is.na(FatherLat),
                                 haversine_distance(OffspringLat, OffspringLong, FatherLat, FatherLong),
                                 NA))

# Calculer les distances m?re-descendant
map_data <- map_data %>%
  mutate(MotherDistance = ifelse(!is.na(MotherLat),
                                 haversine_distance(OffspringLat, OffspringLong, MotherLat, MotherLong),
                                 NA))

# Distance minimale entre p?re et m?re
map_data <- map_data %>%
  mutate(ParentDistance = pmin(FatherDistance, MotherDistance, na.rm = TRUE))

# Calculer le nombre total de liens parent-enfant
total_links <- sum(!is.na(map_data$FatherDistance)) + sum(!is.na(map_data$MotherDistance))

# R?sum? des distances et nombre total de liens
distance_summary <- map_data %>%
  summarize(
    MeanDistance = mean(ParentDistance, na.rm = TRUE),
    MinDistance = min(ParentDistance, na.rm = TRUE),
    MaxDistance = max(ParentDistance, na.rm = TRUE)
  ) %>%
  mutate(TotalLinks = total_links)

# Afficher le r?sum? des distances et des liens
print(distance_summary)






################# Relation DBH et nombre d'offspring #####################
#########################################################################
library(tidyr)


# Combiner les informations sur les parents dans un seul tableau
parent_offspring <- map_data %>%
  select(FatherID, MotherID) %>%
  pivot_longer(cols = c(FatherID, MotherID), names_to = "ParentType", values_to = "ParentID") %>%
  filter(!is.na(ParentID)) %>%
  group_by(ParentID) %>%
  summarize(NumOffspring = n()) %>%
  left_join(tree_info, by = c("ParentID" = "ID"))

# V?rifier les donn?es
head(parent_offspring)

# R?gression lin?aire : NumOffspring ~ DBH
model <- lm(NumOffspring ~ dbh, data = parent_offspring)

# R?sum? du mod?le
summary(model)

# Visualisation de la relation
#700 x 500
ggplot(parent_offspring, aes(x = dbh, y = NumOffspring)) +
  geom_point(size = 3, alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", color = "red", fill = "lightpink", se = TRUE) +
  labs(title = "Regina",
       x = "DBH",
       y = "Offspring Number") +
  theme_minimal()

# Save CSV
output_path <- "C:/Users/bonni/OneDrive/Universit?/Th?se/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.4-parentage_with_colony/Regina/Fig/parent_offspring.csv"
write.csv(parent_offspring, output_path, row.names = FALSE)






###### Stats on relationship results #######
############################################

# Load necessary libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(geosphere)
library(tidyr)



# Read the CSV file (make sure to adjust the path)
file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.4-parentage_with_colony/Regina/Results_Regina/Regina_Colony.BestConfig.csv"
parent_data <- read.csv(file_path, sep=",", header=TRUE, stringsAsFactors=FALSE)
# Load the geographic data
geo_file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Data_initial/Regina/geo_inds_REG.csv"
geo_data <- read.csv(geo_file_path, sep=",")
head(parent_data)


########### Sampled and inferred parents###########
# Total number of offspring
total_offspring <- nrow(parent_data)
total_offspring


# Identify inferred parents (starting with "#") for both FatherID and MotherID
parent_data$Father_status <- ifelse(str_detect(parent_data$FatherID, "^#"), "Inferred", "Sampled")
parent_data$Mother_status <- ifelse(str_detect(parent_data$MotherID, "^#"), "Inferred", "Sampled")

# Count the total number of unique parents (combining FatherID and MotherID)
total_parents <- length(unique(c(parent_data$FatherID, parent_data$MotherID)))
total_parents

# Count total inferred and sampled parents (without duplication)
inferred_parents <- length(unique(c(parent_data$FatherID[parent_data$Father_status == "Inferred"], 
                                    parent_data$MotherID[parent_data$Mother_status == "Inferred"])))

sampled_parents <- length(unique(c(parent_data$FatherID[parent_data$Father_status == "Sampled"], 
                                   parent_data$MotherID[parent_data$Mother_status == "Sampled"])))

# Count offspring linked to inferred or sampled parents
offspring_with_inferred_parents <- sum(parent_data$Father_status == "Inferred") + 
  sum(parent_data$Mother_status == "Inferred")

offspring_with_sampled_parents <- sum(parent_data$Father_status == "Sampled") + 
  sum(parent_data$Mother_status == "Sampled")

# Create a summary table
results_table <- data.frame(
  Category = c("Total Offspring",
               "Total Parents",
               "Total Inferred Parents", "Total Sampled Parents",
               "Offspring with Inferred Parents",
               "Offspring with Sampled Parents"),
  
  Count = c(total_offspring,
            total_parents,
            inferred_parents, sampled_parents,
            offspring_with_inferred_parents,
            offspring_with_sampled_parents)
)

# Print results
print(results_table)


############ Offsprings from selfing ###########
################################################

# Identify offspring from selfing (where FatherID == MotherID)
selfing_data <- parent_data %>%
  filter(FatherID == MotherID & FatherID != "")  # Ensure valid parent names

# Count the number of selfed offspring per unique parent
selfing_summary <- selfing_data %>%
  group_by(FatherID) %>%
  summarise(Offspring_Count = n(), .groups = "drop")

# Merge geographic coordinates for parents
selfing_geo <- selfing_data %>%
  left_join(geo_data, by = c("FatherID" = "ID")) %>%
  rename(Parent_lat = lat, Parent_long = long) %>%
  left_join(geo_data, by = c("OffspringID" = "ID")) %>%
  rename(Offspring_lat = lat, Offspring_long = long)

# Calculate distance between parent and offspring
selfing_geo <- selfing_geo %>%
  mutate(Distance_m = ifelse(!is.na(Parent_lat) & !is.na(Offspring_lat),
                             distHaversine(cbind(Parent_long, Parent_lat),
                                           cbind(Offspring_long, Offspring_lat)),
                             NA))

# Calculate average distance per selfing parent
selfing_avg_distance <- selfing_geo %>%
  group_by(FatherID) %>%
  summarise(Average_Distance_m = mean(Distance_m, na.rm = TRUE), .groups = "drop")

# Merge offspring count and distance data
selfing_results <- selfing_summary %>%
  left_join(selfing_avg_distance, by = "FatherID") %>%
  arrange(desc(Offspring_Count))

print(selfing_results)


# Create merged graph with dual visualization
#1300X600
ggplot(selfing_results, aes(x = reorder(FatherID, Offspring_Count))) +
  # Bar chart for number of offspring
  geom_bar(aes(y = Offspring_Count, fill = Offspring_Count), stat = "identity") +
  geom_text(aes(y = Offspring_Count, label = Offspring_Count), color = "white", hjust = 1.2, fontface = "bold", size = 5) +
  
  # Scatter plot for average distance (shown as points on the same axis)
  geom_point(aes(y = Average_Distance_m / 10), color = "red", size = 3) +  # Scale distance for visibility
  geom_text(aes(y = Average_Distance_m / 10, label = ifelse(!is.na(Average_Distance_m), 
                                                            paste0(round(Average_Distance_m, 2), " m"), "")), 
            hjust = -0.2, size = 4, fontface = "italic", color = "red") +
  
  # Formatting
  coord_flip() +
  labs(title = "Selfing Analysis: Offspring Count & Average Distance in Regina",
       x = "Parent ID",
       y = "Number of Offspring (bars) / Average Distance (scaled points)") +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "mediumorchid4") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))



############# Unique parent pairs ############
# Count the number of offspring for each unique parent pair
sibling_groups <- parent_data%>%
  group_by(FatherID, MotherID) %>%
  summarise(Offspring_Count = n(), .groups = "drop") %>%
  arrange(desc(Offspring_Count))  # Sort by highest number of offspring

# Display the first few rows of the result
print(sibling_groups)

# Create a label combining Father and Mother IDs
sibling_groups$Parent_Pair <- paste(sibling_groups$FatherID, sibling_groups$MotherID, sep = " & ")

# Select top 20 parent pairs with the most offspring for better visualization
top_sibling_groups <- sibling_groups %>%
  slice_max(order_by = Offspring_Count, n = 18)

# Plot the bar chart with labels
ggplot(top_sibling_groups, aes(x = reorder(Parent_Pair, Offspring_Count), y = Offspring_Count, fill = Offspring_Count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Offspring_Count), hjust = -0.2, size = 5) +  # Add text labels
  coord_flip() +  # Flip for better readability
  labs(title = "Top 20 Parent Pairs with Most Offspring in Sparouine",
       x = "Parent Pair",
       y = "Number of Offspring") +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "aquamarine3") +
  theme(axis.text.y = element_text(size = 10),  # Adjust text size for readability
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


# Merge geographic coordinates for fathers and mothers
parent_geo <- parent_data%>%
  left_join(geo_data, by = c("FatherID" = "ID")) %>%
  rename(Father_lat = lat, Father_long = long) %>%
  left_join(geo_data, by = c("MotherID" = "ID")) %>%
  rename(Mother_lat = lat, Mother_long = long)

# Calculate distances between fathers and mothers in meters
parent_geo <- parent_geo %>%
  mutate(Distance_m = ifelse(!is.na(Father_lat) & !is.na(Mother_lat),
                             distHaversine(cbind(Father_long, Father_lat),
                                           cbind(Mother_long, Mother_lat)),
                             NA))

# Summarize offspring count and average distance for each parent pair
sibling_groups <- parent_geo %>%
  group_by(FatherID, MotherID) %>%
  summarise(Offspring_Count = n(),
            Avg_Distance_m = mean(Distance_m, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(Offspring_Count))

# Select top 20 parent pairs for visualization
top_sibling_groups <- sibling_groups %>%
  slice_max(order_by = Offspring_Count, n = 15)

# Plot the bar chart with distances displayed as labels only if not NA
#1100X700
ggplot(top_sibling_groups, aes(x = reorder(paste(FatherID, MotherID, sep = " & "), Offspring_Count),
                               y = Offspring_Count, fill = Offspring_Count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(!is.na(Avg_Distance_m), 
                               paste0("Dist: ", round(Avg_Distance_m, 2), " m"), 
                               "")),
            hjust = -0.2, size = 4) +
  coord_flip() +
  labs(title = "Top 20 Parent Pairs with Most Offspring and Distances in Regina",
       x = "Parent Pair (Father & Mother)",
       y = "Number of Offspring") +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "mediumorchid4") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))



############# Half sib-parents ############
# Merge geographic coordinates for offspring and both parents
parent_geo <- parent_data %>%
  left_join(geo_data, by = c("OffspringID" = "ID")) %>%
  rename(Offspring_lat = lat, Offspring_long = long) %>%
  left_join(geo_data, by = c("FatherID" = "ID")) %>%
  rename(Parent1_lat = lat, Parent1_long = long) %>%
  left_join(geo_data, by = c("MotherID" = "ID")) %>%
  rename(Parent2_lat = lat, Parent2_long = long)

# Calculate distances between offspring and parents
parent_geo <- parent_geo %>%
  mutate(
    Distance_To_Parent1_m = ifelse(!is.na(Parent1_lat) & !is.na(Offspring_lat),
                                   distHaversine(cbind(Parent1_long, Parent1_lat),
                                                 cbind(Offspring_long, Offspring_lat)),
                                   NA),
    Distance_To_Parent2_m = ifelse(!is.na(Parent2_lat) & !is.na(Offspring_lat),
                                   distHaversine(cbind(Parent2_long, Parent2_lat),
                                                 cbind(Offspring_long, Offspring_lat)),
                                   NA)
  )

# Combine Parent1 and Parent2 into a single generic grouping for analysis
parent_geo_long <- parent_geo %>%
  select(OffspringID, Parent1 = FatherID, Parent2 = MotherID, Distance_To_Parent1_m, Distance_To_Parent2_m) %>%
  pivot_longer(cols = starts_with("Distance_To_Parent"), names_to = "Parent", values_to = "Distance_m") %>%
  mutate(ParentID = ifelse(Parent == "Distance_To_Parent1_m", Parent1, Parent2)) %>%
  filter(!is.na(ParentID))

# Summarize offspring counts and distances for each parent
half_sib_parents <- parent_geo_long %>%
  group_by(ParentID) %>%
  summarise(
    Offspring_Count = n(),
    Avg_Distance_m = mean(Distance_m, na.rm = TRUE),
    Offspring_IDs = paste(OffspringID, collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(Offspring_Count > 1) %>%
  arrange(desc(Offspring_Count))

half_sib_parents

# Save the combined results
output_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.4-parentage_with_colony/R_results_on_plots/"
output_file <- file.path(output_path, "Half_Siblings_Parents_with_Distances_REG.csv")
write.csv(half_sib_parents, output_file, row.names = FALSE)



##### Visualization: Offspring Distribution and Parent Distances #####
# Filter the top 20 parents by Offspring_Count
top_20_half_sib_parents <- half_sib_parents %>%
  slice_max(order_by = Offspring_Count, n = 20)

# Bar chart with offspring count inside the bars and average distance outside (only if not NA)
#1000x700
ggplot(top_20_half_sib_parents, aes(x = reorder(ParentID, Offspring_Count), y = Offspring_Count, fill = Offspring_Count)) +
  geom_bar(stat = "identity") +
  # Add offspring count inside the bars
  geom_text(aes(label = Offspring_Count), 
            color = "white", 
            hjust = 1.2, 
            size = 5) +
  # Add average distance outside the bars (only if not NA)
  geom_text(aes(label = ifelse(!is.na(Avg_Distance_m), 
                               paste0(round(Avg_Distance_m, 2), " m"), 
                               "")), 
            hjust = -0.2, 
            size = 4) +
  coord_flip() +
  labs(title = "Top 20 parents by number of offspring (half-sib) and average distance in Regina",
       x = "Parent ID",
       y = "Number of Offspring") +
  theme_minimal() +
  scale_fill_gradient(low = "lightblue", high = "mediumorchid4") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

