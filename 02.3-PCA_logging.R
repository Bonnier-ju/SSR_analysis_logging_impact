################################################################################
##################### Principal component analysis #############################
################################################################################

library(adegenet)
library(ggplot2)
library(dplyr)
library(readr)

# File paths
genepop_path <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Data/HKO50_PAI74/PAI74_HKO50_logged.gen"
treeinfo_path <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Data/HKO50_PAI74/PAI74_HKO50_tree_info.csv"
output_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/02-populations_structure/02.3-PCA"


genind_obj <- read.genepop(genepop_path, ncode = 3)
tree_info <- read_csv(treeinfo_path)

# Assign population names in genind object based on metadata
ind_names <- indNames(genind_obj)
meta_filtered <- tree_info %>% filter(ID %in% ind_names)

# Check that order matches
meta_filtered <- meta_filtered[match(ind_names, meta_filtered$ID), ]

# Set pop (required for adegenet methods)
pop(genind_obj) <- meta_filtered$plot

# Add DBH class
meta_filtered <- meta_filtered %>%
  mutate(
    cohort = case_when(
      dbh < 5 ~ "SED",
      dbh < 30 ~ "INT",
      dbh >= 30 ~ "ADL"
    )
  )

# Perform PCA
X <- tab(genind_obj, NA.method = "mean")
pca <- dudi.pca(X, scannf = FALSE, nf = 3)
summary(pca)

# Extract PCA scores and merge metadata
pca_df <- as.data.frame(pca$li)
pca_df$ID <- ind_names
pca_df <- left_join(pca_df, meta_filtered, by = "ID")

# Define cohort color and plot shape
cohort_colors <- c(SED = "blue", INT = "green", ADL = "red")
plot_shapes <- c(HKO50 = 16, PAI74 = 17)  # 16 = circle, 17 = triangle

# ------------------ PCA colored by plot and cohorts ------------------
plot1 <- ggplot(pca_df, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = cohort, shape = plot), size = 3, alpha = 0.8) +
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  scale_shape_manual(values = plot_shapes, name = "Plot") +
  labs(title = "PCA of SSR data - HKO50 vs PAI74",
       x = paste0("PC1 (", round(pca$eig[1] / sum(pca$eig) * 100, 1), "%)"),
       y = paste0("PC2 (", round(pca$eig[2] / sum(pca$eig) * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")

plot1

ggsave(file.path(output_dir, "PCA_by_plot_and_cohorts.png"), plot1, width = 8, height = 6, dpi = 300, bg = "white")


# ------------------ PCA colored by plot ------------------
plot2 <- ggplot(pca_df, aes(x = Axis1, y = Axis2, color = plot, shape = plot)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("HKO50" = "darkblue", "PAI74" = "orange")) +
  scale_shape_manual(values = c("HKO50" = 16, "PAI74" = 17)) +
  labs(title = "PCA - colored by plot",
       x = paste0("Axis 1 (", round(pca$eig[1] / sum(pca$eig) * 100, 1), "%)"),
       y = paste0("Axis 2 (", round(pca$eig[2] / sum(pca$eig) * 100, 1), "%)")) +
  theme_minimal()

plot2

ggsave(file.path(output_dir, "PCA_by_plot.png"), plot2, width = 8, height = 6, dpi = 300, bg = "white")


