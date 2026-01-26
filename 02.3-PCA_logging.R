# ============================================================
# Script: 02.3-PCA_logging.R
# Purpose: PCA on SSR genotypes (HKO50 vs PAI74) using adegenet + ade4, with plots by plot and cohort.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: data/HKO50_PAI74/PAI74_HKO50_logged.gen ; data/HKO50_PAI74/PAI74_HKO50_tree_info.csv
# Outputs: figures/02_pca/PCA_by_plot_and_cohorts.png ; figures/02_pca/PCA_by_plot.png
# Figures/Tables: PCA figures (plot/cohort)
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
base_path <- "."

genepop_path <- file.path(base_path, "data", "HKO50_PAI74", "PAI74_HKO50_logged.gen")
treeinfo_path <- file.path(base_path, "data", "HKO50_PAI74", "PAI74_HKO50_tree_info.csv")

fig_dir <- file.path("figures", "02_pca")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Plot aesthetics
cohort_colors <- c(SED = "blue", INT = "green", ADL = "red")
plot_shapes <- c(HKO50 = 16, PAI74 = 17) # 16 = circle, 17 = triangle
plot_colors <- c(HKO50 = "darkblue", PAI74 = "orange")
## ======================================================

## ===================== LIBRARIES ======================
library(adegenet)
library(ade4)
library(ggplot2)
library(dplyr)
library(readr)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(genepop_path))
stopifnot(file.exists(treeinfo_path))

genind_obj <- read.genepop(genepop_path, ncode = 3)
tree_info <- read_csv(treeinfo_path, show_col_types = FALSE)
## ======================================================

## ===================== METADATA MATCH =================
ind_names <- indNames(genind_obj)
meta_filtered <- tree_info %>% filter(ID %in% ind_names)

# Ensure same order as genind
meta_filtered <- meta_filtered[match(ind_names, meta_filtered$ID), ]

# Safety check: all individuals matched
stopifnot(!any(is.na(meta_filtered$ID)))

# Set population factor (required for some adegenet methods)
pop(genind_obj) <- meta_filtered$plot

# Add DBH cohort class
meta_filtered <- meta_filtered %>%
  mutate(
    cohort = case_when(
      dbh < 5 ~ "SED",
      dbh < 30 ~ "INT",
      dbh >= 30 ~ "ADL",
      TRUE ~ NA_character_
    )
  )
## ======================================================

## ===================== PCA ============================
X <- tab(genind_obj, NA.method = "mean")
pca <- dudi.pca(X, scannf = FALSE, nf = 3)
print(summary(pca))

# Percent variance explained
pc1_var <- round(pca$eig[1] / sum(pca$eig) * 100, 1)
pc2_var <- round(pca$eig[2] / sum(pca$eig) * 100, 1)

# PCA scores + metadata
pca_df <- as.data.frame(pca$li)
pca_df$ID <- ind_names
pca_df <- left_join(pca_df, meta_filtered, by = "ID")
## ======================================================

## ===================== PLOTS ==========================
# PCA colored by plot and cohort
plot1 <- ggplot(pca_df, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = cohort, shape = plot), size = 3, alpha = 0.8) +
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  scale_shape_manual(values = plot_shapes, name = "Plot") +
  labs(
    title = "PCA of SSR data - HKO50 vs PAI74",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(
  filename = file.path(fig_dir, "PCA_by_plot_and_cohorts.png"),
  plot = plot1, width = 8, height = 6, dpi = 300, bg = "white"
)

# PCA colored by plot only
plot2 <- ggplot(pca_df, aes(x = Axis1, y = Axis2, color = plot, shape = plot)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = plot_colors) +
  scale_shape_manual(values = plot_shapes) +
  labs(
    title = "PCA - colored by plot",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, "PCA_by_plot.png"),
  plot = plot2, width = 8, height = 6, dpi = 300, bg = "white"
)
## ======================================================
