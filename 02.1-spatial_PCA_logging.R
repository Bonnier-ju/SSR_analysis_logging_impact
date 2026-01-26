# ============================================================
# Script: 02.1-spatial_PCA_logging.R
# Purpose: Run sPCA (adegenet) for HKO50, PAI74, and combined dataset, then export key diagnostic plots
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: data/Full_data_HKO50_nSSR.csv ; data/Regina_full_data.csv ; optional isolines shapefile
# Outputs: figures/02_spca/*.png ; results/02_spca_tests/*.txt (optional)
# Figures/Tables: sPCA eigenvalues, screeplot, score plot, lagged score map, loadings
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
plot_name <- "ALL"  # Options: "HKO50", "PAI74", "ALL"

# Input files (use relative paths for GitHub reproducibility)
path_HKO <- file.path("data", "Full_data_HKO50_nSSR.csv")
path_PAI <- file.path("data", "Regina_full_data.csv")

# Output directories
fig_dir <- file.path("figures", "02_spca")
res_dir <- file.path("results", "02_spca_tests")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# sPCA parameters
nb_neighbors <- 30
spca_type <- 6
nperm_tests <- 9999

# Plot colors
plot_colors <- list(
  HKO50 = "#CDAD00",
  PAI74 = "mediumorchid4",
  ALL  = "dodgerblue3"
)

# Optional: isolines shapefile (set paths per dataset if you have them)
# If the file is missing, the map will be produced without isolines.
isolines_shp <- list(
  HKO50 = NA_character_,
  PAI74 = NA_character_,
  ALL  = NA_character_  # you can set e.g. Regina isolines here if you want
)

# Bounding box expansion (degrees, for WGS84)
bbox_buffer <- 0.001

# Reproducibility
set.seed(123)
## ======================================================

## ===================== LIBRARIES ======================
library(adegenet)
library(adespatial)
library(ggplot2)

# For optional isolines + interpolation map
library(sf)
library(akima)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(path_HKO))
stopifnot(file.exists(path_PAI))

HKO <- read.csv(path_HKO, header = TRUE, sep = ";", stringsAsFactors = FALSE)
PAI <- read.csv(path_PAI, header = TRUE, sep = ";", stringsAsFactors = FALSE)

# Add plot labels
HKO$plot <- "HKO50"
PAI$plot <- "PAI74"

# Reorder columns to put 'plot' in 6th position
HKO <- HKO[, c(1:5, ncol(HKO), 6:(ncol(HKO) - 1))]
PAI <- PAI[, c(1:5, ncol(PAI), 6:(ncol(PAI) - 1))]

# Combine data
data_all <- rbind(HKO, PAI)

# Extract genotypes (all columns after 'plot')
geno_cols <- (6 + 1):ncol(data_all)
geno_HKO <- HKO[, geno_cols]
geno_PAI <- PAI[, geno_cols]
geno_all <- data_all[, geno_cols]

# Convert to genind objects
genind_HKO <- df2genind(geno_HKO, ploidy = 2, sep = ",", NA.char = "NA")
genind_PAI <- df2genind(geno_PAI, ploidy = 2, sep = ",", NA.char = "NA")
genind_all <- df2genind(geno_all, ploidy = 2, sep = ",", NA.char = "NA")

# Coordinates (ensure numeric)
coords_HKO <- HKO[, c("lat", "long")]
coords_PAI <- PAI[, c("lat", "long")]
coords_all <- data_all[, c("lat", "long")]

coords_HKO$lat <- as.numeric(coords_HKO$lat)
coords_HKO$long <- as.numeric(coords_HKO$long)
coords_PAI$lat <- as.numeric(coords_PAI$lat)
coords_PAI$long <- as.numeric(coords_PAI$long)
coords_all$lat <- as.numeric(coords_all$lat)
coords_all$long <- as.numeric(coords_all$long)
## ======================================================

## ===================== sPCA RUN =======================
# Run sPCA on each dataset
spca_HKO <- spca(genind_HKO, xy = coords_HKO, ask = FALSE, type = spca_type, nb = nb_neighbors, scannf = FALSE)
spca_PAI <- spca(genind_PAI, xy = coords_PAI, ask = FALSE, type = spca_type, nb = nb_neighbors, scannf = FALSE)
spca_all <- spca(genind_all, xy = coords_all, ask = FALSE, type = spca_type, nb = nb_neighbors, scannf = FALSE)
## ======================================================

## ===================== SELECT DATASET =================
if (plot_name == "HKO50") {
  spca_plot <- spca_HKO
  coords <- coords_HKO
  genind_data <- genind_HKO
  plot_color <- plot_colors$HKO50
} else if (plot_name == "PAI74") {
  spca_plot <- spca_PAI
  coords <- coords_PAI
  genind_data <- genind_PAI
  plot_color <- plot_colors$PAI74
} else if (plot_name == "ALL") {
  spca_plot <- spca_all
  coords <- coords_all
  genind_data <- genind_all
  plot_color <- plot_colors$ALL
} else {
  stop("plot_name must be one of: 'HKO50', 'PAI74', 'ALL'")
}
## ======================================================

## ===================== EIGENVALUES ====================
png(file.path(fig_dir, paste0("eigenvalues_", plot_name, ".png")), width = 800, height = 600)
barplot(
  spca_plot$eig,
  main = paste("sPCA Eigenvalues -", plot_name),
  col = adegenet::spectral(length(spca_plot$eig))
)
legend("topright", fill = adegenet::spectral(2), legend = c("Global", "Local"))
abline(h = 0, col = "red")
dev.off()

png(file.path(fig_dir, paste0("screeplot_", plot_name, ".png")), width = 800, height = 600)
screeplot(spca_plot, main = paste("Moran's I vs Variance -", plot_name))
dev.off()
## ======================================================

## ===================== GLOBAL/LOCAL TESTS =============
X <- scaleGen(genind_data, NA.method = "mean")

g_test <- global.rtest(X, spca_plot$lw, nperm = nperm_tests)
l_test <- local.rtest(X, spca_plot$lw, nperm = nperm_tests)

# Save tests to a text file for traceability
sink(file.path(res_dir, paste0("spca_tests_", plot_name, ".txt")))
cat("Global test:\n")
print(g_test)
cat("\nLocal test:\n")
print(l_test)
sink()
## ======================================================

## ===================== SCORE PLOT ======================
score_df <- data.frame(
  sPCA1 = spca_plot$ls[, 1],
  sPCA2 = spca_plot$ls[, 2],
  lat = coords$lat,
  long = coords$long
)

score_plot <- ggplot(score_df, aes(x = sPCA1, y = sPCA2)) +
  geom_point(size = 3, color = plot_color) +
  labs(
    title = paste("sPCA Scores -", plot_name),
    x = "sPCA Axis 1",
    y = "sPCA Axis 2"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, paste0("scores_", plot_name, ".png")),
  plot = score_plot, width = 7, height = 6, bg = "white"
)
## ======================================================

## =========== LAGGED SCORE MAP + OPTIONAL ISOLINES ======
# Prepare bounding box from points
coords_sf <- st_as_sf(data.frame(coords), coords = c("long", "lat"), crs = 4326)
bbox <- st_bbox(coords_sf)
bbox_exp <- bbox
bbox_exp["xmin"] <- bbox["xmin"] - bbox_buffer
bbox_exp["ymin"] <- bbox["ymin"] - bbox_buffer
bbox_exp["xmax"] <- bbox["xmax"] + bbox_buffer
bbox_exp["ymax"] <- bbox["ymax"] + bbox_buffer

# Interpolation of lagged score (1st axis)
interp_data <- akima::interp(coords$long, coords$lat, spca_plot$ls[, 1], duplicate = "mean")
interp_df <- expand.grid(x = interp_data$x, y = interp_data$y)
interp_df$z <- as.vector(interp_data$z)

# Read and crop isolines if available
isolines_crop <- NULL
this_isolines <- isolines_shp[[plot_name]]

if (!is.na(this_isolines) && file.exists(this_isolines)) {
  isolines <- st_read(this_isolines, quiet = TRUE)
  isolines_wgs84 <- st_transform(isolines, crs = 4326)
  isolines_crop <- st_crop(isolines_wgs84, bbox_exp)
}

lagged_map <- ggplot() +
  geom_raster(data = interp_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(
    colors = c(plot_colors$PAI74, "white", plot_colors$HKO50),
    name = "Lagged\nScore 1"
  ) +
  { if (!is.null(isolines_crop)) geom_sf(data = isolines_crop, color = "grey30", size = 0.3, alpha = 0.8) } +
  geom_point(data = coords, aes(x = long, y = lat), color = "black", size = 2) +
  coord_sf(
    xlim = c(bbox_exp["xmin"], bbox_exp["xmax"]),
    ylim = c(bbox_exp["ymin"], bbox_exp["ymax"]),
    expand = FALSE
  ) +
  labs(
    title = paste("Lagged Score 1 (Interpolated) -", plot_name),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, paste0("lagged_isolines_", plot_name, ".png")),
  plot = lagged_map, width = 10, height = 7, dpi = 300, bg = "white"
)
## ======================================================

## ===================== LOADINGS =======================
# Allele contributions (Axis 1)
allele_contrib <- spca_plot$c1[, 1]^2
names(allele_contrib) <- rownames(spca_plot$c1)

png(file.path(fig_dir, paste0("loadingplot_", plot_name, ".png")), width = 1000, height = 600)
loadingplot(
  allele_contrib,
  xlab = "Alleles",
  ylab = "Weight",
  main = paste("Allele Contributions -", plot_name)
)
dev.off()

# Marker contributions (boxplot by locus)
png(file.path(fig_dir, paste0("marker_boxplot_", plot_name, ".png")), width = 1800, height = 800)
boxplot(
  allele_contrib ~ genind_data$loc.fac,
  las = 3,
  col = "grey",
  ylab = "Contribution",
  xlab = "Locus",
  main = paste("Marker Contributions -", plot_name)
)
dev.off()
## ======================================================
