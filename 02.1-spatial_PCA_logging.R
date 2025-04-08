###########################################################
#          Spatial PCA on Regina logging plots           #
###########################################################

library(adegenet)
library(spdep)
library(interp)
library(ggplot2)
library(adespatial)

########### Pre-processing of data #############

# Define file paths
path_HKO <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/Full_data_HKO50_nSSR.csv"
path_PAI <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Data_initial/Regina/Regina_full_data.csv"

out_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/02-populations_structure/02.2-spatial_PCA/Graphs"

# Load data
HKO <- read.csv(path_HKO, header = TRUE, sep = ";")
PAI <- read.csv(path_PAI, header = TRUE, sep = ";")
# Add plot labels
HKO$plot <- "HKO50"
PAI$plot <- "PAI74"
# Reorder columns to put 'plot' in 6th position
HKO <- HKO[, c(1:5, ncol(HKO), 6:(ncol(HKO) - 1))]
PAI <- PAI[, c(1:5, ncol(PAI), 6:(ncol(PAI) - 1))]

# Combine data
data_all <- rbind(HKO, PAI)

# Extract genotypes (all columns after plot)
geno_cols <- (6 + 1):ncol(data_all)  # after plot
geno_HKO <- HKO[, geno_cols]
geno_PAI <- PAI[, geno_cols]
geno_all <- data_all[, geno_cols]

# Convert to genind objects
genind_HKO <- df2genind(geno_HKO, ploidy = 2, sep = ",", NA.char = "NA")
genind_PAI <- df2genind(geno_PAI, ploidy = 2, sep = ",", NA.char = "NA")
genind_all <- df2genind(geno_all, ploidy = 2, sep = ",", NA.char = "NA")

# Coordinates
coords_HKO <- HKO[, c("lat", "long")]
coords_PAI <- PAI[, c("lat", "long")]
coords_all <- data_all[, c("lat", "long")]

# sPCA on each dataset
# type = 6, number of neighour = 30
spca_HKO <- spca(genind_HKO, xy = coords_HKO, ask = FALSE, type = 6, scannf = FALSE)
spca_PAI <- spca(genind_PAI, xy = coords_PAI, ask = FALSE, type = 6, scannf = FALSE)
spca_all <- spca(genind_all, xy = coords_all, ask = FALSE, type = 6, scannf = FALSE)


################### Analysis #################

# Choose which plot to analyze
# Options: "HKO50", "PAI74", "ALL"
plot_name <- "HKO50"  # Change this line to switch dataset

if (plot_name == "HKO50") {
  spca_plot <- spca_HKO
  coords <- coords_HKO
  genind_data <- genind_HKO
  original_data <- HKO
  plot_color <- "#CDAD00" 
} else if (plot_name == "PAI74") {
  spca_plot <- spca_PAI
  coords <- coords_PAI
  genind_data <- genind_PAI
  original_data <- PAI
  plot_color <- "mediumorchid4"  
} else if (plot_name == "ALL") {
  spca_plot <- spca_all
  coords <- coords_all
  genind_data <- genind_all
  original_data <- data_all
  plot_color <- "dodgerblue3"  
}


# Eigenvalues 
png(file.path(out_dir, paste0("eigenvalues_", plot_name, ".png")), width = 800, height = 600)
barplot(spca_plot$eig,
        main = paste("sPCA Eigenvalues -", plot_name),
        col = spectral(length(spca_plot$eig)))
legend("topright", fill = spectral(2), legend = c("Global", "Local"))
abline(h = 0, col = "red")
dev.off()

### Screeplot
png(file.path(out_dir, paste0("screeplot_", plot_name, ".png")), width = 800, height = 600)
screeplot(spca_plot, main = paste("Moran's I vs Variance -", plot_name))
dev.off()

### Global & Local structure tests
X <- scaleGen(genind_data, NA.method = "mean")

g_test <- global.rtest(X, spca_plot$lw, nperm = 9999)
l_test <- local.rtest(X, spca_plot$lw, nperm = 9999)

# Optional: print results to console
print(g_test)
print(l_test)

# Optional: export test results to text
capture.output(list(Global_Test = g_test, Local_Test = l_test),
               file = file.path(out_dir, paste0("structure_tests_", plot_name, ".txt")))

### sPCA Scores Plot
score_df <- data.frame(
  sPCA1 = spca_plot$ls[, 1],
  sPCA2 = spca_plot$ls[, 2],
  lat = coords$lat,
  long = coords$long
)

score_plot <- ggplot(score_df, aes(x = sPCA1, y = sPCA2)) +
  geom_point(size = 3, color = "darkblue") +
  labs(title = paste("sPCA Scores -", plot_name),
       x = "sPCA Axis 1", y = "sPCA Axis 2") +
  theme_minimal()

ggsave(file.path(out_dir, paste0("scores_", plot_name, ".png")),
       plot = score_plot, width = 7, height = 6)

### Interpolated lagged scores (Axis 1)
png(file.path(out_dir, paste0("interpolation_", plot_name, ".png")), width = 800, height = 700)
interp_data <- interp(coords$long, coords$lat, spca_plot$ls[, 1], duplicate = "mean")
image(interp_data, col = azur(100))
points(coords$long, coords$lat)
dev.off()

### Allele contributions (Axis 1)
allele_contrib <- spca_plot$c1[, 1]^2
names(allele_contrib) <- rownames(spca_plot$c1)

png(file.path(out_dir, paste0("loadingplot_", plot_name, ".png")), width = 1000, height = 600)
loadingplot(allele_contrib, xlab = "Alleles", ylab = "Weight",
            main = paste("Allele Contributions -", plot_name))
dev.off()

### Marker contributions (Boxplot by locus)
png(file.path(out_dir, paste0("marker_boxplot_", plot_name, ".png")), width = 1800, height = 800)
boxplot(allele_contrib ~ genind_data$loc.fac, las = 3, col = "grey",
        ylab = "Contribution", xlab = "Locus",
        main = paste("Marker Contributions -", plot_name))
dev.off()