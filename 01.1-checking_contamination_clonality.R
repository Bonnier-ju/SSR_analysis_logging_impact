# ============================================================
# Script: 1.1-checking_contamination_clonality.R
# Purpose: Identify potential clonality/contamination using RClone (MLG/MLL, psex, distance distributions).
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: data/Full_data_<PLOT>.csv (SSRSeq loci)
# Outputs: results/01_qc_rclone/suspected_contaminations_<PLOT>.csv
#          results/01_qc_rclone/genetic_distance_hist_<PLOT>.png (optional)
# Figures/Tables: QC histogram + MLL list
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
plot_name <- "HKO50"  # change when needed

# File paths (no hard-coded absolute paths)
input_file <- file.path("data", paste0("Full_data_", plot_name, ".csv"))
output_dir <- file.path("results", "01_qc_rclone")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_csv <- file.path(output_dir, paste0("suspected_contaminations_", plot_name, ".csv"))
output_png <- file.path(output_dir, paste0("genetic_distance_hist_", plot_name, ".png"))

# MLL thresholds
alpha2_mll_main   <- 14  # for clonality/diversity indices (as used originally)
alpha2_mll_export <- 5   # for exporting suspect MLL groups (stricter)

# Reproducibility (for resampling/simulations)
set.seed(123)
## ======================================================

## ===================== LIBRARIES ======================
# NOTE: Install packages outside this script (README / setup script)
library(RClone)
library(dplyr)
library(pegas)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(input_file))
data <- read.csv(input_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)

indiv_ids <- data$ID
markers <- data %>% dplyr::select(starts_with("SSRSeq"))
markers[is.na(markers)] <- "000"
## ======================================================

## ===================== MLG CHECKS =====================
# Identify Unique Multilocus Genotypes (MLGs)
list_all_tab(markers)

MLG_tab(markers, vecpop = NULL)

MLGlist <- MLG_list(markers)
for (i in seq_along(MLGlist)) {
  cat("MLG", i, "contains individuals (IDs):", indiv_ids[MLGlist[[i]]], "\n")
}

# Allelic frequencies per locus
freq_RR(markers)
## ======================================================

## ===================== LOCI/UNITS RESAMPLING ==========
# Evaluate loci informativeness (genotype accumulation curves + diversity indices)
res <- sample_loci(markers, nbrepeat = 100, He = TRUE, graph = TRUE, bar = TRUE)

res$res_MLG
res$res_alleles

sample_units(markers, nbrepeat = 100)
## ======================================================

## ===================== MLL INFERENCE (pgen / psex) =====
# Probability of observing a given multilocus genotype by chance
pgen(markers)

# Probability that a genotype appears more than once purely by sexual reproduction
res_psex <- psex(markers, RR = TRUE, nbrepeat = 100)
res_psex
## ======================================================

## ===================== GENETIC DISTANCES ===============
# Observed pairwise genetic distances (allele differences)
respop <- genet_dist(markers)

# Theoretical distribution of genetic distances based on simulations
ressim <- genet_dist_sim(markers, nbrepeat = 100)

# Simulated distribution excluding self-fertilization
ressimWS <- genet_dist_sim(markers, genet = TRUE, nbrepeat = 100)

# Common x limit
limx <- max(
  max(respop$distance_matrix),
  max(ressim$distance_matrix),
  max(ressimWS$distance_matrix)
)

# Save histogram (optional but recommended for reproducibility)
png(filename = output_png, width = 700, height = 600)

hist(
  respop$distance_matrix,
  freq = FALSE,
  col = "mediumorchid4",
  main = paste0("Genetic distance distributions - ", plot_name),
  xlab = "Genetic distances",
  breaks = seq(0, limx + 1, 1),
  xlim = c(0, 100),
  ylim = c(0, 0.20)
)

hist(
  ressim$distance_matrix,
  freq = FALSE,
  add = TRUE,
  col = rgb(0.7, 0.9, 1, 0.5),
  breaks = seq(0, limx + 1, 1)
)

hist(
  ressimWS$distance_matrix,
  freq = FALSE,
  add = TRUE,
  col = rgb(0.9, 0.5, 1, 0.3),
  breaks = seq(0, limx + 1, 1)
)

legend(
  x = limx * 0.1,
  y = 0.20,
  legend = c("Observed data", "Simulated data", "Simulated (no selfing)"),
  fill = c("mediumorchid4", rgb(0.7, 0.9, 1, 0.5), rgb(0.9, 0.5, 1, 0.3)),
  bg = "white",
  box.lwd = 1,
  cex = 0.8
)

dev.off()

# Frequency table: number of pairs sharing a given genetic distance
table(respop$distance_matrix)
## ======================================================

## ===================== MLL GENERATION ==================
# Main MLLs (used for diversity/clonality indices)
MLLlist_main <- MLL_generator(markers, alpha2 = alpha2_mll_main)
## ======================================================

## ===================== CLONALITY INDICES ===============
# Genotypic diversity and clonality indices (MLG-based)
clonal_index(markers)

# Same indices using MLLs
clonal_index(markers, listMLL = MLLlist_main)
## ======================================================

## ===================== PARETO INDEX ====================
Pareto_index(markers)
Pareto_index(markers, listMLL = MLLlist_main)
Pareto_index(markers, full = TRUE, graph = TRUE, legends = 2)
## ======================================================

## ===================== EXPORT SUSPECT MLLs ==============
# Stricter MLL grouping for exporting suspicious clusters
MLLlist_export <- MLL_generator(markers, alpha2 = alpha2_mll_export)

# Assign each individual to an MLL ID
MLL_vector <- rep(NA_integer_, length(indiv_ids))

for (i in seq_along(MLLlist_export)) {
  individuals_in_MLL <- MLLlist_export[[i]]
  MLL_vector[individuals_in_MLL] <- i
}

MLL_df <- data.frame(ID = indiv_ids, MLL = MLL_vector)

# Identify MLL groups with multiple individuals
suspected_MLLs <- MLL_df %>%
  group_by(MLL) %>%
  filter(!is.na(MLL) & n() > 1) %>%
  arrange(MLL)

print(suspected_MLLs)

# Export
write.csv(suspected_MLLs, output_csv, row.names = FALSE)

cat("\nSaved suspected MLL groups to:", output_csv, "\n")
cat("Saved histogram to:", output_png, "\n")
## ======================================================
