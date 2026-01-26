# ============================================================
# Script: 03.1-private_allelic_richness.R
# Purpose: Compute private allelic richness proxy (count of private alleles) for a given plot and DBH cohort.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: data/Full_plot_infos.csv (must contain ID, plot, dbh and allele columns ending with _1/_2)
# Outputs: console output (optionally results/03_privateAR/privateAR_<plot>_<dbhmin>-<dbhmax>.csv)
# Figures/Tables: Private allelic richness results (used in diversity section / table)
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
base_path <- "."

input_file <- file.path(base_path, "data", "Full_plot_infos.csv")

plot_name <- "HKO50"       # Plot to analyze
dbh_range <- c(0, 2)       # DBH range to analyze: c(min, max), max excluded

# Optional output file (set to TRUE to export)
export_results <- FALSE
output_dir <- file.path("results", "03_privateAR")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_csv <- file.path(output_dir, paste0("privateAR_", plot_name, "_", dbh_range[1], "-", dbh_range[2], ".csv"))
## ======================================================

## ===================== LIBRARIES ======================
library(dplyr)
library(adegenet)
library(poppr)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(input_file))
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Basic checks
required_cols <- c("ID", "plot", "dbh")
stopifnot(all(required_cols %in% names(data)))
## ======================================================

## ===================== FILTERING ======================
filtered_data <- data %>%
  filter(plot == plot_name, dbh >= dbh_range[1], dbh < dbh_range[2])

if (nrow(filtered_data) == 0) {
  stop("No individuals found after filtering. Check plot_name and dbh_range.")
}

allele_cols <- grep("_[12]$", names(filtered_data), value = TRUE)
if (length(allele_cols) == 0) {
  stop("No allele columns detected. Expected columns ending with '_1' and '_2'.")
}
## ======================================================

## ===================== GENIND =========================
genind_obj <- df2genind(
  filtered_data[, allele_cols],
  ploidy = 2,
  sep = ",",
  ind.names = filtered_data$ID,
  pop = as.factor(rep(paste0(plot_name, "_", paste(dbh_range, collapse = "-")), nrow(filtered_data)))
)
## ======================================================

## ===================== PRIVATE ALLELES ================
# private_alleles() returns a table of allele counts per population
pa <- private_alleles(
  genind_obj,
  report = "table",
  level = "population",
  count.alleles = TRUE,
  drop = FALSE
)

# Here: total number of alleles private to this population (presence/absence across loci)
total_pa <- sum(pa > 0, na.rm = TRUE)
## ======================================================

## ===================== OUTPUT =========================
cat("Plot:", plot_name, "\n")
cat("DBH range:", paste(dbh_range, collapse = "–"), "\n")
cat("Private alleles (count):", total_pa, "\n")

if (export_results) {
  out_df <- data.frame(
    plot = plot_name,
    dbh_min = dbh_range[1],
    dbh_max = dbh_range[2],
    private_alleles_count = total_pa
  )
  write.csv(out_df, output_csv, row.names = FALSE)
  cat("Saved:", output_csv, "\n")
}
## ======================================================
