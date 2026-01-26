# ============================================================
# Script: 02.1-STRUCTURE_results.R
# Purpose: Visualize STRUCTURE Harvester outputs (Evanno plots) and generate STRUCTURE barplots (K clusters)
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: evanno.txt ; K*.indfile ; STRUCTURE .str file ; metadata (tree_info.csv)
# Outputs: figures/02_structure/*.png
# Figures/Tables: Evanno ΔK plot, Ln'(K) plot, STRUCTURE barplot (K selected)
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
# Root of the project (use relative paths for GitHub reproducibility)
base_path <- "."

# Evanno / Harvester paths
evanno_file <- file.path(
  base_path,
  "Analysis", "02-populations_structure", "02.1-population_structure", "02.1.1-STRUCTURE",
  "results_harvester_PAI74_HKO50", "evanno.txt"
)

harvester_dir <- file.path(
  base_path,
  "Analysis", "02-populations_structure", "02.1-population_structure", "02.1.1-STRUCTURE",
  "results_harvester_PAI74_HKO50"
)

# STRUCTURE barplot parameters
k_value <- 8

indfile_path <- file.path(harvester_dir, paste0("K", k_value, ".indfile"))

str_path <- file.path(
  base_path,
  "Data", "HKO50_PAI74", "PAI74_HKO50.str.txt"
)

meta_path <- file.path(
  base_path,
  "Data", "HKO50_PAI74", "PAI74_HKO50_tree_info.csv"
)

# Output folders
fig_dir_evanno <- file.path("figures", "02_structure", "evanno")
fig_dir_barplot <- file.path("figures", "02_structure", "barplots")
dir.create(fig_dir_evanno, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir_barplot, showWarnings = FALSE, recursive = TRUE)
## ======================================================

## ===================== LIBRARIES ======================
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
## ======================================================

## ===================== CHECK INPUTS ===================
stopifnot(file.exists(evanno_file))
stopifnot(file.exists(indfile_path))
stopifnot(file.exists(str_path))
stopifnot(file.exists(meta_path))
## ======================================================

################################################################################
############################ Evanno results visualization #######################
################################################################################

# Read and clean data
col_names <- c(
  "K", "Reps", "Mean_LnP_K", "Stdev_LnP_K",
  "Ln_prime_K", "Abs_Ln_doubleprime_K", "Delta_K"
)

evanno_data <- read_tsv(evanno_file, comment = "#", col_names = col_names, show_col_types = FALSE)

# Plot Delta K
delta_plot <- ggplot(evanno_data, aes(x = K, y = Delta_K)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Evanno Method - ΔK vs K", x = "K", y = "ΔK") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(
  filename = file.path(fig_dir_evanno, "evanno_plot_delta_K.png"),
  plot = delta_plot, width = 7, height = 5, dpi = 300, bg = "white"
)

# Plot Ln′(K)
lnprime_plot <- ggplot(evanno_data, aes(x = K, y = Ln_prime_K)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Evanno Method - Ln′(K) vs K", x = "K", y = "Ln′(K)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(
  filename = file.path(fig_dir_evanno, "lnprime_plot.png"),
  plot = lnprime_plot, width = 7, height = 5, dpi = 300, bg = "white"
)

################################################################################
########################### STRUCTURE barplot visualization ####################
################################################################################

# Read STR file to extract individual IDs
# (STRUCTURE format: two lines per individual => take odd lines)
str_data <- read_table(str_path, skip = 1, col_names = FALSE, show_col_types = FALSE)
individual_ids <- str_data$X1[seq(1, nrow(str_data), by = 2)]
n_ind <- length(individual_ids)

# Read Harvester .indfile (first replicate only, aligned to individuals)
ind_lines <- read_lines(indfile_path)
stopifnot(length(ind_lines) >= n_ind)

indfile_firstrep <- ind_lines[1:n_ind]

ind_probs <- indfile_firstrep %>%
  str_trim() %>%
  str_remove("^\\d+\\s+\\d+\\s+\\(\\d+\\)\\s+\\d+\\s+:\\s+") %>%
  str_split_fixed(" ", k_value) %>%
  as_tibble() %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ . / rowSums(across(everything())))) %>%
  mutate(ID = individual_ids)

# Rename cluster columns
colnames(ind_probs)[1:k_value] <- paste0("Cluster", 1:k_value)

# Read metadata and assign cohort based on DBH
meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(cohort = case_when(
    dbh < 5 ~ "SED",
    dbh < 30 ~ "INT",
    dbh >= 30 ~ "ADL",
    TRUE ~ NA_character_
  ))

# Add dominant cluster info
ind_probs <- ind_probs %>%
  rowwise() %>%
  mutate(
    max_prob = max(c_across(starts_with("Cluster")), na.rm = TRUE),
    dominant_cluster = names(pick(starts_with("Cluster")))[which.max(c_across(starts_with("Cluster")))]
  ) %>%
  ungroup()

# Merge and sort individuals by plot, cohort, dominant cluster, and descending max prob
plot_data <- ind_probs %>%
  left_join(meta, by = "ID") %>%
  arrange(plot, cohort, dominant_cluster, desc(max_prob)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Probability") %>%
  mutate(Cluster = factor(Cluster, levels = paste0("Cluster", 1:k_value)))

# Plot
plot_structure <- ggplot(plot_data, aes(x = ID, y = Probability, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ plot + cohort, scales = "free_x", space = "free") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = paste("STRUCTURE barplot - K =", k_value),
    x = "Individuals", y = "Assignment Probability"
  )

ggsave(
  filename = file.path(fig_dir_barplot, paste0("structure_barplot_PAI74_HKO50_K", k_value, ".png")),
  plot = plot_structure, width = 19, height = 5, dpi = 300, bg = "white"
)
