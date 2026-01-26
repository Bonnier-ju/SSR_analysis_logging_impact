# ============================================================
# Script: 03.3-bootstrapping_diversity_index_by_plots.R
# Purpose: Bootstrap comparisons of genetic indices (AR, He, Ho, Fi), selfing (S), and SGS slope (b-log / Sp) across plots and cohorts.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: Bootstrap CSV files from SPAGeDi (diversity, selfing, SGS)
# Outputs: results/03_bootstrap/*.csv ; figures/03_bootstrap/*.png
# Figures/Tables: Bootstrap comparison plots (inter-plot and intra-plot)
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
base_path <- "."

# ---- Diversity bootstrap inputs (per plot) ----
file_pai <- file.path(base_path, "results", "spagedi", "PAI74", "results_PAI74_full_bootstrap.csv")
file_hko <- file.path(base_path, "results", "spagedi", "HKO50", "results_HKO50_full_bootstrap.csv")

# ---- Selfing bootstrap input (both plots together) ----
file_selfing <- file.path(base_path, "results", "spagedi", "PAI74_HKO50", "results_HKO50_PAI74_bootstrap_selfing.csv")

# ---- SGS bootstrap input (b-log / Sp comparisons) ----
file_sgs <- file.path(base_path, "results", "spagedi", "PAI74_HKO50", "results_HKO50_PAI74_Sp_sign_btw_plot.csv")

# ---- Main output dirs ----
out_dir <- file.path("results", "03_bootstrap")
fig_dir <- file.path("figures", "03_bootstrap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Output files ----
out_diversity_csv <- file.path(out_dir, "bootstrap_comparisons_diversity_inter_intra_plots.csv")
out_selfing_csv <- file.path(out_dir, "final_selfing_comparisons.csv")
out_sgs_csv <- file.path(out_dir, "bootstrap_b_log_comparisons.csv")

# Bootstrap iterations
n_iter <- 10000
set.seed(123)

# Indicators (diversity)
indicators <- c("AR", "He", "Ho", "Fi")

# Biological cohorts to keep for intra-plot filtering
biological_cats <- c("SED", "INT", "ADL")

# Plot aesthetics
fill_colors <- c("TRUE" = "tomato", "FALSE" = "grey70")
## ======================================================

## ===================== LIBRARIES ======================
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)
## ======================================================

## ===================== CHECK INPUTS ===================
stopifnot(file.exists(file_pai))
stopifnot(file.exists(file_hko))
stopifnot(file.exists(file_selfing))
stopifnot(file.exists(file_sgs))
## ======================================================

################################################################################
########################## 1) BOOTSTRAP DIVERSITY ##############################
################################################################################

pai <- read_csv(file_pai, show_col_types = FALSE) %>% mutate(Plot = "PAI74")
hko <- read_csv(file_hko, show_col_types = FALSE) %>% mutate(Plot = "HKO50")
data_all <- bind_rows(pai, hko)

bootstrap_diff <- function(df1, df2, var, n_iter = 10000) {
  merged <- inner_join(
    df1 %>% select(Markers, value1 = !!sym(var)),
    df2 %>% select(Markers, value2 = !!sym(var)),
    by = "Markers"
  ) %>%
    mutate(Diff = value2 - value1) %>%
    filter(!is.na(Diff))
  
  if (nrow(merged) < 5) {
    return(data.frame(Obs_Diff = NA, CI_lower = NA, CI_upper = NA, Significant = NA))
  }
  
  diffs_boot <- replicate(n_iter, mean(sample(merged$Diff, replace = TRUE)))
  ci <- quantile(diffs_boot, c(0.025, 0.975), na.rm = TRUE)
  obs_diff <- mean(merged$Diff)
  
  data.frame(
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = !(0 >= ci[1] & 0 <= ci[2])
  )
}

plots <- unique(data_all$Plot)
categories <- unique(data_all$Category)
results_list <- list()

# Inter-plot comparisons (same category between plots)
for (cat in categories) {
  for (var in indicators) {
    g1 <- data_all %>% filter(Plot == "HKO50", Category == cat)
    g2 <- data_all %>% filter(Plot == "PAI74", Category == cat)
    
    res <- bootstrap_diff(g1, g2, var, n_iter = n_iter)
    res$Comparison <- "PAI74 vs HKO50"
    res$Category <- cat
    res$Variable <- var
    res$Type <- "Inter-plot"
    results_list[[length(results_list) + 1]] <- res
  }
}

# Intra-plot comparisons (pairwise within each plot)
for (plot in plots) {
  sub_data <- data_all %>% filter(Plot == plot)
  combs <- combn(unique(sub_data$Category), 2, simplify = FALSE)
  
  for (pair in combs) {
    for (var in indicators) {
      g1 <- sub_data %>% filter(Category == pair[1])
      g2 <- sub_data %>% filter(Category == pair[2])
      
      res <- bootstrap_diff(g1, g2, var, n_iter = n_iter)
      res$Comparison <- paste0(plot, ": ", pair[2], " vs ", pair[1])
      res$Category <- paste(pair[1], "vs", pair[2])
      res$Variable <- var
      res$Type <- "Intra-plot"
      results_list[[length(results_list) + 1]] <- res
    }
  }
}

final_results_div <- bind_rows(results_list) %>%
  select(Type, Comparison, Category, Variable, Obs_Diff, CI_lower, CI_upper, Significant)

valid_intra_biological <- function(category_string) {
  parts <- unlist(strsplit(category_string, " vs "))
  length(parts) == 2 && all(parts %in% biological_cats)
}

filtered_results_div <- final_results_div %>%
  filter(
    Type == "Inter-plot" |
      (Type == "Intra-plot" & sapply(Category, valid_intra_biological))
  ) %>%
  mutate(
    ComparisonGroup = case_when(
      Type == "Inter-plot" ~ "HKO50 vs PAI74",
      str_starts(Comparison, "HKO50") ~ "HKO50 only",
      str_starts(Comparison, "PAI74") ~ "PAI74 only",
      TRUE ~ "Other"
    )
  )

write.csv(filtered_results_div, out_diversity_csv, row.names = FALSE)
cat("Saved diversity bootstrap results:", out_diversity_csv, "\n")

# ---- Plotting (diversity) ----
plot_comparisons <- function(df, title_label, use_comparison_col = TRUE) {
  df$Variable <- factor(df$Variable, levels = c("AR", "He", "Ho", "Fi"))
  facet_var <- if (use_comparison_col) "Comparison_detail" else "Comparison"
  
  ggplot(df, aes(x = Variable, y = Obs_Diff, fill = as.character(Significant))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x", ncol = 2) +
    scale_fill_manual(values = fill_colors, name = "Significant") +
    labs(title = title_label, y = "Observed difference", x = "Genetic indicator") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

data_inter <- filtered_results_div %>%
  filter(ComparisonGroup == "HKO50 vs PAI74") %>%
  mutate(Comparison_detail = paste(Category, "-", Comparison))

data_pai <- filtered_results_div %>% filter(ComparisonGroup == "PAI74 only")
data_hko <- filtered_results_div %>% filter(ComparisonGroup == "HKO50 only")

plot1 <- plot_comparisons(data_inter, "PAI74 vs HKO50 comparisons", use_comparison_col = TRUE)
plot2 <- plot_comparisons(data_pai, "Intra-plot comparisons in PAI74", use_comparison_col = FALSE)
plot3 <- plot_comparisons(data_hko, "Intra-plot comparisons in HKO50", use_comparison_col = FALSE)

ggsave(file.path(fig_dir, "bootstrap_diversity_interplot.png"), plot1, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "bootstrap_diversity_intraplot_PAI74.png"), plot2, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "bootstrap_diversity_intraplot_HKO50.png"), plot3, width = 10, height = 7, dpi = 300, bg = "white")

################################################################################
########################## 2) BOOTSTRAP SELFING ###############################
################################################################################

selfing_data <- read_csv(file_selfing, show_col_types = FALSE)

# Marker columns: everything except metadata columns
marker_cols <- setdiff(names(selfing_data), c("Plot", "Population", "S(selfing rate)"))

selfing_data <- selfing_data %>%
  mutate(Group = paste(Plot, Population, sep = "_"))

groups <- unique(selfing_data$Group)
group_comparisons <- combn(groups, 2, simplify = FALSE)

results_selfing <- list()

for (comp in group_comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]
  
  m1 <- selfing_data %>% filter(Group == g1)
  m2 <- selfing_data %>% filter(Group == g2)
  
  # Skip if one group is missing or has multiple rows
  if (nrow(m1) != 1 || nrow(m2) != 1) next
  
  boot_diffs <- replicate(n_iter, {
    sampled_cols <- sample(marker_cols, length(marker_cols), replace = TRUE)
    mean(as.numeric(m2[1, sampled_cols]), na.rm = TRUE) -
      mean(as.numeric(m1[1, sampled_cols]), na.rm = TRUE)
  })
  
  obs_diff <- mean(as.numeric(m2[1, marker_cols]), na.rm = TRUE) -
    mean(as.numeric(m1[1, marker_cols]), na.rm = TRUE)
  
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  significant <- !(0 >= ci[1] & 0 <= ci[2])
  
  results_selfing[[length(results_selfing) + 1]] <- data.frame(
    Group1 = g1,
    Group2 = g2,
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = significant
  )
}

final_results_selfing <- bind_rows(results_selfing)
write.csv(final_results_selfing, out_selfing_csv, row.names = FALSE)
cat("Saved selfing bootstrap results:", out_selfing_csv, "\n")

################################################################################
########################## 3) BOOTSTRAP SGS (b-log) ############################
################################################################################

sgs_data <- read_csv(file_sgs, show_col_types = FALSE) %>%
  mutate(Group = paste(Plot, category, sep = "_"))

groups_sgs <- unique(sgs_data$Group)
group_comparisons_sgs <- combn(groups_sgs, 2, simplify = FALSE)

results_sgs <- list()

for (comp in group_comparisons_sgs) {
  g1 <- comp[1]
  g2 <- comp[2]
  
  b1 <- sgs_data %>% filter(Group == g1) %>% pull(b_log)
  b2 <- sgs_data %>% filter(Group == g2) %>% pull(b_log)
  
  if (length(b1) == 0 || length(b2) == 0) next
  
  obs_diff <- mean(b2, na.rm = TRUE) - mean(b1, na.rm = TRUE)
  
  boot_diffs <- replicate(n_iter, {
    mean(sample(b2, size = length(b2), replace = TRUE)) -
      mean(sample(b1, size = length(b1), replace = TRUE))
  })
  
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  significant <- !(0 >= ci[1] & 0 <= ci[2])
  
  results_sgs[[length(results_sgs) + 1]] <- data.frame(
    Group1 = g1,
    Group2 = g2,
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = significant
  )
}

final_results_sgs <- bind_rows(results_sgs)
write.csv(final_results_sgs, out_sgs_csv, row.names = FALSE)
cat("Saved SGS bootstrap results:", out_sgs_csv, "\n")

################################################################################
########################## 4) SYNTHETIC PANEL PLOT #############################
################################################################################

# Merge diversity + SGS + selfing into one display table (optional)
# Here we reload diversity CSV to keep workflow simple and reproducible
df_main <- read_csv(out_diversity_csv, show_col_types = FALSE)

# (Optional) add percent scale for visual comparison
df_main <- df_main %>%
  mutate(
    Percent_Diff = 100 * Obs_Diff / (abs(Obs_Diff) + abs(CI_upper) + 1e-6),
    Variable = factor(Variable, levels = c("AR", "He", "Ho", "Fi")),
    Category = factor(Category)
  )

plot_percent <- function(data, title_label) {
  ggplot(data, aes(x = Variable, y = Percent_Diff, fill = as.character(Significant))) +
    geom_col(width = 0.7) +
    facet_wrap(~ Category, scales = "free_x") +
    scale_fill_manual(
      values = fill_colors,
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = "Significance"
    ) +
    labs(
      title = title_label,
      x = "Genetic indicator",
      y = "Observed difference (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

df_pai74 <- df_main %>% filter(Type == "Intra-plot" & str_starts(Comparison, "PAI74"))
df_hko50 <- df_main %>% filter(Type == "Intra-plot" & str_starts(Comparison, "HKO50"))

plot2p <- plot_percent(df_pai74, "Intra-plot comparisons within PAI74")
plot3p <- plot_percent(df_hko50, "Intra-plot comparisons within HKO50")

final_panel <- plot2p / plot3p
ggsave(file.path(fig_dir, "bootstrap_synthetic_panel.png"), final_panel, width = 12, height = 10, dpi = 300, bg = "white")

cat("Saved synthetic panel:", file.path(fig_dir, "bootstrap_synthetic_panel.png"), "\n")
