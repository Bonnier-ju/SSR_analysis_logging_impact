################################################################################
##################### Bootstrap Analysis of Genetic Indices ####################
################################################################################

# Load required packages
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load and combine input data
# ------------------------------------------------------------------------------

# Define input file paths
file_pai <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/PAI74/results_PAI74_full_bootstrap.csv"
file_hko <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/HKO50/results_HKO50_full_bootstrap.csv"

# Read data
pai <- read_csv(file_pai, show_col_types = FALSE) %>% mutate(Plot = "PAI74")
hko <- read_csv(file_hko, show_col_types = FALSE) %>% mutate(Plot = "HKO50")

# Combine into one dataset
data_all <- bind_rows(pai, hko)

# ------------------------------------------------------------------------------
# 2. Bootstrap function
# ------------------------------------------------------------------------------

bootstrap_diff <- function(df1, df2, var) {
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
  
  diffs_boot <- replicate(10000, mean(sample(merged$Diff, replace = TRUE)))
  ci <- quantile(diffs_boot, c(0.025, 0.975), na.rm = TRUE)
  obs_diff <- mean(merged$Diff)
  
  data.frame(
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = !(0 >= ci[1] & 0 <= ci[2])
  )
}

# ------------------------------------------------------------------------------
# 3. Compute bootstrap comparisons
# ------------------------------------------------------------------------------

set.seed(123)
indicators <- c("AR", "He", "Ho", "Fi")
plots <- unique(data_all$Plot)
categories <- unique(data_all$Category)
results_list <- list()

# Inter-plot comparisons (same category between plots)
for (cat in categories) {
  for (var in indicators) {
    g1 <- data_all %>% filter(Plot == "HKO50", Category == cat)
    g2 <- data_all %>% filter(Plot == "PAI74", Category == cat)
    
    res <- bootstrap_diff(g1, g2, var)
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
      
      res <- bootstrap_diff(g1, g2, var)
      res$Comparison <- paste0(plot, ": ", pair[2], " vs ", pair[1])
      res$Category <- paste(pair[1], "vs", pair[2])
      res$Variable <- var
      res$Type <- "Intra-plot"
      results_list[[length(results_list) + 1]] <- res
    }
  }
}

# Combine all results
final_results <- bind_rows(results_list) %>%
  select(Type, Comparison, Category, Variable, Obs_Diff, CI_lower, CI_upper, Significant)

# ------------------------------------------------------------------------------
# 4. Filter biologically relevant comparisons
# ------------------------------------------------------------------------------

biological_cats <- c("SED", "INT", "ADL")

valid_intra_biological <- function(category_string) {
  parts <- unlist(strsplit(category_string, " vs "))
  length(parts) == 2 && all(parts %in% biological_cats)
}

filtered_results <- final_results %>%
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

# Export results
write.csv(
  filtered_results,
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/bootstrap_comparisons_inter_intra_plots.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Prepare data for plotting
# ------------------------------------------------------------------------------

# Split into 3 comparison types
data_inter <- filtered_results %>% filter(ComparisonGroup == "HKO50 vs PAI74") %>%
  mutate(Comparison_detail = paste(Category, "-", Comparison))
data_pai <- filtered_results %>% filter(ComparisonGroup == "PAI74 only")
data_hko <- filtered_results %>% filter(ComparisonGroup == "HKO50 only")

# ------------------------------------------------------------------------------
# 6. Plotting function
# ------------------------------------------------------------------------------

plot_comparisons <- function(df, title_label, use_comparison_col = TRUE) {
  df$Variable <- factor(df$Variable, levels = c("AR", "He", "Ho", "Fi"))
  facet_var <- if (use_comparison_col) "Comparison_detail" else "Comparison"
  
  ggplot(df, aes(x = Variable, y = Obs_Diff, fill = Significant)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey70")) +
    labs(
      title = title_label,
      y = "Observed Difference",
      x = "Genetic Indicator"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------------------
# 7. Generate and display plots
# ------------------------------------------------------------------------------

plot1 <- plot_comparisons(data_inter, "PAI74 vs HKO50 comparisons", use_comparison_col = TRUE)
plot2 <- plot_comparisons(data_pai, "Intra-plot comparisons in PAI74", use_comparison_col = FALSE)
plot3 <- plot_comparisons(data_hko, "Intra-plot comparisons in HKO50", use_comparison_col = FALSE)

print(plot1)
print(plot2)
print(plot3)





################################################################################
######################## Bootstrap Analysis of selfing #########################
################################################################################

# Load required libraries
library(dplyr)
library(readr)
library(tidyr)

# Set file path
file_path <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/PAI74_HKO50/results_HKO50_full_bootstrap_selfing.csv"

# Read data
selfing_data <- read_csv(file_path, show_col_types = FALSE)


marker_cols <- setdiff(names(selfing_data), c("Plot", "Population", "S(selfing rate)"))
head(marker_cols)


selfing_data <- selfing_data %>%
  mutate(Group = paste(Plot, Population, sep = "_"))
groups <- unique(selfing_data$Group)
group_comparisons <- combn(groups, 2, simplify = FALSE)


# Number of bootstrap iterations
n_iter <- 10000
set.seed(123)
results <- list()

# Loop over all pairwise group comparisons
for (comp in group_comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]
  
  # Extract selfing data for each group
  m1 <- selfing_data %>% filter(Group == g1)
  m2 <- selfing_data %>% filter(Group == g2)
  
  # Skip comparison if one group is missing or has multiple rows
  if (nrow(m1) != 1 | nrow(m2) != 1) next
  
  # Perform bootstrap by resampling marker columns with replacement
  boot_diffs <- replicate(n_iter, {
    sampled_cols <- sample(marker_cols, length(marker_cols), replace = TRUE)
    mean(as.numeric(m2[1, sampled_cols]), na.rm = TRUE) -
      mean(as.numeric(m1[1, sampled_cols]), na.rm = TRUE)
  })
  
  # Compute observed multilocus selfing difference
  obs_diff <- mean(as.numeric(m2[1, marker_cols]), na.rm = TRUE) -
    mean(as.numeric(m1[1, marker_cols]), na.rm = TRUE)
  
  # Compute 95% confidence interval and test for significance
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  significant <- !(0 >= ci[1] & 0 <= ci[2])
  
  # Store results
  results[[length(results) + 1]] <- data.frame(
    Group1 = g1,
    Group2 = g2,
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = significant
  )
}


# Combine all results into a single dataframe
final_results <- bind_rows(results)

print(final_results)

# Export results to CSV
write.csv(
  final_results,
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/final_selfing_comparisons.csv",
  row.names = FALSE
)


################################################################################
########################## Bootstrap Analysis of SGS ###########################
################################################################################

library(dplyr)
library(readr)

# Load your CSV file
file_path <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/PAI74_HKO50/results_HKO50_PAI74_Sp_sign_btw_plot.csv"
data <- read_csv(file_path, show_col_types = FALSE)

# Create a unique group identifier
data <- data %>%
  mutate(Group = paste(Plot, category, sep = "_"))

# Get list of all unique group names
groups <- unique(data$Group)
group_comparisons <- combn(groups, 2, simplify = FALSE)

# Initialize bootstrap parameters
n_iter <- 10000
set.seed(123)
results <- list()

# Loop over all group comparisons
for (comp in group_comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]
  
  b1 <- data %>% filter(Group == g1) %>% pull(b_log)
  b2 <- data %>% filter(Group == g2) %>% pull(b_log)
  
  # Skip if one group is empty
  if (length(b1) == 0 | length(b2) == 0) next
  
  # Compute observed difference
  obs_diff <- mean(b2, na.rm = TRUE) - mean(b1, na.rm = TRUE)
  
  # Bootstrap
  boot_diffs <- replicate(n_iter, {
    mean(sample(b2, size = length(b2), replace = TRUE)) -
      mean(sample(b1, size = length(b1), replace = TRUE))
  })
  
  # Confidence interval and significance
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  significant <- !(0 >= ci[1] & 0 <= ci[2])
  
  # Store result
  results[[length(results) + 1]] <- data.frame(
    Group1 = g1,
    Group2 = g2,
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = significant
  )
}

# Combine all results
final_results <- bind_rows(results)

# View results
print(final_results)

# Optional: export to CSV
write.csv(
  final_results,
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/bootstrap_b_log_comparisons.csv",
  row.names = FALSE
)















