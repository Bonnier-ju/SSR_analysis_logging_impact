################################################################################
###################### Bootstrapping genetic index #############################
################################################################################

library(dplyr)
library(readr)
library(purrr)

# Define file paths
file_pai <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/PAI74/results_PAI74_full_bootstrap.csv"
file_hko <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/HKO50/results_HKO50_full_bootstrap.csv"

# Load CSV files
pai <- read_csv(file_pai, show_col_types = FALSE)
hko <- read_csv(file_hko, show_col_types = FALSE)

# Add plot column
pai$Plot <- "PAI74"
hko$Plot <- "HKO50"

# Combine all data
data_all <- bind_rows(pai, hko)

# Define variables to test
indicators <- c("AR", "He", "Ho", "Fi")
plots <- unique(data_all$Plot)
categories <- unique(data_all$Category)
n_iter <- 10000
set.seed(123)

# Bootstrap function
bootstrap_diff <- function(df1, df2, var) {
  merged <- inner_join(
    df1 %>% select(Markers, value1 = !!sym(var)),
    df2 %>% select(Markers, value2 = !!sym(var)),
    by = "Markers"
  ) %>%
    mutate(Diff = value2 - value1) %>%
    filter(!is.na(Diff))
  
  if (nrow(merged) < 5) {
    return(data.frame(
      Obs_Diff = NA, CI_lower = NA, CI_upper = NA, Significant = NA
    ))
  }
  
  diffs_boot <- replicate(n_iter, {
    mean(sample(merged$Diff, replace = TRUE))
  })
  
  ci <- quantile(diffs_boot, c(0.025, 0.975), na.rm = TRUE)
  obs_diff <- mean(merged$Diff)
  
  data.frame(
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = !(0 >= ci[1] & 0 <= ci[2])
  )
}

# Store results
all_results <- list()

# 1. Inter-plot comparisons: same category between PAI74 and HKO50
for (cat in categories) {
  for (var in indicators) {
    group1 <- data_all %>% filter(Plot == "HKO50", Category == cat)
    group2 <- data_all %>% filter(Plot == "PAI74", Category == cat)
    
    result <- bootstrap_diff(group1, group2, var)
    result$Comparison <- paste0("PAI74 vs HKO50")
    result$Category <- cat
    result$Variable <- var
    result$Type <- "Inter-plot"
    
    all_results[[length(all_results) + 1]] <- result
  }
}

# 2. Intra-plot comparisons: different categories within each plot
for (plot in plots) {
  sub_data <- data_all %>% filter(Plot == plot)
  combs <- combn(unique(sub_data$Category), 2, simplify = FALSE)
  
  for (comb in combs) {
    cat1 <- comb[1]
    cat2 <- comb[2]
    
    for (var in indicators) {
      group1 <- sub_data %>% filter(Category == cat1)
      group2 <- sub_data %>% filter(Category == cat2)
      
      result <- bootstrap_diff(group1, group2, var)
      result$Comparison <- paste0(plot, ": ", cat2, " vs ", cat1)
      result$Category <- paste(cat1, "vs", cat2)
      result$Variable <- var
      result$Type <- "Intra-plot"
      
      all_results[[length(all_results) + 1]] <- result
    }
  }
}

# Combine all results
final_results <- bind_rows(all_results) %>%
  select(Type, Comparison, Category, Variable, Obs_Diff, CI_lower, CI_upper, Significant)

# Print final result
print(final_results)

# Export to CSV in specified folder
write.csv(
  final_results,
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/bootstrap_comparisons_inter_intra_plots.csv",
  row.names = FALSE
)
