################################################################################
###################### Bootstrapping genetic index #############################
################################################################################

library(dplyr)
library(ggplot2)

# Define file paths
file_pai <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/PAI74/results_PAI74_full_bootstrap.csv"
file_hko <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/03-diversity_and_SGS_analysis/spagedie/Results/HKO50/results_HKO50_full_bootstrap.csv"

# Read the input files
pai <- read.csv(file_pai)
hko <- read.csv(file_hko)

# Add a column identifying the plot
pai$Plot <- "PAI74"
hko$Plot <- "HKO50"

# Combine the two datasets
data_all <- bind_rows(pai, hko)

# Check available column names
print(colnames(data_all))

# List of genetic diversity indicators to compare
indicators <- c("AR", "He", "Ho", "Fi")

# Create a list to store the bootstrap results
bootstrap_results <- list()

# Number of bootstrap iterations
n_iter <- 10000
set.seed(123)  # Ensure reproducibility

# Loop through each indicator
for (ind in indicators) {
  
  # Merge the two plots by marker for the current indicator
  merged <- inner_join(
    hko %>% select(Marker, !!sym(ind)) %>% rename(HKO50 = !!sym(ind)),
    pai %>% select(Marker, !!sym(ind)) %>% rename(PAI74 = !!sym(ind)),
    by = "Marker"
  )
  
  # Check that we have enough shared markers
  if (nrow(merged) < 5) {
    warning(paste("Not enough shared markers for", ind))
    next
  }
  
  # Compute the difference (PAI74 - HKO50) per marker
  merged <- merged %>%
    mutate(Diff = PAI74 - HKO50)
  
  # Perform bootstrap resampling: resample the marker differences with replacement
  diffs_boot <- replicate(n_iter, {
    sample(merged$Diff, size = nrow(merged), replace = TRUE) %>% mean()
  })
  
  # Compute 95% confidence interval of the bootstrap distribution
  ci <- quantile(diffs_boot, probs = c(0.025, 0.975))
  obs_diff <- mean(merged$Diff)
  
  # Store the result for this indicator
  bootstrap_results[[ind]] <- data.frame(
    Indicator = ind,
    Obs_Diff = obs_diff,
    CI_lower = ci[1],
    CI_upper = ci[2],
    Significant = !(0 >= ci[1] & 0 <= ci[2])  # TRUE if 0 is not in the CI
  )
}

# Combine all indicator results into a final summary table
final_results <- bind_rows(bootstrap_results)

# Print the results
print(final_results)

# Optionally save to CSV
write.csv(final_results, "bootstrap_diff_HKO50_vs_PAI74.csv", row.names = FALSE)
