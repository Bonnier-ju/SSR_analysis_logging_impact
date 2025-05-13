################################################################################
###################### Bootstrapping genetic index #############################
################################################################################

library(dplyr)
library(readr)
library(purrr)
library(stringr)

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



# Define valid biological categories
biological_cats <- c("SED", "INT", "ADL")

# Function to check if a comparison is between two biological categories only
valid_intra_biological <- function(category_string) {
  parts <- unlist(strsplit(category_string, " vs "))
  length(parts) == 2 && all(parts %in% biological_cats)
}

# Filter the results
filtered_results <- final_results %>%
  filter(
    # Keep all inter-plot comparisons
    Type == "Inter-plot" |
      
      # Keep intra-plot comparisons only if both compared groups are in SED, INT, ADL
      (Type == "Intra-plot" & sapply(Category, valid_intra_biological))
  )


# Add a new column to classify comparison group
filtered_results <- filtered_results %>%
  mutate(
    ComparisonGroup = case_when(
      Type == "Inter-plot" ~ "HKO50 vs PAI74",
      Type == "Intra-plot" & str_starts(Comparison, "HKO50") ~ "HKO50 only",
      Type == "Intra-plot" & str_starts(Comparison, "PAI74") ~ "PAI74 only",
      TRUE ~ "Other"
    )
  )




# --- 4. Prepare three datasets ---
# Graph 1: Inter-plot comparisons
data_inter <- filtered_results %>%
  filter(Type == "Inter-plot")

# Graph 2: Intra-plot comparisons for PAI74
data_pai <- filtered_results %>%
  filter(Type == "Intra-plot" & str_starts(Comparison, "PAI74"))

# Graph 3: Intra-plot comparisons for HKO50
data_hko <- filtered_results %>%
  filter(Type == "Intra-plot" & str_starts(Comparison, "HKO50"))

# --- 5. Define plotting function ---
plot_comparisons <- function(df, title_label) {
  df$Variable <- factor(df$Variable, levels = c("AR", "He", "Ho", "Fi"))
  
  ggplot(df, aes(x = Variable, y = Obs_Diff, fill = Significant)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    facet_wrap(~Comparison, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey70")) +
    labs(
      title = title_label,
      y = "Observed Difference",
      x = "Genetic Indicator"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create a more informative label for inter-plot comparisons
data_inter <- data_inter %>%
  mutate(Comparison_detail = paste(Category, "-", Comparison))

plot1 <- ggplot(data_inter, aes(x = Variable, y = Obs_Diff, fill = Significant)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  facet_wrap(~Comparison_detail, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey70")) +
  labs(
    title = "PAI74 vs HKO50 comparisons by category",
    y = "Observed Difference",
    x = "Genetic Indicator"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- 6. Create plots ---
plot1 <- plot_comparisons(data_inter, "PAI74 vs HKO50 comparisons")
plot2 <- plot_comparisons(data_pai, "Intra-plot comparisons in PAI74")
plot3 <- plot_comparisons(data_hko, "Intra-plot comparisons in HKO50")

# --- 7. Display plots ---
print(plot1)
print(plot2)
print(plot3)
