################################################################################
############################ STRUCTURE results visualization ###################
################################################################################

library(readr)
library(ggplot2)

# File paths 
evanno_file <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/results_harvester_PAI74_HKO50/evanno.txt"
output_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/results_harvester_PAI74_HKO50"

# Read and clean data 
col_names <- c("K", "Reps", "Mean_LnP_K", "Stdev_LnP_K", "Ln_prime_K", "Abs_Ln_doubleprime_K", "Delta_K")
evanno_data <- read_tsv(evanno_file, comment = "#", col_names = col_names)

# Plot Delta K 
delta_plot <- ggplot(evanno_data, aes(x = K, y = Delta_K)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 3, color = "steelblue") +
  theme_minimal() +
  labs(title = "Evanno Method - ΔK vs K",
       x = "K", y = "ΔK") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(filename = file.path(output_dir, "evanno_plot_delta_K.png"), 
       plot = delta_plot, width = 7, height = 5, dpi = 300, bg = "white")

# Plot Ln′(K)
lnprime_plot <- ggplot(evanno_data, aes(x = K, y = Ln_prime_K)) +
  geom_line(color = "darkorange", size = 1) +
  geom_point(size = 3, color = "darkorange") +
  theme_minimal() +
  labs(title = "Evanno Method - Ln′(K) vs K",
       x = "K", y = "Ln′(K)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(filename = file.path(output_dir, "ln_plot.png"), 
       plot = lnprime_plot, width = 7, height = 5, dpi = 300, bg = "white")



################################################################################
########################### Barplots visualisation #############################
################################################################################

# Parameter
k_value <- 8
base_path <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact"
indfile_path <- file.path(base_path, "Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/results_harvester_PAI74_HKO50", paste0("K", k_value, ".indfile"))
str_path <- file.path(base_path, "Data/HKO50_PAI74/PAI74_HKO50.str.txt")
meta_path <- file.path(base_path, "Data/HKO50_PAI74/PAI74_HKO50_tree_info.csv")
output_dir <- file.path(base_path, "Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/plots_structure_PAI74_HKO50")
dir.create(output_dir, showWarnings = FALSE)

# Read STR file to extract individual IDs
str_data <- read_table(str_path, skip = 1, col_names = FALSE)
individual_ids <- str_data$X1[seq(1, nrow(str_data), by = 2)]

# Read Harvester .indfile (first replicate only)
indfile <- read_lines(indfile_path)[1:length(individual_ids)]
ind_probs <- indfile %>%
  str_trim() %>%
  str_remove("^\\d+\\s+\\d+\\s+\\(\\d+\\)\\s+\\d+\\s+:\\s+") %>%
  str_split_fixed(" ", k_value) %>%
  as_tibble() %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ . / rowSums(across(everything())))) %>%
  mutate(ID = individual_ids)


# Rename cluster columns
colnames(ind_probs)[1:k_value] <- paste0("Cluster", 1:k_value)

# Read metadata and assign size class
meta <- read_csv(meta_path)
meta <- meta %>%
  mutate(cohort = case_when(
    dbh < 5 ~ "SED",
    dbh < 30 ~ "INT",
    dbh >= 30 ~ "ADL"
  ))

# Add dominant cluster info
ind_probs <- ind_probs %>%
  rowwise() %>%
  mutate(
    max_prob = max(c_across(starts_with("Cluster"))),
    dominant_cluster = names(.)[which.max(c_across(starts_with("Cluster")))]
  ) %>%
  ungroup()

# Merge and sort individuals by plot, cohort, then dominant cluster and descending max prob
plot_data <- ind_probs %>%
  left_join(meta, by = "ID") %>%
  arrange(plot, cohort, dominant_cluster, desc(max_prob)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Probability") %>%
  mutate(Cluster = factor(Cluster, levels = paste0("Cluster", 1:k_value)))  # fixed stack order

# Plot
plot_structure <- ggplot(plot_data, aes(x = ID, y = Probability, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ plot + cohort, scales = "free_x", space = "free") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "top") +
  labs(title = paste("STRUCTURE barplot - K =", k_value),
       x = "Individuals", y = "Assignment Probability")
plot_structure

ggsave(filename = file.path(output_dir, paste0("structure_barplotPAI74_HKO50_K", k_value, ".png")),
       plot = plot_structure, width = 19, height = 5, dpi = 300, bg = "white")





