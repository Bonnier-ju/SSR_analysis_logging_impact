################################################################################
############################ STRUCTURE results visualization ###################
################################################################################

library(readr)
library(ggplot2)

# -------- File paths --------
evanno_file <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/results_harvester_PAI74_HKO50/evanno.txt"
output_dir <- "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/02-populations_structure/02.1-population_structure/02.1.1-STRUCTURE/results_harvester_PAI74_HKO50"

# -------- Read and clean data --------
col_names <- c("K", "Reps", "Mean_LnP_K", "Stdev_LnP_K", "Ln_prime_K", "Abs_Ln_doubleprime_K", "Delta_K")
evanno_data <- read_tsv(evanno_file, comment = "#", col_names = col_names)

# -------- Plot Delta K --------
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

ggsave(filename = file.path(output_dir, "evanno_plot.png"), 
       plot = delta_plot, width = 7, height = 5, dpi = 300, bg = "white")

# -------- Plot Ln′(K) --------
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

ggsave(filename = file.path(output_dir, "lnprime_plot.png"), 
       plot = lnprime_plot, width = 7, height = 5, dpi = 300, bg = "white")
