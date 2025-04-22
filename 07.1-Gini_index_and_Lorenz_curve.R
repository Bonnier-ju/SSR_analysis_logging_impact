###########################################################################
##################### Gini index and Lorenz curve #########################
###########################################################################

library(tidyverse)
library(ineq)

plot_code <- "NOU"                            
plot_fullname <- "Nouragues"      

# Define input and output paths
input_file <- file.path(
  "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype",
  paste0("comparison_data_", plot_code, ".csv")
)

output_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/07-post_analysis_on_dispersal/07.1-gini_index_and_lorenz_curve"


# === COLOR SETTINGS ===
plot_colors <- list(
  "Sparouine" = "aquamarine3",
  "Paracou"   = "chocolate3",
  "NOU"       = "#548B54",       
  "PAI74"     = "mediumorchid4",
  "HKO50"     = "#CDAD00"
)

plot_color <- plot_colors[[plot_code]]

# === LOAD DATA ===
data <- read.csv(input_file, sep = ";")

# === EXTRACT ALL PARENTS ===
all_parents <- data %>%
  transmute(
    Parent1 = Parent1_Cervus,
    Parent2 = Parent2_Cervus
  ) %>%
  pivot_longer(cols = everything(), values_to = "ParentID") %>%
  filter(!is.na(ParentID))

# === COUNT OFFSPRING PER PARENT ===
parent_counts <- table(all_parents$ParentID)

# === CALCULATE GINI INDEX ===
gini_value <- Gini(parent_counts)
cat("Gini index for", plot_code, "=", round(gini_value, 4), "\n")

# === PLOT & SAVE LORENZ CURVE ===
output_file <- file.path(output_dir, paste0("lorenz_curve_", plot_code, ".png"))
png(filename = output_file, width = 800, height = 600)

plot(
  Lc(parent_counts),
  main = paste("Lorenz Curve -", plot_fullname),
  xlab = "Cumulative % of Parents",
  ylab = "Cumulative % of Offspring",
  col = plot_color,
  lwd = 2
)

dev.off()
