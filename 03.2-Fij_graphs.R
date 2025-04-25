
#################### Visualization of Fij by distance class ###################

library(tidyverse)
library(ggplot2)
library(dplyr)


### Files parameters ###
site_id <- "PAI74"         
site_title <- "PAI74"   
suffix <- "BIG"         
line_color <- "mediumorchid4"

base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/03-diversity_and_SGS_analysis"


input_file <- file.path(base_path, paste0("Fij_", site_id, "_", suffix, ".csv"))
output_file <- file.path(base_path, paste0("Fij_", site_id, "_", suffix, ".png"))

#input_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/03-diversity_and_SGS_analysis/Fij_HKO50_FULL.csv"
#output_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/03-diversity_and_SGS_analysis/PAI74_fij_graphs.png"

### Import data ###
df <- read.csv(input_file, sep = ";", header = FALSE)

b_log <- round(as.numeric(df[6, 14]), 4)
distance_classes <- c("0", "30", "60", "90", "130", "170", "220", "300", "600")
x_plot <- c(1, 3, 4.3, 5.3, 6.3, 7.1, 8, 9, 11.2)

obs_val <- as.numeric(df[6, 2:10])
CI_inf  <- as.numeric(df[9, 2:10])
CI_sup  <- as.numeric(df[10, 2:10])

kinship_df <- data.frame(
  distance_class = distance_classes,
  obs_val = obs_val,
  CI_inf = CI_inf,
  CI_sup = CI_sup
) %>%
  mutate(
    distance_m = as.numeric(as.character(distance_class)),
    x_plot = x_plot
  )

### Create plot ###
plot <- ggplot() +
  geom_ribbon(data = subset(kinship_df, distance_m != 0),
              aes(x = x_plot, ymin = CI_inf, ymax = CI_sup),
              fill = "#8B8B83", alpha = 0.4) +
  geom_line(data = subset(kinship_df, distance_m != 0),
            aes(x = x_plot, y = obs_val),
            color = line_color, size = 1) +
  geom_point(data = subset(kinship_df, distance_m != 0),
             aes(x = x_plot, y = obs_val),
             size = 3, color = line_color) +
  geom_point(data = subset(kinship_df, distance_m == 0),
             aes(x = x_plot, y = obs_val),
             size = 4, shape = 18, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(
    breaks = kinship_df$x_plot,
    labels = kinship_df$distance_class
  ) +
  coord_cartesian(ylim = c(-0.08, 0.08)) +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("b-log = ", b_log)),
    hjust = 1.5, vjust = 3, size = 4
  ) +
  labs(
    title = "PAI74 - BIG",
    x = "Distance class (m)",
    y = "Kinship coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

print(plot)

### Save plot to file ###
ggsave(
  filename = output_file,
  plot = plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)
