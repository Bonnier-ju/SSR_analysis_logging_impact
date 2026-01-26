# ============================================================
# Script: 03.2-Fij_graphs.R
# Purpose: Plot kinship coefficient (Fij) by distance class with CI ribbon and b-log annotation.
# Manuscript: Logging impact (HKO50 vs PAI74) – SSR analyses
# Inputs: results CSV exported from SGS/kinship analysis (Fij_<site>_<suffix>.csv)
# Outputs: figures/03_fij/Fij_<site>_<suffix>.png
# Figures/Tables: Kinship-by-distance plots (SGS section / supplementary)
# Last update: 2026-01-26
# ============================================================

## ===================== PARAMETERS =====================
site_id <- "PAI74"
site_title <- "PAI74"
suffix <- "SED_inopening"   # e.g. "FULL", "SED", "INT", "ADL", "SED_inopening"...
line_color <- "brown2"

base_path <- "."  # project root

input_dir <- file.path(base_path, "results", "03_fij")   # adapt if your CSVs are elsewhere
fig_dir <- file.path("figures", "03_fij")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(input_dir, paste0("Fij_", site_id, "_", suffix, ".csv"))
output_file <- file.path(fig_dir, paste0("Fij_", site_id, "_", suffix, ".png"))

# Distance classes and plotting positions (pseudo-log spacing)
distance_classes <- c("0", "30", "60", "90", "130", "170", "220", "300", "600")
x_plot <- c(1, 3, 4.3, 5.3, 6.3, 7.1, 8, 9, 11.2)

# Y limits
ylim_vals <- c(-0.09, 0.09)
## ======================================================

## ===================== LIBRARIES ======================
library(dplyr)
library(ggplot2)
## ======================================================

## ===================== INPUT ==========================
stopifnot(file.exists(input_file))
df <- read.csv(input_file, sep = ",", header = FALSE, stringsAsFactors = FALSE)
## ======================================================

## ===================== EXTRACT VALUES =================
# b-log position in your exported file (keep as in your original format)
b_log <- round(as.numeric(df[6, 14]), 4)

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
    distance_m = as.numeric(distance_class),
    x_plot = x_plot
  )
## ======================================================

## ===================== PLOT ===========================
p <- ggplot() +
  geom_ribbon(
    data = subset(kinship_df, distance_m != 0),
    aes(x = x_plot, ymin = CI_inf, ymax = CI_sup),
    fill = "#8B8B83", alpha = 0.4
  ) +
  geom_line(
    data = subset(kinship_df, distance_m != 0),
    aes(x = x_plot, y = obs_val),
    color = line_color, size = 1
  ) +
  geom_point(
    data = subset(kinship_df, distance_m != 0),
    aes(x = x_plot, y = obs_val),
    size = 3, color = line_color
  ) +
  geom_point(
    data = subset(kinship_df, distance_m == 0),
    aes(x = x_plot, y = obs_val),
    size = 4, shape = 18, color = "red"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(
    breaks = kinship_df$x_plot,
    labels = kinship_df$distance_class
  ) +
  coord_cartesian(ylim = ylim_vals) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("b-log = ", b_log),
    hjust = 1.2, vjust = 2.2, size = 5
  ) +
  labs(
    title = paste(site_title, suffix, sep = " - "),
    x = "Distance class (m)",
    y = "Kinship coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15)
  )

print(p)

ggsave(
  filename = output_file,
  plot = p,
  width = 7, height = 5, dpi = 300, bg = "white"
)

cat("Saved figure:", output_file, "\n")
## ======================================================
