#############################################################################
####################### Gini index and Lorenz curve #########################
#############################################################################

# Load required packages
library(tidyverse)
library(ineq)

# === SETUP ===
plots <- list(
  list(code = "HKO50", fullname = "HKO50 (unlogged)"),
  list(code = "PAI74", fullname = "PAI74 (logged)"),
  list(code = "PAR", fullname = "Paracou"),
  list(code = "SPR", fullname = "Sparouine"),
  list(code = "NOU", fullname = "Nouragues (unlogged)")
)

plot_colors <- list(
  "Sparouine" = "aquamarine3",
  "Paracou"   = "chocolate3",
  "NOU"       = "#548B54",
  "PAI74"     = "mediumorchid4",
  "HKO50"     = "#CDAD00"
)

input_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/05-Parentage_analysis/05.6-parentage_with_haplotype"
output_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/07-post_analysis_on_dispersal/07.1-gini_index_and_lorenz_curve"


# === STEP 1: IMPORT DATA FOR ALL PLOTS ===
all_data <- list()

for (p in plots) {
  file_path <- file.path(input_dir, paste0("comparison_data_", p$code, ".csv"))
  df <- read.csv(file_path, sep = ";")
  all_data[[p$code]] <- df
}

# === FUNCTION TO COMPUTE GINI, IC, AND COUNTS ===
compute_gini_info <- function(data, plot_code) {
  all_parents <- data %>%
    transmute(
      Parent1 = Parent1_Cervus,
      Parent2 = Parent2_Cervus
    ) %>%
    pivot_longer(cols = everything(), values_to = "ParentID") %>%
    filter(!is.na(ParentID))
  
  parent_counts <- table(all_parents$ParentID)
  gini_value <- Gini(parent_counts)
  
  # Bootstrap IC
  set.seed(123)
  gini_boot <- replicate(1000, Gini(sample(parent_counts, replace = TRUE)))
  ci <- quantile(gini_boot, c(0.025, 0.975))
  
  return(list(
    code = plot_code,
    counts = parent_counts,
    gini = gini_value,
    ci_low = ci[1],
    ci_high = ci[2]
  ))
}

# === STEP 2: COMPUTE GINI FOR EACH PLOT ===
gini_results <- list()
for (p in plots) {
  gini_results[[p$code]] <- compute_gini_info(all_data[[p$code]], p$code)
}

# Export Gini summary
gini_summary <- tibble(
  Plot = sapply(plots, function(p) p$code),
  Gini = sapply(gini_results, function(x) x$gini),
  CI_low = sapply(gini_results, function(x) x$ci_low),
  CI_high = sapply(gini_results, function(x) x$ci_high)
)

write.csv(gini_summary, file.path(output_dir, "gini_summary.csv"), row.names = FALSE)

# === STEP 3: PAIRWISE PERMUTATION TESTS ===
pairwise_results <- expand.grid(Plot1 = names(gini_results), Plot2 = names(gini_results)) %>%
  filter(as.character(Plot1) < as.character(Plot2))

set.seed(123)

pairwise_results <- pairwise_results %>%
  rowwise() %>%
  mutate(
    p_value = {
      g1 <- gini_results[[Plot1]]
      g2 <- gini_results[[Plot2]]
      
      combined <- c(g1$counts, g2$counts)
      group_labels <- c(rep(1, length(g1$counts)), rep(2, length(g2$counts)))
      obs_diff <- abs(g1$gini - g2$gini)
      
      perm_diffs <- replicate(1000, {
        shuffled <- sample(group_labels)
        G1 <- Gini(combined[shuffled == 1])
        G2 <- Gini(combined[shuffled == 2])
        abs(G1 - G2)
      })
      
      mean(perm_diffs >= obs_diff)
    }
  )

write.csv(pairwise_results, file.path(output_dir, "gini_pairwise_pvalues.csv"), row.names = FALSE)

# === STEP 4: PLOT ALL LORENZ CURVES ===
png(filename = file.path(output_dir, "lorenz_curves_all_plots.png"), width = 1000, height = 800)

first_plot <- gini_results[[1]]
plot(
  Lc(first_plot$counts),
  main = "Lorenz Curves – Reproductive Success (All Plots)",
  xlab = "Cumulative % of Parents",
  ylab = "Cumulative % of Offspring",
  col = plot_colors[[as.character(first_plot$code)]],
  lwd = 2,
  lty = 1,
  ylim = c(0, 1),
  xlim = c(0, 1)
)

for (i in 2:length(gini_results)) {
  code <- names(gini_results)[i]
  lines(
    Lc(gini_results[[code]]$counts),
    col = plot_colors[[as.character(code)]],
    lwd = 2,
    lty = i
  )
}

legend("topleft",
       legend = sapply(plots, function(p) p$fullname),
       col = sapply(plots, function(p) plot_colors[[p$code]]),
       lwd = 2,
       lty = 1:length(plots),
       bty = "n",
       cex = 0.9)

dev.off()

