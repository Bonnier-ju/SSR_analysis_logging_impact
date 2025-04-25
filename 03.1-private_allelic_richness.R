#############################################################################
########  Private allelic richness (pAR) – by plot and by DBH cohort ########
#############################################################################
# Private allelic richness (pAR) is the average number of alleles found exclusively in one group 
# (i.e. not shared with others), corrected for unequal sample sizes through rarefaction. 
# It represents the unique genetic reservoir of that group.

library(tidyverse)
library(poppr)
library(adegenet)


#Parameter
input_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Analysis/03-diversity_and_SGS_analysis/privateAR/Full_plot_infos.csv"

plot_name <- "HKO50"       # Nom du plot à analyser
dbh_range <- c(0, 2)       # Intervalle DBH à analyser (ex : 5–30)

data <- read.csv(input_file)

#filtering data
filtered_data <- data %>%
  filter(plot == plot_name, dbh >= dbh_range[1], dbh < dbh_range[2])

allele_cols <- grep("_[12]$", names(filtered_data), value = TRUE)

# Convert to genind
genind_obj <- df2genind(
  filtered_data[, allele_cols],
  ploidy = 2,
  sep = ",",
  ind.names = filtered_data$ID,
  pop = as.factor(rep(paste0(plot_name, "_", paste(dbh_range, collapse = "-")), nrow(filtered_data)))
)

# Computing private allele
pa <- private_alleles(genind_obj, report = "table", level = "population", count.alleles = TRUE, drop = FALSE)
total_pa <- sum(pa > 0, na.rm = TRUE)



# === 7. AFFICHAGE DU RÉSULTAT ===
cat("Plot :", plot_name, "\n")
cat("DBH range :", paste(dbh_range, collapse = "–"), "\n")
cat("Private Allelic Richness (pAR) :", total_pa, "\n")
