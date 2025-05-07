#################################################################
################ GLM - dbh - alt - nb offspring #################
#################################################################


# Load packages
library(tidyverse)

# === SETUP ===
plot_code <- "PAI74"

parentage_file <- file.path(
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo",
  paste0("comparison_data_", plot_code, ".csv")
)

treeinfo_file <- file.path(
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Data",
  plot_code,
  paste0(plot_code, "_tree_info.csv")
)

# === STEP 1: Import data ===
parentage <- read.csv(parentage_file, sep = ";")
treeinfo <- read.csv(treeinfo_file)

# === STEP 2: Extract all parent IDs ===
parents_all <- parentage %>%
  select(Parent1_Cervus, Parent2_Cervus) %>%
  pivot_longer(cols = everything(), values_to = "ParentID") %>%
  filter(!is.na(ParentID)) %>%
  count(ParentID, name = "N_offspring")

# === STEP 3: Merge with tree info ===
final_table <- parents_all %>%
  left_join(treeinfo, by = c("ParentID" = "ID")) %>%
  select(ParentID, N_offspring, dbh, alt)

# === STEP 4: Export (optional) ===
write.csv(final_table,
          file = file.path(dirname(parentage_file), paste0("parent_summary_", plot_code, ".csv")),
          row.names = FALSE)

# === OPTIONAL: View result
print(final_table)
