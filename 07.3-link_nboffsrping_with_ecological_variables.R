#################################################################
################ GLM - dbh - alt - nb offspring #################
#################################################################

library(MASS)
library(ggplot2)
library(tidyverse)

# === SETUP ===
plot_code <- "Nouragues"

parentage_file <- file.path(
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Analysis/04-parentage_analysis/04.3-inferring_mother_haplo",
  paste0("comparison_data_", plot_code, ".csv")
)

treeinfo_file <- file.path(
  "C:/Users/bonni/OneDrive/University/Thesis/Dicorynia/Article-Logging_impact/Data",
  plot_code,
  paste0(plot_code, "_tree_info.csv")
)

# Import data
parentage <- read.csv(parentage_file)
treeinfo <- read.csv(treeinfo_file)

#Extract all parent IDs
parents_all <- parentage %>%
  select(Parent1_Cervus, Parent2_Cervus) %>%
  pivot_longer(cols = everything(), values_to = "ParentID") %>%
  filter(!is.na(ParentID)) %>%
  count(ParentID, name = "N_offspring")

#Merge with tree info
final_table <- parents_all %>%
  left_join(treeinfo, by = c("ParentID" = "ID")) %>%
  select(ParentID, N_offspring, dbh, alt)


# Descriptive stats
summary(final_table)
str(final_table)

#Visual check of distributions
par(mfrow = c(1, 3))
hist(final_table$N_offspring, main = "Offspring count", col = "lightblue")
hist(final_table$dbh, main = "DBH", col = "lightgreen")
hist(final_table$alt, main = "Altitude", col = "lightpink")
par(mfrow = c(1, 1))

#Check mean vs variance (Poisson assumption) 
cat("Mean and variance of offspring count:\n")
print(mean(final_table$N_offspring))
print(var(final_table$N_offspring))


# Fit negative binomial GLM
glm_nb <- glm.nb(N_offspring ~ dbh + alt, data = final_table)
summary(glm_nb)


# Interpret coefficients
cat("\nExponentiated coefficients (effect on offspring count per unit increase):\n")
exp_coef <- exp(coef(glm_nb))
print(round(exp_coef, 3))

cat("\nConfidence intervals (95%) for exponentiated coefficients:\n")
ci <- confint(glm_nb)
exp_ci <- exp(ci)
print(round(exp_ci, 3))


# DBH effect plot
ggplot(final_table, aes(x = dbh, y = N_offspring)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = MASS::glm.nb, formula = y ~ x, se = TRUE, color = "darkblue") +
  labs(title = paste0("Effect of DBH on reproductive success (", plot_code, ")"),
       x = "DBH (cm)", y = "Number of offspring") +
  theme_minimal()

# Altitude effect plot
ggplot(final_table, aes(x = alt, y = N_offspring)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = MASS::glm.nb, formula = y ~ x, se = TRUE, color = "darkgreen") +
  labs(title = paste0("Effect of altitude on reproductive success (", plot_code, ")"),
       x = "Altitude (m)", y = "Number of offspring") +
  theme_minimal()


