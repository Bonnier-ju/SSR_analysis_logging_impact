
################################################################################
############### Checking contamination with Rclone #############################
################################################################################

#Load and install packages
install.packages("devtools")
library(devtools)

devtools::install_github("dbailleul/RClone")

library(RClone)
library(dplyr)
library(pegas)

#Load dataset
data <- read.csv("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Logging impact/Data/Full_data_HKO50.csv", sep=";", header=TRUE)

indiv_ids <- data$ID
markers <- data %>% select(starts_with("SSRSeq"))
markers[is.na(markers)] <- "000"

# Identify Unique Multilocus Genotypes (MLGs)
# A Multilocus Genotype (MLG) is a unique genetic profile across multiple SSR loci.
# If two individuals share the exact same MLG, they are considered clones or duplicates.
# If many individuals have identical MLGs, it may indicate contamination or sample misidentification

list_all_tab(markers)

MLG_tab(markers, vecpop = NULL)

MLGlist <- MLG_list(markers)
for (i in seq_along(MLGlist)) {
  cat("MLG", i, "contains individuals (IDs):", indiv_ids[MLGlist[[i]]], "\n")
}

# Allelic frequencies per locus
freq_RR(markers)

# Evaluate the reliability and informativeness of the loci (markers)
# Perform a resampling procedure to measure how diversity indices (e.g., heterozygosity, number of multilocus genotypes) 
# change as the number of sampled loci increases.
# also present results with a genotype accumulation curves 
# plot : 800x600
res <- sample_loci(markers, nbrepeat = 100, He = TRUE, graph = TRUE, bar = TRUE)

res$res_MLG

res$res_alleles

sample_units(markers, nbrepeat = 100)

################## Determination of MLL (MultiLocus Lineages) ##################
# Multilocus Lineage, which refers to groups of individuals sharing very similar multilocus genotypes, 
# but not necessarily exactly identical ones

# probability of observing a given multilocus genotype (MLG) by chance, given the observed allele frequencies in your dataset
pgen(markers)

# the probability that a given genotype appears more than once purely by sexual reproduction (without clonal reproduction)
# helps to determine if the repeated occurrence of a genotype indicates: Clonal reproduction, Sample contamination, Sampling errors
# Very low psex values mean the repeated occurrence of a genotype is unlikely by chance alone, suggesting clonal reproduction or contamination

res_psex <- psex(markers, RR = TRUE, nbrepeat = 100)
res_psex


############ Tests for MLLs occurrence and assessment of their memberships ############
# Calculate genetic distances between individuals to evaluate genetic similarity.
# Define Multilocus Lineages (MLLs) by grouping individuals genetically similar beyond a defined threshold (parameter alpha2).
# Identify groups of individuals (MLLs) potentially resulting from clonality, contamination, or genotyping errors.

#genetic distances computation, distance on allele differences:
# observed pairwise genetic distances between individuals based on allele differences
respop <- genet_dist(markers)

#  theoretical distribution of genetic distances based on simulations
ressim <- genet_dist_sim(markers, nbrepeat = 100) #theoretical distribution

# simulated theoretical distribution but explicitly excluding self-fertilization
ressimWS <- genet_dist_sim(markers, genet = TRUE, nbrepeat = 100) #idem, without selfing

#graph prep.:
p1 <- hist(respop$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), main = "Regina",
           xlab = "Genetic distances", breaks = seq(0, max(respop$distance_matrix)+1, 1))


p2 <- hist(ressim$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), main = "Regina",
           xlab = "Genetic distances", breaks = seq(0, max(ressim$distance_matrix)+1, 1))


p3 <- hist(ressimWS$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3),
           main = "popSRWS", xlab = "Genetic distances",
           breaks = seq(0, max(ressimWS$distance_matrix)+1, 1))


# Définir limite commune sur X
limx <- max(max(respop$distance_matrix), max(ressim$distance_matrix), max(ressimWS$distance_matrix))

# Tracer les histogrammes avec ylim étendu jusqu'à 0.1
hist(respop$distance_matrix, freq = FALSE, col = "mediumorchid4",
     main = "Genetic distance distributions - Regina",
     xlab = "Genetic distances",
     breaks = seq(0, limx+1, 1), 
     xlim = c(0, 100),
     ylim = c(0, 0.20))  # étendre l'axe Y ici à 0.1

hist(ressim$distance_matrix, freq = FALSE, add = TRUE, 
     col = rgb(0.7, 0.9, 1, 0.5), breaks = seq(0, limx+1, 1))

hist(ressimWS$distance_matrix, freq = FALSE, add = TRUE, 
     col = rgb(0.9, 0.5, 1, 0.3), breaks = seq(0, limx+1, 1))

# Placer manuellement la légende avec des coordonnées pour éviter la superposition :
legend(x = limx*0.1, y = 0.20,  # Ajuste précisément ces valeurs en fonction du graphique
       legend = c("Observed data", "Simulated data", "Simulated (no selfing)"),
       fill = c("mediumorchid4", rgb(0.7, 0.9, 1, 0.5), rgb(0.9, 0.5, 1, 0.3)),
       bg = "white", box.lwd = 1, cex = 0.8)

# plot 700x600
# frequency table displaying the number of pairs of individuals in your dataset that share a 
# specific genetic distance (number of allelic differences)
# The first row lists observed genetic distances (number of allele differences).
# The second row shows the number of individual pairs having each specific distance
table(respop$distance_matrix)

# identify MLLs (Multilocus Lineages)
# alpha2 = 4 indicates that individuals with ≤ 4 allelic differences are grouped into the same MLL 
MLLlist <- MLL_generator(markers, alpha2 = 14)

###################### Genotypic diversity, richness and evenness indices calculation ################

# computes genotypic diversity and clonality indices
clonal_index(markers)

# computes genotypic diversity and clonality indices on multilocus lineage
clonal_index(markers, listMLL = MLLlist)

# N (total individuals): 367
# G (MLL) = 364 (three fewer than total, indicating three cases where individuals were grouped into the same MLL, meaning slight genetic similarity)
# R (Genotypic richness) = ~0.992
# Indicates high richness, but slightly lower than 1 due to the few grouped individuals (possible clones or contamination).
# H'' (Shannon-Wiener diversity) = 5.89403 (high value), indicating very high genotypic diversity.
# J' (Evenness index) ~0.99947, and D (Simpson's complement) ≈ 1.0
# Indicates very even and high diversity.
# Hill number (~364) approaches the actual number of unique MLLs, again showing high diversity.


###################################  Pareto index #######################################

Pareto_index(markers)

Pareto_index(markers, listMLL = MLLlist)

Pareto_index(markers, full = TRUE, graph = TRUE, legends = 2)


#################### Saving suspect individuals ###################################

# Generate MLLs using alpha2 = 5
MLLlist <- MLL_generator(markers, alpha2 = 5)

# Create an empty vector to store MLL assignments
MLL_vector <- rep(NA, length(indiv_ids))

# Properly fill MLL assignments based on MLLlist
for(i in seq_along(MLLlist)){
  individuals_in_MLL <- MLLlist[[i]]
  MLL_vector[individuals_in_MLL] <- i
}

# Check the correct assignment
MLL_df <- data.frame(ID = indiv_ids, MLL = MLL_vector)
head(MLL_df)

# Clearly identify MLL groups with multiple individuals (potential contamination)
library(dplyr)
suspected_MLLs <- MLL_df %>%
  group_by(MLL) %>%
  filter(n() > 1) %>% 
  arrange(MLL)

# View explicitly which individuals share MLL (potential contamination)
print(suspected_MLLs)

# Export to CSV clearly
write.csv(suspected_MLLs, 
          "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - SSR Populations/Analysis/01-Pre-traitements/01.7-checking_contamination_Rclone/suspected_contaminations_Regina.csv", 
          row.names = FALSE)


