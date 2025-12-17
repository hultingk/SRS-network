# -------------------------------------- #
#### Diversity metrics analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# requires: 03-network-analysis.R
# -------------------------------------- #
# end result - summerized dataframe ready for analysis

# loading libraries
librarian::shelf(tidyverse, vegan, glmmTMB, iNEXT)

# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv")) # plant-pollinator interaction data
network_metrics <- read.csv(file = file.path("data", "network_metrics.csv")) # network results for joining at the end


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 

#### Pollinator abundance ####
# total pollinator visition
abundance <- pollinator %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance = n) %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  dplyr::select(!c("n"))

# pollinator visitation excuding Apis mellifera
abundance_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_noApis = n) %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  dplyr::select(!c("n"))


#### Flower diversity ####
# Hill numbers #
floral_inext <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each flower species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

floral_hill <- iNEXT(floral_inext, q = c(0,1,2), datatype = "abundance")
floral_est <- floral_hill$iNextEst[["coverage_based"]]
floral_est


floral_estimate <- estimateD(floral_inext, q = c(0,1,2), datatype = "abundance", level=NULL)
floral_estimate <- floral_estimate %>%
  dplyr::select(Assemblage, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  separate(Assemblage, into = c("block", "patch")) %>%
  rename(floral_hill_0 = `0`, floral_hill_1 = `1`, floral_hill_2 = `2`)

m1 <- glmmTMB(floral_hill_2 ~ patch + (1|block),
              data = floral_estimate)
summary(m1)


## flower diversity
floral_wider <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_diversity <- diversity(floral_wider, "shannon") # calculating diversity of flowers that are interacted with

# floral diversity with no Apis
floral_wider_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_diversity_noApis <- diversity(floral_wider_noApis, "shannon") # calculating diversity of flowers that are interacted with

# floral diversity no dominant pollinators. (Apis and Poecilognathus)
floral_wider_noDominant <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera", "Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_diversity_noDominant <- diversity(floral_wider_noDominant, "shannon") # calculating diversity of flowers that are interacted with




#### Pollinator diversity ####
## pollinator diversity
pollinator_wider <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_diversity <- diversity(pollinator_wider, "shannon") 

# pollinator diversity no apis
pollinator_wider_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_diversity_noApis <- diversity(pollinator_wider_noApis, "shannon") 


# pollinator diversity no dominant pollinators (Apis and Poecilognathus)
pollinator_wider_noDominant <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera", "Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_diversity_noDominant <- diversity(pollinator_wider_noDominant, "shannon") 

# Hill numbers #
pollinator_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera", "Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for pollinator species by patch matrix
  column_to_rownames("pollinator_species")

pollinator_hill <- iNEXT(pollinator_inext, q = c(0,1,2), datatype = "abundance")
pollinator_est <- pollinator_hill$iNextEst[["coverage_based"]]
pollinator_est


pollinator_estimate <- estimateD(pollinator_inext, q = c(0,1,2), datatype = "abundance", level=NULL)
pollinator_estimate <- pollinator_estimate %>%
  dplyr::select(Assemblage, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  separate(Assemblage, into = c("block", "patch")) %>%
  rename(pollinator_hill_0 = `0`, pollinator_hill_1 = `1`, pollinator_hill_2 = `2`)

m1 <- glmmTMB(pollinator_hill_2 ~ patch + (1|block),
              data = pollinator_estimate)
summary(m1)


# adding all together
diversity_metrics <- abundance %>%
  left_join(abundance_noApis, by = c("block", "patch"))
diversity_metrics$floral_diversity <- floral_diversity
diversity_metrics$floral_diversity_noApis <- floral_diversity_noApis
diversity_metrics$floral_diversity_noDominant <- floral_diversity_noDominant


diversity_metrics$pollinator_diversity <- pollinator_diversity
diversity_metrics$pollinator_diversity_noApis <- pollinator_diversity_noApis
diversity_metrics$pollinator_diversity_noDominant <- pollinator_diversity_noDominant



m1 <- glmmTMB(floral_diversity ~ patch + (1|block),
              data = diversity_metrics)
summary(m1)
m2 <- glmmTMB(floral_diversity_noApis ~ patch + (1|block),
              data = diversity_metrics)
summary(m2)
m2b <- glmmTMB(floral_diversity_noDominant ~ patch + (1|block),
              data = diversity_metrics)
summary(m2b)
m2c <- glmmTMB(floral_rare ~ patch + (1|block),
               data = diversity_metrics)
summary(m2c)
m2d <- glmmTMB(floral_rare_noApis ~ patch + (1|block),
               data = diversity_metrics)
summary(m2d)




m3 <- glmmTMB(pollinator_diversity ~ patch + (1|block),
              data = diversity_metrics)
summary(m3)
m4 <- glmmTMB(pollinator_diversity_noApis ~ patch + (1|block),
              data = diversity_metrics)
summary(m4)
m4b <- glmmTMB(pollinator_diversity_noDominant ~ patch + (1|block),
              data = diversity_metrics)
summary(m4b)
m4c <- glmmTMB(pollinator_rare ~ patch + (1|block),
               data = diversity_metrics)
summary(m4c)
m4d <- glmmTMB(pollinator_rare_noApis ~ patch + (1|block),
               data = diversity_metrics)
summary(m4d)


pollinator %>%
  filter(order == "Hymenoptera") %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  count(patch, pollinator_species)
