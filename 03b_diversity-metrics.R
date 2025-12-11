# -------------------------------------- #
#### Diversity metrics analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# requires: 03-network-analysis.R
# -------------------------------------- #
# end result - summerized dataframe ready for analysis

# loading libraries
librarian::shelf(tidyverse, vegan)

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


#### Flower and pollinator species diversity ####
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
floral_div_noApis <- diversity(floral_wider_noApis, "shannon") # calculating diversity of flowers that are interacted with


# floral richness - rarefied
spAbund <- rowSums(floral_wider) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare <- rarefy(floral_wider, 160) # now use function rarefy

# floral richness - rarefied no apis
spAbund_noApis <- rowSums(floral_wider_noApis) # calculating minimum # of observation 
min(spAbund_noApis) # 7 is the fewest interactions observed per patch
sRare_noApis <- rarefy(floral_wider_noApis, 141) # now use function rarefy


# floral richness - not rarefied
floral_wider[floral_wider>0] <- 1 
fl.rich <- rowSums(floral_wider)



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
pollinator_div_noApis <- diversity(pollinator_wider_noApis, "shannon") 



# rarified pollinator richness
spAbund <- rowSums(pollinator_wider) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare_pollinator <- rarefy(pollinator_wider, 160) # now use function rarefy

# rarified pollinator richness no apis
spAbund <- rowSums(pollinator_wider_noApis) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare_pollinator_noApis <- rarefy(pollinator_wider_noApis, 141) # now use function rarefy
