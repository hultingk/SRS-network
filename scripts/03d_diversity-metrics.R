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


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 

#### Pollinator abundance ####
# total pollinator visitation
abundance <- pollinator %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# pollinator visitation excluding Apis mellifera
abundance_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_noApis = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# pollinator visitation excluding Poecilognathus sulphureus
abundance_noPoe <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_noPoe = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# pollinator visitation excluding Poecilognathus sulphureus and Apis mellifera
abundance_noApis_Poe <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_noApis_Poe = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# butterfly visitation 
abundance_lep <- pollinator %>%
  filter(order %in% c("Lepidoptera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_lep = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# bee visitation 
abundance_bee <- pollinator %>%
  filter(order %in% c("Hymenoptera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_bee = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))

# fly visitation 
abundance_fly <- pollinator %>%
  filter(order %in% c("Diptera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_fly = n) %>%
  separate(unique_ID, into = c("block", "patch"), remove = F) %>%
  dplyr::select(!c("n"))


#### Flower diversity ####
# Hill numbers #
floral_inext <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each flower species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

min(colSums(floral_inext))
floral_hill <- iNEXT(floral_inext, q = c(0,1,2), datatype = "abundance", size = c(160)) # estimating at abundance of smallest sample size (160 observations)
floral_est <- floral_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
floral_est <- floral_est %>% # subsetting estimates at abundance = 160
  filter(m == 160) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_0 = `0`, floral_1 = `1`, floral_2 = `2`)


# No Apis - Hill numbers #
floral_noApis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each floral species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

min(colSums(floral_noApis_inext))
floral_noApis_hill <- iNEXT(floral_noApis_inext, q = c(0,1,2), datatype = "abundance", size = c(141)) # estimating at abundance of smallest sample size (160 observations)
floral_noApis_est <- floral_noApis_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
floral_noApis_est <- floral_noApis_est %>% # subsetting estimates at abundance = 160
  filter(m == 141) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_noApis_0 = `0`, floral_noApis_1 = `1`, floral_noApis_2 = `2`)


# No Poecilognathus sulphureus - Hill numbers #
floral_noPoe_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each floral species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

min(colSums(floral_noPoe_inext))
floral_noPoe_hill <- iNEXT(floral_noPoe_inext, q = c(0,1,2), datatype = "abundance", size = c(106)) # estimating at abundance of smallest sample size (160 observations)
floral_noPoe_est <- floral_noPoe_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
floral_noPoe_est <- floral_noPoe_est %>% # subsetting estimates at abundance = 160
  filter(m == 106) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>% 
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_noPoe_0 = `0`, floral_noPoe_1 = `1`, floral_noPoe_2 = `2`)


# No Poecilognathus sulphureus or Apis mellifera - Hill numbers #
floral_noPoe_Apis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID, flower_species) %>% # counting # of each floral species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("flower_species")

min(colSums(floral_noPoe_Apis_inext))
floral_noPoe_Apis_hill <- iNEXT(floral_noPoe_Apis_inext, q = c(0,1,2), datatype = "abundance", size = c(104)) # estimating at abundance of smallest sample size (160 observations)
floral_noPoe_Apis_est <- floral_noPoe_Apis_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
floral_noPoe_Apis_est <- floral_noPoe_Apis_est %>% # subsetting estimates at abundance = 160
  filter(m == 104) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, floral_noPoe_Apis_0 = `0`, floral_noPoe_Apis_1 = `1`, floral_noPoe_Apis_2 = `2`)



#### Pollinator diversity ####
# All species - Hill numbers #
pollinator_inext <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(pollinator_inext))
pollinator_hill <- iNEXT(pollinator_inext, q = c(0,1,2), datatype = "abundance", size = c(160)) # estimating at abundance of smallest sample size (160 observations)
pollinator_est <- pollinator_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
pollinator_est <- pollinator_est %>% # subsetting estimates at abundance = 160
  filter(m == 160) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_0 = `0`, pollinator_1 = `1`, pollinator_2 = `2`)


# No Apis - Hill numbers #
pollinator_noApis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(pollinator_noApis_inext))
pollinator_noApis_hill <- iNEXT(pollinator_noApis_inext, q = c(0,1,2), datatype = "abundance", size = c(141)) # estimating at abundance of smallest sample size (160 observations)
pollinator_noApis_est <- pollinator_noApis_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
pollinator_noApis_est <- pollinator_noApis_est %>% # subsetting estimates at abundance = 160
  filter(m == 141) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_noApis_0 = `0`, pollinator_noApis_1 = `1`, pollinator_noApis_2 = `2`)



# No Poecilognathus sulphureus - Hill numbers #
pollinator_noPoe_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(pollinator_noPoe_inext))
pollinator_noPoe_hill <- iNEXT(pollinator_noPoe_inext, q = c(0,1,2), datatype = "abundance", size = c(106)) # estimating at abundance of smallest sample size (160 observations)
pollinator_noPoe_est <- pollinator_noPoe_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
pollinator_noPoe_est <- pollinator_noPoe_est %>% # subsetting estimates at abundance = 160
  filter(m == 106) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_noPoe_0 = `0`, pollinator_noPoe_1 = `1`, pollinator_noPoe_2 = `2`)


# No Poecilognathus sulphureus or Apis mellifera - Hill numbers #
pollinator_noPoe_Apis_inext <- pollinator %>%
  filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(pollinator_noPoe_Apis_inext))
pollinator_noPoe_Apis_hill <- iNEXT(pollinator_noPoe_Apis_inext, q = c(0,1,2), datatype = "abundance", size = c(104)) # estimating at abundance of smallest sample size (160 observations)
pollinator_noPoe_Apis_est <- pollinator_noPoe_Apis_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
pollinator_noPoe_Apis_est <- pollinator_noPoe_Apis_est %>% # subsetting estimates at abundance = 160
  filter(m == 104) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, pollinator_noPoe_Apis_0 = `0`, pollinator_noPoe_Apis_1 = `1`, pollinator_noPoe_Apis_2 = `2`)



# Leps- Hill numbers #
lep_inext <- pollinator %>%
  filter(order %in% c("Lepidoptera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(lep_inext))
lep_hill <- iNEXT(lep_inext, q = c(0,1,2), datatype = "abundance", size = c(32)) # estimating at abundance of smallest sample size (160 observations)
lep_est <- lep_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
lep_est <- lep_est %>% # subsetting estimates at abundance = 160
  filter(m == 32) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, lep_0 = `0`, lep_1 = `1`, lep_2 = `2`)

# Bees - Hill numbers #
bee_inext <- pollinator %>%
  filter(order %in% c("Hymenoptera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(bee_inext))
bee_hill <- iNEXT(bee_inext, q = c(0,1,2), datatype = "abundance", size = c(46)) # estimating at abundance of smallest sample size (160 observations)
bee_est <- bee_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
bee_est <- bee_est %>% # subsetting estimates at abundance = 160
  filter(m == 46) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, bee_0 = `0`, bee_1 = `1`, bee_2 = `2`)


# fly - Hill numbers #
fly_inext <- pollinator %>%
  filter(order %in% c("Diptera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>% # counting # of each pollinator species per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) %>% # pivoting into wider format for flower species by patch matrix
  column_to_rownames("pollinator_species")

min(colSums(fly_inext))
fly_hill <- iNEXT(fly_inext, q = c(0,1,2), datatype = "abundance", size = c(20)) # estimating at abundance of smallest sample size (160 observations)
fly_est <- fly_hill$iNextEst[["coverage_based"]] # subsetting coverage based estimates
fly_est <- fly_est %>% # subsetting estimates at abundance = 160
  filter(m == 20) %>%
  separate(Assemblage, into = c("block", "patch"), remove = F) %>%
  mutate(Order.q = as.factor(Order.q)) %>%
  dplyr::select(Assemblage, block, patch, Order.q, qD) %>%
  pivot_wider(names_from = Order.q, values_from = qD) %>%
  dplyr::rename(unique_ID = Assemblage, fly_0 = `0`, fly_1 = `1`, fly_2 = `2`)





# putting all together
diversity_metrics <- abundance %>%
  left_join(abundance_noApis, by = c("unique_ID", "block", "patch")) %>%
  left_join(abundance_noPoe, by = c("unique_ID", "block", "patch")) %>%
  left_join(abundance_noApis_Poe, by = c("unique_ID", "block", "patch")) %>%
  left_join(abundance_bee, by = c("unique_ID", "block", "patch")) %>%
  left_join(abundance_fly, by = c("unique_ID", "block", "patch")) %>%
  left_join(abundance_lep, by = c("unique_ID", "block", "patch")) %>%
  left_join(floral_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(floral_noApis_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(floral_noPoe_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(floral_noPoe_Apis_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_noApis_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_noPoe_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(pollinator_noPoe_Apis_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(lep_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(bee_est, by = c("unique_ID", "block", "patch")) %>%
  left_join(fly_est, by = c("unique_ID", "block", "patch")) #%>%
  #dplyr::select(!ends_with("0"))
  
write.csv(diversity_metrics, file = file.path("data", "diversity_metrics.csv"))

