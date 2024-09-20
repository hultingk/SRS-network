# -------------------------------------- #
#### Data cleaning script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libaries
library(tidyverse)

# loading data 
pollinator <- read.csv(file = file.path("data", "2024-SRS-plant-pollinator.csv"))

# removing wasps, lost insects, unknowns, no observations
pollinator <- pollinator %>%
  filter(!pollinator_common %in% c("Wasp", "unknown butterfly", "unknown skipper", "No associated insect", "0", "WASP", "LOST", "NOT A POLLINATOR"))


# --------------------------- #
#### flower species cleaning ####
# --------------------------- #
# combining Eupatorium glaucescens and Eupatorium linearifolium as one species
pollinator$flower_species <- str_replace(pollinator$flower_species, "Eupatorium glaucescens", "Eupatorium linearifolium")

# checking flower species
pollinator %>%
  count(flower_species)

# --------------------------- #
#### pollinator species cleaning ####
# --------------------------- #
# grouping Erynnis species 
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis horatius", "Erynnis sp.")
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis zarucco", "Erynnis sp.")

# checking pollinator species
pollinator %>%
  count(pollinator_species)


# --------------------------- #
#### writing cleaned file ####
# --------------------------- #
write.csv(pollinator, file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


