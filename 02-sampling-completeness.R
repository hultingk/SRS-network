# -------------------------------------- #
#### Sampling completeness script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libaries
library(tidyverse)
library(iNEXT)
library(bipartite)

# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


sc_pollinator <- pollinator %>%
  filter(!pollinator_species %in% c("Milesia sp.", "Exoprosopa sp.", "", " ")) # excluding non-identified pollinators for now

# getting plant-pollinator data into correct format for sampling completeness estimation
sc_pollinator <- sc_pollinator %>%
  mutate(interaction = paste(pollinator_species, flower_species, sep = "-")) %>% # combining pollinator and flower into one column - creating interaction
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating unique ID for each patch
  count(unique_ID, interaction) %>% # counting # of each interaction per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) # pivoting into wider format for interaction by patch matrix

# excluding patch from matrix
sc_nosite_pollinator <- as.data.frame(sc_pollinator[, 2:15])

# estimating sampling completeness
out <- iNEXT(sc_nosite_pollinator, q=c(0), datatype="abundance", endpoint=NULL, se = TRUE) 

# mean sampling completeness across all patches
mean(out[["DataInfo"]][["SC"]])

# visualizing sampling completeness
ggiNEXT(out, type=1, facet.var="Order.q") # Sample-size-based R/E curve
ggiNEXT(out, type=2, facet.var="Order.q") # Sample completeness curve
ggiNEXT(out, type=3, facet.var="Order.q") # Coverage-based R/E curve
