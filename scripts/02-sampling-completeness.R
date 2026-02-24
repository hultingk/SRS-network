# -------------------------------------- #
#### Sampling completeness script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, iNEXT, bipartite)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))

# excluding non-identified pollinators for now
# sc_pollinator <- pollinator %>%
#   filter(!pollinator_species %in% c(" ", "")) 

# getting plant-pollinator data into correct format for sampling completeness estimation
sc_pollinator <- pollinator %>%
  mutate(interaction = paste(pollinator_species, flower_species, sep = "-")) %>% # combining pollinator and flower into one column - creating interaction
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating unique ID for each patch
  dplyr::count(unique_ID, interaction) %>% # counting # of each interaction per patch
  pivot_wider(names_from = unique_ID, values_from = n, values_fill = 0) # pivoting into wider format for interaction by patch matrix

# excluding patch from interaction matrix
sc_nosite_pollinator <- as.data.frame(sc_pollinator[, 2:15])

# estimating sampling completeness
out <- iNEXT(sc_nosite_pollinator, q=c(1), datatype="abundance", endpoint=NULL, se = TRUE) 

# mean sampling completeness across all patches
mean(out[["DataInfo"]][["SC"]])

# visualizing sampling completeness
ggiNEXT(out, type=1, facet.var="Order.q") # Sample-size-based R/E curve
ggiNEXT(out, type=2, facet.var="Order.q") # Sample completeness curve
ggiNEXT(out, type=3, facet.var="Order.q") # Coverage-based R/E curve


sampling_plot <- ggiNEXT(out, type=2, facet.var="Order.q") # Sample completeness curve
sampling_plot <- sampling_plot + theme_classic()
sampling_plot
