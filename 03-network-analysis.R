# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libaries
library(tidyverse)
library(bipartite)
library(vegan)



# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))

# removing non-IDed species, creating unique ID
pollinator <- pollinator %>%
  #filter(!pollinator_analysis %in% c("Apis mellifera")) %>%
  filter(!pollinator_species == "") %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "pollinator_analysis", "flower_species")) 



#### network analysis ####
pollinator_split <- pollinator %>%
  count(unique_ID, pollinator_analysis, flower_species) %>%
  group_by(unique_ID) %>%
  group_split() 
  
# function to get data into correct format
prepare_matrix <- function(df) {
  df_wide <- df %>% 
    pivot_wider(names_from = pollinator_analysis, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    column_to_rownames("flower_species") #convert years to rownames
}

# getting each network into correct format
webs <- pollinator_split %>%
  lapply(prepare_matrix)
  
# adding patch names to each network
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")
names(webs) <- webs.names


# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'nestedness') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity <- lapply(webs, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species', level = "higher") 



# Make null models for all sites using the vaznull null
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) 


# Null distribution function for nestedness 
net.null.nest = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'nestedness'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}

# Null distribution function for specialization 
net.null.h2 = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'H2'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}

# Null distribution function for interaction diversity 
net.null.diversity = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'Shannon diversity'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}

# Null distribution function for links per species 
net.null.links = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'links per species', level = "higher"))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}


vaz.nest <- net.null.nest(net.nulls.vaz)
vaz.h2 <- net.null.h2(net.nulls.vaz)
vaz.diversity <- net.null.diversity(net.nulls.vaz)
vaz.links <- net.null.links(net.nulls.vaz)

# Z-score function for comparing different networks
net.zscore = function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
} 

# Function that perform z-score calculation of nestedness using the observed and null networks
nest.zscore = function(nulltype){
  net.nest.zscore <- list() 
  for(i in 1:length(net.metrics.nest)){
    net.nest.zscore[[i]] = net.zscore(net.metrics.nest[[i]]['nestedness'], 
                                      nulltype[[i]][ ,'nestedness'])
  }
  names(net.nest.zscore) <- webs.names
  return(net.nest.zscore)
}

# Function that perform z-score calculation of specialization using the observed and null networks
h2.zscore = function(nulltype){
  net.h2.zscore <- list() 
  for(i in 1:length(net.metrics.h2)){
    net.h2.zscore[[i]] = net.zscore(net.metrics.h2[[i]]['H2'], 
                                       nulltype[[i]][ ,'H2'])
  }
  names(net.h2.zscore) <- webs.names
  return(net.h2.zscore)
}

# Function that perform z-score calculation of interaction diversity using the observed and null networks
diversity.zscore = function(nulltype){
  net.diversity.zscore <- list() 
  for(i in 1:length(net.metrics.diversity)){
    net.diversity.zscore[[i]] = net.zscore(net.metrics.diversity[[i]]['Shannon diversity'], 
                                    nulltype[[i]][ ,'Shannon diversity'])
  }
  names(net.diversity.zscore) <- webs.names
  return(net.diversity.zscore)
}

# Function that perform z-score calculation of interaction diversity using the observed and null networks
links.zscore = function(nulltype){
  net.links.zscore <- list() 
  for(i in 1:length(net.metrics.links)){
    net.links.zscore[[i]] = net.zscore(net.metrics.links[[i]]['links per species'], 
                                           nulltype[[i]][ ,'links per species'])
  }
  names(net.links.zscore) <- webs.names
  return(net.links.zscore)
}


vaz.nest.zscore <- nest.zscore(vaz.nest)
vaz.h2.zscore <- h2.zscore(vaz.h2)
vaz.diversity.zscore <- diversity.zscore(vaz.diversity)
vaz.links.zscore <- links.zscore(vaz.links) # CAN'T USE - NO SD FOR NULL MODEL, USE REAL VALUES IN ANALYSIS

# creating dataframes
nestedness <- as.data.frame(do.call('rbind', vaz.nest.zscore) )
nestedness <- nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into columns
  
h2 <- as.data.frame(do.call('rbind', vaz.h2.zscore) )
h2 <- h2 %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

shannon <- as.data.frame(do.call('rbind', vaz.diversity.zscore) )
shannon <- shannon %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

links <- as.data.frame(do.call('rbind', net.metrics.links))
links <- links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column





#### Abundance ####
abundance <- pollinator %>%
  count(unique_ID) %>%
  mutate(abundance = n) %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  dplyr::select(!c("n"))

#### Flower and pollinator species diversity ####
## flower diversity
floral_wider <- pollinator %>%
  count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_diversity <- diversity(floral_wider, "shannon") # calculating diversity of flowers that are interacted with

floral_wider_noApis <- pollinator %>%
  filter(!pollinator_analysis %in% c("Apis mellifera")) %>%
  count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_div_noApis <- diversity(floral_wider_noApis, "shannon") # calculating diversity of flowers that are interacted with


## pollinator diversity
pollinator_wider <- pollinator %>%
  count(unique_ID, pollinator_analysis) %>%
  pivot_wider(names_from = pollinator_analysis, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_diversity <- diversity(pollinator_wider, "shannon") 

pollinator_wider_noApis <- pollinator %>%
  filter(!pollinator_analysis %in% c("Apis mellifera")) %>%
  count(unique_ID, pollinator_analysis) %>%
  pivot_wider(names_from = pollinator_analysis, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_div_noApis <- diversity(pollinator_wider, "shannon") 


#### exporting csv ####
# adding all together
network_metrics <- nestedness %>%
  left_join(h2, by = c("block", "patch")) %>%
  left_join(shannon, by = c("block", "patch")) %>%
  left_join(links, by = c("block", "patch")) %>%
  left_join(abundance, by = c("block", "patch"))
network_metrics$floral_diversity <- floral_diversity
network_metrics$floral_div_noApis <- floral_div_noApis
network_metrics$pollinator_diversity <- pollinator_diversity
network_metrics$pollinator_div_noApis <- pollinator_div_noApis


write.csv(network_metrics, file = file.path("data", "network_metrics.csv"))










