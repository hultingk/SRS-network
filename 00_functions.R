# -------------------------------------- #
#### Functions script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# function to get data into correct format for network analysis
prepare_matrix <- function(df) {
  df_wide <- df %>% 
    pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    column_to_rownames("flower_species") #convert years to rownames
}

####### NULL DISTRIBUTION FUNCTIONS #######
# Null distribution function for nestedness 
net.null.nest = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'NODF'))
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

# Null distribution function for linkage density
net.null.density = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'linkage density'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}




####### Z SCORE FUNCTIONS #########
# Z-score function for comparing different networks
net.zscore = function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
} 

# Function that perform z-score calculation of nestedness using the observed and null networks
nest.zscore = function(nulltype){
  net.nest.zscore <- list() 
  for(i in 1:length(net.metrics.nest)){
    net.nest.zscore[[i]] = net.zscore(net.metrics.nest[[i]]['NODF'], 
                                      nulltype[[i]][ ,'NODF'])
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

# Function that perform z-score calculation of linkage density using the observed and null networks
density.zscore = function(nulltype){
  net.density.zscore <- list() 
  for(i in 1:length(net.metrics.density)){
    net.density.zscore[[i]] = net.zscore(net.metrics.density[[i]]['linkage density'], 
                                         nulltype[[i]][ ,'linkage density'])
  }
  names(net.density.zscore) <- webs.names
  return(net.density.zscore)
}