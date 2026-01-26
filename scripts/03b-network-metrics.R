# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# developed based on code from https://fukamilab.github.io/BIO202/09-B-networks.html 
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, parallel, data.table)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))
load(file = file.path("data", "nulls.RData"))


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 


#### network analysis: WITH ALL SPECIES ####
pollinator_split <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 

# getting each network into correct format
webs <- pollinator_split %>%
  lapply(prepare_matrix)

# adding patch names to each network
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")
names(webs) <- webs.names

# Calculate network connectance - unweighted
net.metrics.connect <- lapply(webs, networklevel, index = 'connectance') 
# Calculate network connectance - weighted
net.metrics.weightconnect <- lapply(webs, networklevel, index = 'weighted connectance') # NOTE this is linkage density divided by number of species in the network
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity <- lapply(webs, networklevel, index = 'Shannon diversity') 
# Calculate network interaction evenness for all plant-pollinator sites
net.metrics.even <- lapply(webs, networklevel, index = 'interaction evenness') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.pol.links <- lapply(webs, networklevel, index = 'links per species', level = "higher") 
# Calculate links per plant species for all plant-pollinator sites
net.metrics.plant.links <- lapply(webs, networklevel, index = 'links per species', level = "lower") 


# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest <- net.null.networklevel(nulls = net.nulls.vaz, metric = "NODF", level = "both", web.names = web.names)
vaz.h2 <- net.null.networklevel(nulls = net.nulls.vaz, metric = "H2", level = "both", web.names = web.names)
vaz.diversity <- net.null.networklevel(nulls = net.nulls.vaz, metric = "Shannon diversity", level = "both", web.names = web.names)
vaz.even <- net.null.networklevel(nulls = net.nulls.vaz, metric = "interaction evenness", level = "both", web.names = web.names)
# vaz.pol.links <- net.null.networklevel(nulls = net.nulls.vaz, metric = "links per species", level = "higher", web.names = web.names)  # can't use -- NAN
# vaz.plant.links <- net.null.networklevel(nulls = net.nulls.vaz, metric = "links per species", level = "lower", web.names = web.names) # can't use -- NAN


# getting z score 
vaz.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = vaz.nest, metric = "NODF")
vaz.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = vaz.h2, metric = "H2")
vaz.diversity.zscore <- zscore_metric(obsval = net.metrics.diversity,
                                      nullval = vaz.diversity, metric = "Shannon diversity")
vaz.even.zscore <- zscore_metric(obsval = net.metrics.even,
                                      nullval = vaz.even, metric = "interaction evenness")
# vaz.pol.links.zscore <- zscore_metric(obsval = net.metrics.pol.links,
#                                   nullval = vaz.pol.links, metric = "links per species") # can't use -- NAN
# vaz.plant.links.zscore <- zscore_metric(obsval = net.metrics.plant.links,
#                                       nullval = vaz.plant.links, metric = "links per species") # can't use -- NAN





# creating dataframes
net.connect <- as.data.frame(do.call('rbind', net.metrics.connect) )
net.connect <- net.connect %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(connectance = connectance)

net.weightconnect <- as.data.frame(do.call('rbind', net.metrics.weightconnect) )
net.weightconnect <- net.weightconnect %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(weightconnectance = `weighted connectance`)

vaz.nestedness <- as.data.frame(do.call('rbind', vaz.nest.zscore) )
vaz.nestedness <- vaz.nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF = NODF)

vaz.h2 <- as.data.frame(do.call('rbind', vaz.h2.zscore) )
vaz.h2 <- vaz.h2 %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.h2 = H2)

vaz.shannon <- as.data.frame(do.call('rbind', vaz.diversity.zscore) )
vaz.shannon <- vaz.shannon %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.shannon = `Shannon diversity`)

vaz.evenness <- as.data.frame(do.call('rbind', vaz.even.zscore) )
vaz.evenness <- vaz.evenness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.evenness = `interaction evenness`)

pol.links <- as.data.frame(do.call('rbind', net.metrics.pol.links)) 
pol.links <- pol.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links = `links per species`)

plant.links <- as.data.frame(do.call('rbind', net.metrics.plant.links)) 
plant.links <- plant.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links = `links per species`)




#### network analysis: NO APIS ####
pollinator_split_noApis <- pollinator %>%
  dplyr::filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split()

# getting each network into correct format
webs_noApis <- pollinator_split_noApis %>%
  lapply(prepare_matrix)

# adding patch names to each network
names(webs_noApis) <- webs.names

# Calculate network connectance - unweighted
net.metrics.connect_noApis <- lapply(webs_noApis, networklevel, index = 'connectance') 
# Calculate network connectance - weighted
net.metrics.weightconnect_noApis <- lapply(webs_noApis, networklevel, index = 'weighted connectance') # NOTE this is linkage density divided by number of species in the network
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest_noApis <- lapply(webs_noApis, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noApis <- lapply(webs_noApis, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity_noApis <- lapply(webs_noApis, networklevel, index = 'Shannon diversity') 
# Calculate network interaction evenness for all plant-pollinator sites
net.metrics.even_noApis <- lapply(webs_noApis, networklevel, index = 'interaction evenness') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.pol.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species', level = "higher") 
# Calculate links per plant species for all plant-pollinator sites
net.metrics.plant.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species', level = "lower") 


# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "NODF", level = "both")
vaz.h2_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "H2", level = "both")
vaz.diversity_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "Shannon diversity", level = "both")
vaz.even_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "interaction evenness", level = "both")
# vaz.pol.links_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "links per species", level = "higher") # can't use -- NAN
# vaz.plant.links_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "links per species", level = "lower") # can't use -- NAN


# getting z score 
vaz.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis,
                                 nullval = vaz.nest_noApis, metric = "NODF")
vaz.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis,
                               nullval = vaz.h2_noApis, metric = "H2")
vaz.diversity.zscore_noApis <- zscore_metric(obsval = net.metrics.diversity_noApis,
                                      nullval = vaz.diversity_noApis, metric = "Shannon diversity")
vaz.even.zscore_noApis <- zscore_metric(obsval = net.metrics.even_noApis,
                                 nullval = vaz.even_noApis, metric = "interaction evenness")
# vaz.pol.links.zscore_noApis <- zscore_metric(obsval = net.metrics.pol.links_noApis,
#                                       nullval = vaz.pol.links_noApis, metric = "links per species") # can't use -- NAN
# vaz.plant.links.zscore_noApis <- zscore_metric(obsval = net.metrics.plant.links_noApis,
#                                         nullval = vaz.plant.links_noApis, metric = "links per species") # can't use -- NAN




# creating dataframes
net.connect_noApis <- as.data.frame(do.call('rbind', net.metrics.connect_noApis) )
net.connect_noApis <- net.connect_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(connectance_noApis = connectance)

net.weightconnect_noApis <- as.data.frame(do.call('rbind', net.metrics.weightconnect_noApis) )
net.weightconnect_noApis <- net.weightconnect_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(weightconnectance_noApis = `weighted connectance`)

vaz.nestedness_noApis <- as.data.frame(do.call('rbind', vaz.nest.zscore_noApis) )
vaz.nestedness_noApis <- vaz.nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF_noApis = NODF)

vaz.h2_noApis <- as.data.frame(do.call('rbind', vaz.h2.zscore_noApis) )
vaz.h2_noApis <- vaz.h2_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.h2_noApis = H2)

vaz.shannon_noApis <- as.data.frame(do.call('rbind', vaz.diversity.zscore_noApis) )
vaz.shannon_noApis <- vaz.shannon_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.shannon_noApis = `Shannon diversity`)

vaz.evenness_noApis <- as.data.frame(do.call('rbind', vaz.even.zscore_noApis) )
vaz.evenness_noApis <- vaz.evenness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.evenness_noApis = `interaction evenness`)

pol.links_noApis <- as.data.frame(do.call('rbind', net.metrics.pol.links_noApis)) 
pol.links_noApis <- pol.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links_noApis = `links per species`)

plant.links_noApis <- as.data.frame(do.call('rbind', net.metrics.plant.links_noApis)) 
plant.links_noApis <- plant.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links_noApis = `links per species`)




#### network analysis: NO Poecilognathus sulphureus ####
pollinator_split_noPoe <- pollinator %>%
  dplyr::filter(!pollinator_species %in% c("Poecilognathus sulphureus")) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split()

# getting each network into correct format
webs_noPoe <- pollinator_split_noPoe %>%
  lapply(prepare_matrix)

# adding patch names to each network
names(webs_noPoe) <- webs.names

# Calculate network connectance - unweighted
net.metrics.connect_noPoe <- lapply(webs_noPoe, networklevel, index = 'connectance') 
# Calculate network connectance - weighted
net.metrics.weightconnect_noPoe <- lapply(webs_noPoe, networklevel, index = 'weighted connectance') # NOTE this is linkage density divided by number of species in the network
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest_noPoe <- lapply(webs_noPoe, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noPoe <- lapply(webs_noPoe, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity_noPoe <- lapply(webs_noPoe, networklevel, index = 'Shannon diversity') 
# Calculate network interaction evenness for all plant-pollinator sites
net.metrics.even_noPoe <- lapply(webs_noPoe, networklevel, index = 'interaction evenness') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.pol.links_noPoe <- lapply(webs_noPoe, networklevel, index = 'links per species', level = "higher") 
# Calculate links per plant species for all plant-pollinator sites
net.metrics.plant.links_noPoe <- lapply(webs_noPoe, networklevel, index = 'links per species', level = "lower") 

# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "NODF", level = "both")
vaz.h2_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "H2", level = "both")
vaz.diversity_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "Shannon diversity", level = "both")
vaz.even_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "interaction evenness", level = "both")
# vaz.pol.links_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "links per species", level = "higher") # can't use -- NAN
# vaz.plant.links_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "links per species", level = "lower") # can't use -- NAN


# getting z score 
vaz.nest.zscore_noPoe <- zscore_metric(obsval = net.metrics.nest_noPoe,
                                        nullval = vaz.nest_noPoe, metric = "NODF")
vaz.h2.zscore_noPoe <- zscore_metric(obsval = net.metrics.h2_noPoe,
                                      nullval = vaz.h2_noPoe, metric = "H2")
vaz.diversity.zscore_noPoe <- zscore_metric(obsval = net.metrics.diversity_noPoe,
                                             nullval = vaz.diversity_noPoe, metric = "Shannon diversity")
vaz.even.zscore_noPoe <- zscore_metric(obsval = net.metrics.even_noPoe,
                                        nullval = vaz.even_noPoe, metric = "interaction evenness")
# vaz.pol.links.zscore_noPoe <- zscore_metric(obsval = net.metrics.pol.links_noPoe,
#                                              nullval = vaz.pol.links_noPoe, metric = "links per species") # can't use -- NAN
# vaz.plant.links.zscore_noPoe <- zscore_metric(obsval = net.metrics.plant.links_noPoe,
#                                                nullval = vaz.plant.links_noPoe, metric = "links per species") # can't use -- NAN




# creating dataframes
net.connect_noPoe <- as.data.frame(do.call('rbind', net.metrics.connect_noPoe) )
net.connect_noPoe <- net.connect_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(connectance_noPoe = connectance)

net.weightconnect_noPoe <- as.data.frame(do.call('rbind', net.metrics.weightconnect_noPoe) )
net.weightconnect_noPoe <- net.weightconnect_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(weightconnectance_noPoe = `weighted connectance`)

vaz.nestedness_noPoe <- as.data.frame(do.call('rbind', vaz.nest.zscore_noPoe) )
vaz.nestedness_noPoe <- vaz.nestedness_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF_noPoe = NODF)

vaz.h2_noPoe <- as.data.frame(do.call('rbind', vaz.h2.zscore_noPoe) )
vaz.h2_noPoe <- vaz.h2_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.h2_noPoe = H2)

vaz.shannon_noPoe <- as.data.frame(do.call('rbind', vaz.diversity.zscore_noPoe) )
vaz.shannon_noPoe <- vaz.shannon_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.shannon_noPoe = `Shannon diversity`)

vaz.evenness_noPoe <- as.data.frame(do.call('rbind', vaz.even.zscore_noPoe) )
vaz.evenness_noPoe <- vaz.evenness_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.evenness_noPoe = `interaction evenness`)

pol.links_noPoe <- as.data.frame(do.call('rbind', net.metrics.pol.links_noPoe)) 
pol.links_noPoe <- pol.links_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links_noPoe = `links per species`)

plant.links_noPoe <- as.data.frame(do.call('rbind', net.metrics.plant.links_noPoe)) 
plant.links_noPoe <- plant.links_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links_noPoe = `links per species`)




# putting all together
network_vaznull <- net.connect %>%
  left_join(net.weightconnect, by = c("block", "patch")) %>%
  left_join(vaz.nestedness, by = c("block", "patch")) %>%
  left_join(vaz.h2, by = c("block", "patch")) %>%
  left_join(vaz.shannon, by = c("block", "patch")) %>%
  left_join(vaz.evenness, by = c("block", "patch")) %>%
  left_join(pol.links, by = c("block", "patch")) %>%
  left_join(plant.links, by = c("block", "patch")) %>%
  left_join(net.connect_noApis, by = c("block", "patch")) %>%
  left_join(net.weightconnect_noApis, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noApis, by = c("block", "patch")) %>%
  left_join(vaz.h2_noApis, by = c("block", "patch")) %>%
  left_join(vaz.shannon_noApis, by = c("block", "patch")) %>%
  left_join(vaz.evenness_noApis, by = c("block", "patch")) %>%
  left_join(pol.links_noApis, by = c("block", "patch")) %>%
  left_join(plant.links_noApis, by = c("block", "patch")) %>%
  left_join(net.connect_noPoe, by = c("block", "patch")) %>%
  left_join(net.weightconnect_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.h2_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.shannon_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.evenness_noPoe, by = c("block", "patch")) %>%
  left_join(pol.links_noPoe, by = c("block", "patch")) %>%
  left_join(plant.links_noPoe, by = c("block", "patch"))

write.csv(network_vaznull, file = file.path("data", "network_vaznull.csv"))









