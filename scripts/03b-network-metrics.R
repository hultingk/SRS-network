# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# developed based on code from https://fukamilab.github.io/BIO202/09-B-networks.html 
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table)

set.seed(100)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))
load(file = file.path("data", "L2_nulls", "nulls.RData"))


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
# Calculate average links per species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species') 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.pol.links <- lapply(webs, grouplevel, index = 'mean number of links', level = "higher") 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.plant.links <- lapply(webs, grouplevel, index = 'mean number of links', level = "lower") 



# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest <- net.null.networklevel(nulls = net.nulls.vaz, fun = networklevel, metric = "NODF", level = "both", web.names = web.names)
vaz.h2 <- net.null.networklevel(nulls = net.nulls.vaz, fun = networklevel, metric = "H2", level = "both", web.names = web.names)
# vaz.links <- net.null.networklevel(nulls = net.nulls.vaz, fun = networklevel, metric = "links per species", level = "both", web.names = web.names)  # can't use -- NAN
vaz.pol.links <- net.null.networklevel(nulls = net.nulls.vaz, fun = grouplevel, metric = "mean number of links", level = "higher", web.names = web.names)
vaz.plant.links <- net.null.networklevel(nulls = net.nulls.vaz, fun = grouplevel, metric = "mean number of links", level = "lower", web.names = web.names)


# getting z score 
vaz.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = vaz.nest, metric = "NODF")
vaz.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = vaz.h2, metric = "H2")
vaz.pol.links.zscore <- zscore_metric(obsval = net.metrics.pol.links,
                                  nullval = vaz.pol.links, metric = "mean.number.of.links.HL") 
vaz.plant.links.zscore <- zscore_metric(obsval = net.metrics.plant.links,
                                      nullval = vaz.plant.links, metric = "mean.number.of.links.LL") 





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

links.per.sp <- as.data.frame(do.call('rbind', net.metrics.links)) 
links.per.sp <- links.per.sp %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(links.per.sp = `links per species`)

vaz.pol.links <- as.data.frame(do.call('rbind', vaz.pol.links.zscore)) 
vaz.pol.links <- vaz.pol.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links = mean.number.of.links.HL)

vaz.plant.links <- as.data.frame(do.call('rbind', vaz.plant.links.zscore)) 
vaz.plant.links <- vaz.plant.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links = mean.number.of.links.LL)




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
# Calculate average links per species for all plant-pollinator sites
net.metrics.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species') 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.pol.links_noApis <- lapply(webs_noApis, grouplevel, index = 'mean number of links', level = "higher") 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.plant.links_noApis <- lapply(webs_noApis, grouplevel, index = 'mean number of links', level = "lower") 


# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, fun = networklevel, metric = "NODF", level = "both")
vaz.h2_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, fun = networklevel, metric = "H2", level = "both")
# vaz.links_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, fun = networklevel, fun = networklevel, metric = "links per species", level = "both", web.names = web.names)  # can't use -- NAN
vaz.pol.links_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, fun = grouplevel, metric = "mean number of links", level = "higher", web.names = web.names)
vaz.plant.links_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, fun = grouplevel, metric = "mean number of links", level = "lower", web.names = web.names)


# getting z score 
vaz.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis,
                                 nullval = vaz.nest_noApis, metric = "NODF")
vaz.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis,
                               nullval = vaz.h2_noApis, metric = "H2")
vaz.pol.links.zscore_noApis <- zscore_metric(obsval = net.metrics.pol.links_noApis,
                                      nullval = vaz.pol.links_noApis, metric = "mean.number.of.links.HL") 
vaz.plant.links.zscore_noApis <- zscore_metric(obsval = net.metrics.plant.links_noApis,
                                        nullval = vaz.plant.links_noApis, metric = "mean.number.of.links.LL") 



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


links.per.sp_noApis <- as.data.frame(do.call('rbind', net.metrics.links_noApis)) 
links.per.sp_noApis <- links.per.sp_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(links.per.sp_noApis = `links per species`)

vaz.pol.links_noApis <- as.data.frame(do.call('rbind', vaz.pol.links.zscore_noApis)) 
vaz.pol.links_noApis <- vaz.pol.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links_noApis = mean.number.of.links.HL)

vaz.plant.links_noApis <- as.data.frame(do.call('rbind', vaz.plant.links.zscore_noApis)) 
vaz.plant.links_noApis <- vaz.plant.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links_noApis = mean.number.of.links.LL)




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
# Calculate average links per species for all plant-pollinator sites
net.metrics.links_noPoe <- lapply(webs_noPoe, networklevel, index = 'links per species') 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.pol.links_noPoe <- lapply(webs_noPoe, grouplevel, index = 'mean number of links', level = "higher") 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.plant.links_noPoe <- lapply(webs_noPoe, grouplevel, index = 'mean number of links', level = "lower") 





# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, fun = networklevel, metric = "NODF", level = "both")
vaz.h2_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, fun = networklevel, metric = "H2", level = "both")
# vaz.links_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, fun = networklevel, fun = networklevel, metric = "links per species", level = "both", web.names = web.names)  # can't use -- NAN
vaz.pol.links_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, fun = grouplevel, metric = "mean number of links", level = "higher", web.names = web.names)
vaz.plant.links_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, fun = grouplevel, metric = "mean number of links", level = "lower", web.names = web.names)


# getting z score 
vaz.nest.zscore_noPoe <- zscore_metric(obsval = net.metrics.nest_noPoe,
                                        nullval = vaz.nest_noPoe, metric = "NODF")
vaz.h2.zscore_noPoe <- zscore_metric(obsval = net.metrics.h2_noPoe,
                                      nullval = vaz.h2_noPoe, metric = "H2")
vaz.pol.links.zscore_noPoe <- zscore_metric(obsval = net.metrics.pol.links_noPoe,
                                             nullval = vaz.pol.links_noPoe, metric = "mean.number.of.links.HL") 
vaz.plant.links.zscore_noPoe <- zscore_metric(obsval = net.metrics.plant.links_noPoe,
                                               nullval = vaz.plant.links_noPoe, metric = "mean.number.of.links.LL") 




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

links.per.sp_noPoe <- as.data.frame(do.call('rbind', net.metrics.links_noPoe)) 
links.per.sp_noPoe <- links.per.sp_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(links.per.sp_noPoe = `links per species`)

vaz.pol.links_noPoe <- as.data.frame(do.call('rbind', vaz.pol.links.zscore_noPoe)) 
vaz.pol.links_noPoe <- vaz.pol.links_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links_noPoe = mean.number.of.links.HL)

vaz.plant.links_noPoe <- as.data.frame(do.call('rbind', vaz.plant.links.zscore_noPoe)) 
vaz.plant.links_noPoe <- vaz.plant.links_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links_noPoe = mean.number.of.links.LL)




#### network analysis: NO Poecilognathus sulphureus OR Apis mellifera ####
pollinator_split_noPoe_Apis <- pollinator %>%
  dplyr::filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split()

# getting each network into correct format
webs_noPoe_Apis <- pollinator_split_noPoe_Apis %>%
  lapply(prepare_matrix)

# adding patch names to each network
names(webs_noPoe_Apis) <- webs.names

# Calculate network connectance - unweighted
net.metrics.connect_noPoe_Apis <- lapply(webs_noPoe_Apis, networklevel, index = 'connectance') 
# Calculate network connectance - weighted
net.metrics.weightconnect_noPoe_Apis <- lapply(webs_noPoe_Apis, networklevel, index = 'weighted connectance') # NOTE this is linkage density divided by number of species in the network
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest_noPoe_Apis <- lapply(webs_noPoe_Apis, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noPoe_Apis <- lapply(webs_noPoe_Apis, networklevel, index = 'H2') 
# Calculate average links per species for all plant-pollinator sites
net.metrics.links_noPoe_Apis <- lapply(webs_noPoe_Apis, networklevel, index = 'links per species') 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.pol.links_noPoe_Apis <- lapply(webs_noPoe_Apis, grouplevel, index = 'mean number of links', level = "higher") 
# Calculate mean links per species for pollinators for all plant-pollinator sites
net.metrics.plant.links_noPoe_Apis <- lapply(webs_noPoe_Apis, grouplevel, index = 'mean number of links', level = "lower") 







# getting vaznull null for each metric - not for connectance or weighted connectance bc connectance is maintained in vaznull
vaz.nest_noPoe_Apis <- net.null.networklevel(nulls = net.nulls.vaz_noPoe_Apis, fun = networklevel, metric = "NODF", level = "both")
vaz.h2_noPoe_Apis <- net.null.networklevel(nulls = net.nulls.vaz_noPoe_Apis, fun = networklevel, metric = "H2", level = "both")
# vaz.links_noPoe_Apis <- net.null.networklevel(nulls = net.nulls.vaz_noPoe_Apis, fun = networklevel, metric = "links per species", level = "both", web.names = web.names)  # can't use -- NAN
vaz.pol.links_noPoe_Apis <- net.null.networklevel(nulls = net.nulls.vaz_noPoe_Apis, fun = grouplevel, metric = "mean number of links", level = "higher", web.names = web.names)
vaz.plant.links_noPoe_Apis <- net.null.networklevel(nulls = net.nulls.vaz_noPoe_Apis, fun = grouplevel, metric = "mean number of links", level = "lower", web.names = web.names)


# getting z score 
vaz.nest.zscore_noPoe_Apis <- zscore_metric(obsval = net.metrics.nest_noPoe_Apis,
                                       nullval = vaz.nest_noPoe_Apis, metric = "NODF")
vaz.h2.zscore_noPoe_Apis <- zscore_metric(obsval = net.metrics.h2_noPoe_Apis,
                                     nullval = vaz.h2_noPoe_Apis, metric = "H2")
vaz.pol.links.zscore_noPoe_Apis <- zscore_metric(obsval = net.metrics.pol.links_noPoe_Apis,
                                            nullval = vaz.pol.links_noPoe_Apis, metric = "mean.number.of.links.HL") 
vaz.plant.links.zscore_noPoe_Apis <- zscore_metric(obsval = net.metrics.plant.links_noPoe_Apis,
                                              nullval = vaz.plant.links_noPoe_Apis, metric = "mean.number.of.links.LL") 



# creating dataframes
net.connect_noPoe_Apis <- as.data.frame(do.call('rbind', net.metrics.connect_noPoe_Apis) )
net.connect_noPoe_Apis <- net.connect_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(connectance_noPoe_Apis = connectance)

net.weightconnect_noPoe_Apis <- as.data.frame(do.call('rbind', net.metrics.weightconnect_noPoe_Apis) )
net.weightconnect_noPoe_Apis <- net.weightconnect_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(weightconnectance_noPoe_Apis = `weighted connectance`)

vaz.nestedness_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.nest.zscore_noPoe_Apis) )
vaz.nestedness_noPoe_Apis <- vaz.nestedness_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF_noPoe_Apis = NODF)

vaz.h2_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.h2.zscore_noPoe_Apis) )
vaz.h2_noPoe_Apis <- vaz.h2_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.h2_noPoe_Apis = H2)

links.per.sp_noPoe_Apis <- as.data.frame(do.call('rbind', net.metrics.links_noPoe_Apis)) 
links.per.sp_noPoe_Apis <- links.per.sp_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(links.per.sp_noPoe_Apis = `links per species`)

vaz.pol.links_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.pol.links.zscore_noPoe_Apis)) 
vaz.pol.links_noPoe_Apis <- vaz.pol.links_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(pol.links_noPoe_Apis = mean.number.of.links.HL)

vaz.plant.links_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.plant.links.zscore_noPoe_Apis)) 
vaz.plant.links_noPoe_Apis <- vaz.plant.links_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(plant.links_noPoe_Apis = mean.number.of.links.LL)






#### putting all together ####
# first - all species
network_vaznull <- net.connect %>%
  left_join(net.weightconnect, by = c("block", "patch")) %>%
  left_join(vaz.nestedness, by = c("block", "patch")) %>%
  left_join(vaz.h2, by = c("block", "patch")) %>%
  left_join(links.per.sp, by = c("block", "patch")) %>%
  left_join(vaz.pol.links, by = c("block", "patch")) %>%
  left_join(vaz.plant.links, by = c("block", "patch")) 

# no apis
network_vaznull <- network_vaznull %>%
  left_join(net.connect_noApis, by = c("block", "patch")) %>%
  left_join(net.weightconnect_noApis, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noApis, by = c("block", "patch")) %>%
  left_join(vaz.h2_noApis, by = c("block", "patch")) %>%
  left_join(links.per.sp_noApis, by = c("block", "patch")) %>%
  left_join(vaz.pol.links_noApis, by = c("block", "patch")) %>%
  left_join(vaz.plant.links_noApis, by = c("block", "patch"))

# no poe
network_vaznull <- network_vaznull %>%
  left_join(net.connect_noPoe, by = c("block", "patch")) %>%
  left_join(net.weightconnect_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.h2_noPoe, by = c("block", "patch")) %>%
  left_join(links.per.sp_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.pol.links_noPoe, by = c("block", "patch")) %>%
  left_join(vaz.plant.links_noPoe, by = c("block", "patch")) 

# no poe or apis
network_vaznull <- network_vaznull %>%
  left_join(net.connect_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(net.weightconnect_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(vaz.h2_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(links.per.sp_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(vaz.pol.links_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(vaz.plant.links_noPoe_Apis, by = c("block", "patch"))



write.csv(network_vaznull, file = file.path("data", "L4_metrics", "network_vaznull.csv"))









