# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# developed based on code from https://fukamilab.github.io/BIO202/09-B-networks.html 
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


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

#web.plot <- plotweb(webs[["53N.B"]], text.rot=90)
#pdf(file = file.path("plots", "web.plot.pdf"), width = 8, height = 6)
#plotweb(webs[["53N.B"]], text.rot=90, y.lim = c(-1,2.5))
#dev.off()

# Calculate network connectance - unweighted
net.metrics.connect <- lapply(webs, networklevel, index = 'connectance') 
# Calculate network connectance - weighted
net.metrics.weightconnect <- lapply(webs, networklevel, index = 'weighted connectance') # NOTE this is linkage density divided by number of species in the network
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'NODF') 
# Calculate network nestedness - weighted
net.metrics.weightnest <- lapply(webs, networklevel, index = 'weighted NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity <- lapply(webs, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density <- lapply(webs, networklevel, index = "linkage density")


# Make null models for all sites 
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) # using the vaznull null - maintains connectance 
net.nulls.r2d <- lapply(webs, nullmodel, method = "r2dtable", N = 500) # using the r2dtable null - maintains row and column sums


# getting r2dtable null for each metric
r2d.connect <- net.null.connect(net.nulls.r2d) # can only use r2dtable null model for connectance -- vaznull maintains connectance
r2d.weightconnect <- net.null.weightconnect(net.nulls.r2d) # can only use r2dtable null model for connectance -- vaznull maintains connectance
r2d.nest <- net.null.nest(net.nulls.r2d)
r2d.weightnest <- net.null.weightnest(net.nulls.r2d)
r2d.h2 <- net.null.h2(net.nulls.r2d)
r2d.diversity <- net.null.diversity(net.nulls.r2d)
r2d.links <- net.null.links(net.nulls.r2d)
r2d.density <- net.null.density(net.nulls.r2d)

# getting vaznull null for each metric
vaz.nest <- net.null.nest(net.nulls.vaz)
vaz.weightnest <- net.null.weightnest(net.nulls.vaz)
vaz.h2 <- net.null.h2(net.nulls.vaz)
vaz.diversity <- net.null.diversity(net.nulls.vaz)
vaz.links <- net.null.links(net.nulls.vaz)
vaz.density <- net.null.density(net.nulls.vaz)


# getting z score - r2dtable
r2d.connect.zscore <- zscore_metric(obsval = net.metrics.connect, 
                                    nullval = r2d.connect, metric = "connectance")
r2d.weightconnect.zscore <- zscore_metric(obsval = net.metrics.weightconnect,
                                          nullval = r2d.weightconnect, metric = "weighted connectance")
r2d.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = r2d.nest, metric = "NODF")
r2d.weightnest.zscore <- zscore_metric(obsval = net.metrics.weightnest,
                                       nullval = r2d.weightnest, metric = "weighted NODF")
r2d.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = r2d.h2, metric = "H2")
r2d.diversity.zscore <- zscore_metric(obsval = net.metrics.diversity,
                                      nullval = r2d.diversity, metric = "Shannon diversity")
r2d.links.zscore <- zscore_metric(obsval = net.metrics.links,
                                  nullval = r2d.links, metric = "links per species")
r2d.density.zscore <- zscore_metric(obsval = net.metrics.density,
                                    nullval = r2d.density, metric = "linkage density")

# getting z score - vaznull
vaz.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = vaz.nest, metric = "NODF")
vaz.weightnest.zscore <- zscore_metric(obsval = net.metrics.weightnest,
                                       nullval = vaz.weightnest, metric = "weighted NODF")
vaz.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = vaz.h2, metric = "H2")
vaz.diversity.zscore <- zscore_metric(obsval = net.metrics.diversity,
                                      nullval = vaz.diversity, metric = "Shannon diversity")
vaz.links.zscore <- zscore_metric(obsval = net.metrics.links,
                                  nullval = vaz.links, metric = "links per species")
vaz.density.zscore <- zscore_metric(obsval = net.metrics.density,
                                    nullval = vaz.density, metric = "linkage density")

# creating dataframes - r2dtable
r2d.connect <- as.data.frame(do.call('rbind', r2d.connect.zscore) )
r2d.connect <- r2d.connect %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.connect = connectance)

r2d.weightconnect <- as.data.frame(do.call('rbind', r2d.weightconnect.zscore) )
r2d.weightconnect <- r2d.weightconnect %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.weightconnect = `weighted connectance`)

r2d.nestedness <- as.data.frame(do.call('rbind', r2d.nest.zscore) )
r2d.nestedness <- r2d.nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.NODF = NODF)

r2d.weightnestedness <- as.data.frame(do.call('rbind', r2d.weightnest.zscore) )
r2d.weightnestedness <- r2d.weightnestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.weightedNODF = `weighted NODF`)

r2d.h2 <- as.data.frame(do.call('rbind', r2d.h2.zscore) )
r2d.h2 <- r2d.h2 %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.h2 = H2)

r2d.shannon <- as.data.frame(do.call('rbind', r2d.diversity.zscore) )
r2d.shannon <- r2d.shannon %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.shannon = `Shannon diversity`)

r2d.links <- as.data.frame(do.call('rbind', r2d.links.zscore)) 
r2d.links <- r2d.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.links = `links per species`)

r2d.density <- as.data.frame(do.call('rbind', r2d.density.zscore))
r2d.density <- r2d.density %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.density = `linkage density`)


# creating dataframes - vaznull
vaz.nestedness <- as.data.frame(do.call('rbind', vaz.nest.zscore) )
vaz.nestedness <- vaz.nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF = NODF)

vaz.weightnestedness <- as.data.frame(do.call('rbind', vaz.weightnest.zscore) )
vaz.weightnestedness <- vaz.weightnestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.weightedNODF = `weighted NODF`)

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

vaz.links <- as.data.frame(do.call('rbind', vaz.links.zscore)) 
vaz.links <- vaz.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.links = `links per species`)

vaz.density <- as.data.frame(do.call('rbind', vaz.density.zscore))
vaz.density <- vaz.density %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.density = `linkage density`)









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
# Calculate network nestedness - weighted
net.metrics.weightnest_noApis <- lapply(webs_noApis, networklevel, index = 'weighted NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noApis <- lapply(webs_noApis, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity_noApis <- lapply(webs_noApis, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density_noApis <- lapply(webs_noApis, networklevel, index = "linkage density")


# Make null models for all sites 
net.nulls.vaz_noApis <- lapply(webs_noApis, nullmodel, method = "vaznull", N = 500) # using the vaznull null - maintains connectance 
net.nulls.r2d_noApis <- lapply(webs_noApis, nullmodel, method = "r2dtable", N = 500) # using the r2dtable null - maintains row and column sums


# getting r2dtable null for each metric
r2d.connect_noApis <- net.null.connect(net.nulls.r2d_noApis) # can only use r2dtable null model for connectance -- vaznull maintains connectance
r2d.weightconnect_noApis <- net.null.weightconnect(net.nulls.r2d_noApis) # can only use r2dtable null model for connectance -- vaznull maintains connectance
r2d.nest_noApis <- net.null.nest(net.nulls.r2d_noApis)
r2d.weightnest_noApis <- net.null.weightnest(net.nulls.r2d_noApis)
r2d.h2_noApis <- net.null.h2(net.nulls.r2d_noApis)
r2d.diversity_noApis <- net.null.diversity(net.nulls.r2d_noApis)
r2d.links_noApis <- net.null.links(net.nulls.r2d_noApis)
r2d.density_noApis <- net.null.density(net.nulls.r2d_noApis)

# getting vaznull null for each metric
vaz.nest_noApis <- net.null.nest(net.nulls.vaz_noApis)
vaz.weightnest_noApis <- net.null.weightnest(net.nulls.vaz_noApis)
vaz.h2_noApis <- net.null.h2(net.nulls.vaz_noApis)
vaz.diversity_noApis <- net.null.diversity(net.nulls.vaz_noApis)
vaz.links_noApis <- net.null.links(net.nulls.vaz_noApis)
vaz.density_noApis <- net.null.density(net.nulls.vaz_noApis)


# getting z score - r2dtable
r2d.connect.zscore_noApis <- zscore_metric(obsval = net.metrics.connect_noApis, 
                                            nullval = r2d.connect_noApis, metric = "connectance")
r2d.weightconnect.zscore_noApis <- zscore_metric(obsval = net.metrics.weightconnect_noApis,
                                                 nullval = r2d.weightconnect_noApis, metric = "weighted connectance")
r2d.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis, 
                                        nullval = r2d.nest_noApis, metric = "NODF")
r2d.weightnest.zscore_noApis <- zscore_metric(obsval = net.metrics.weightnest_noApis, 
                                              nullval = r2d.weightnest_noApis, metric = "weighted NODF")
r2d.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis, 
                                      nullval = r2d.h2_noApis, metric = "H2")
r2d.diversity.zscore_noApis <- zscore_metric(obsval = net.metrics.diversity_noApis, 
                                             nullval = r2d.diversity_noApis, metric = "Shannon diversity")
r2d.links.zscore_noApis <- zscore_metric(obsval = net.metrics.links_noApis, 
                                         nullval = r2d.links_noApis, metric = "links per species") 
r2d.density.zscore_noApis <- zscore_metric(obsval = net.metrics.density_noApis, 
                                           nullval = r2d.density_noApis, metric = "linkage density") 

# getting z score - vaznull
vaz.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis, 
                                      nullval = vaz.nest_noApis, metric = "NODF")
vaz.weightnest.zscore_noApis <- zscore_metric(obsval = net.metrics.weightnest_noApis, 
                                                  nullval = vaz.weightnest_noApis, metric = "weighted NODF")
vaz.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis, 
                                  nullval = vaz.h2_noApis, metric = "H2")
vaz.diversity.zscore_noApis <- zscore_metric(obsval = net.metrics.diversity_noApis, 
                                                nullval = vaz.diversity_noApis, metric = "Shannon diversity")
vaz.links.zscore_noApis <- zscore_metric(obsval = net.metrics.links_noApis, 
                                        nullval = vaz.links_noApis, metric = "links per species") # CAN'T USE - NO SD FOR NULL MODEL, USE REAL VALUES IN ANALYSIS
vaz.density.zscore_noApis <- zscore_metric(obsval = net.metrics.density_noApis, 
                                            nullval = vaz.density_noApis, metric = "linkage density") 

# creating dataframes - r2dtable
r2d.connect_noApis <- as.data.frame(do.call('rbind', r2d.connect.zscore_noApis) )
r2d.connect_noApis <- r2d.connect_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.connect_noApis = connectance)

r2d.weightconnect_noApis <- as.data.frame(do.call('rbind', r2d.weightconnect.zscore_noApis) )
r2d.weightconnect_noApis <- r2d.weightconnect_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.weightconnect_noApis = `weighted connectance`)

r2d.nestedness_noApis <- as.data.frame(do.call('rbind', r2d.nest.zscore_noApis) )
r2d.nestedness_noApis <- r2d.nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.NODF_noApis = NODF)

r2d.weightnestedness_noApis <- as.data.frame(do.call('rbind', r2d.weightnest.zscore_noApis) )
r2d.weightnestedness_noApis <- r2d.weightnestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.weightedNODF_noApis = `weighted NODF`)

r2d.h2_noApis <- as.data.frame(do.call('rbind', r2d.h2.zscore_noApis) )
r2d.h2_noApis <- r2d.h2_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.h2_noApis = H2)

r2d.shannon_noApis <- as.data.frame(do.call('rbind', r2d.diversity.zscore_noApis) )
r2d.shannon_noApis <- r2d.shannon_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.shannon_noApis = `Shannon diversity`)

r2d.links_noApis <- as.data.frame(do.call('rbind', r2d.links.zscore_noApis)) 
r2d.links_noApis <- r2d.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.links_noApis = `links per species`)

r2d.density_noApis <- as.data.frame(do.call('rbind', r2d.density.zscore_noApis))
r2d.density_noApis <- r2d.density_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(r2d.density_noApis = `linkage density`)


# creating dataframes - vaznull
vaz.nestedness_noApis <- as.data.frame(do.call('rbind', vaz.nest.zscore_noApis) )
vaz.nestedness_noApis <- vaz.nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF_noApis = NODF)

vaz.weightnestedness_noApis <- as.data.frame(do.call('rbind', vaz.weightnest.zscore_noApis) )
vaz.weightnestedness_noApis <- vaz.weightnestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.weightedNODF_noApis = `weighted NODF`)

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

vaz.links_noApis <- as.data.frame(do.call('rbind', vaz.links.zscore_noApis)) 
vaz.links_noApis <- vaz.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.links_noApis = `links per species`)

vaz.density_noApis <- as.data.frame(do.call('rbind', vaz.density.zscore_noApis))
vaz.density_noApis <- vaz.density_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.density_noApis = `linkage density`)




#### exporting csv ####
# adding all together
network_metrics <- r2d.connect %>%
  left_join(r2d.weightconnect, by = c("block", "patch")) %>%
  left_join(r2d.nestedness, by = c("block", "patch")) %>%
  left_join(r2d.weightnestedness, by = c("block", "patch")) %>%
  left_join(r2d.h2, by = c("block", "patch")) %>%
  left_join(r2d.shannon, by = c("block", "patch")) %>%
  left_join(r2d.links, by = c("block", "patch")) %>%
  left_join(r2d.density, by = c("block", "patch")) %>%
  left_join(vaz.nestedness, by = c("block", "patch")) %>%
  left_join(vaz.weightnestedness, by = c("block", "patch")) %>%
  left_join(vaz.h2, by = c("block", "patch")) %>%
  left_join(vaz.shannon, by = c("block", "patch")) %>%
  left_join(vaz.links, by = c("block", "patch")) %>%
  left_join(vaz.density, by = c("block", "patch")) %>%
  left_join(r2d.connect_noApis, by = c("block", "patch")) %>%
  left_join(r2d.weightconnect_noApis, by = c("block", "patch")) %>%
  left_join(r2d.nestedness_noApis, by = c("block", "patch")) %>%
  left_join(r2d.weightnestedness_noApis, by = c("block", "patch")) %>%
  left_join(r2d.h2_noApis, by = c("block", "patch")) %>%
  left_join(r2d.shannon_noApis, by = c("block", "patch")) %>%
  left_join(r2d.links_noApis, by = c("block", "patch")) %>%
  left_join(r2d.density_noApis, by = c("block", "patch")) %>%
  left_join(vaz.nestedness_noApis, by = c("block", "patch")) %>%
  left_join(vaz.weightnestedness_noApis, by = c("block", "patch")) %>%
  left_join(vaz.h2_noApis, by = c("block", "patch")) %>%
  left_join(vaz.shannon_noApis, by = c("block", "patch")) %>%
  left_join(vaz.links_noApis, by = c("block", "patch")) %>%
  left_join(vaz.density_noApis, by = c("block", "patch")) 

# network_metrics <- nestedness %>%
#   left_join(h2, by = c("block", "patch")) %>%
#   left_join(shannon, by = c("block", "patch")) %>%
#   left_join(links, by = c("block", "patch")) %>%
#   left_join(density, by = c("block", "patch")) %>%
#   left_join(nestedness_noApis, by = c("block", "patch")) %>%
#   left_join(h2_noApis, by = c("block", "patch")) %>%
#   left_join(shannon_noApis, by = c("block", "patch")) %>%
#   left_join(links_noApis, by = c("block", "patch")) %>%
#   left_join(density_noApis, by = c("block", "patch")) #%>%
  # left_join(abundance, by = c("block", "patch")) %>%
  # left_join(abundance_noApis, by = c("block", "patch"))
# network_metrics$floral_diversity <- floral_diversity
# network_metrics$floral_div_noApis <- floral_div_noApis
# network_metrics$pollinator_diversity <- pollinator_diversity
# network_metrics$pollinator_div_noApis <- pollinator_div_noApis
# network_metrics$pollinator.rich.rare <- sRare_pollinator
# network_metrics$pollinator.rich.rare.noApis <- sRare_pollinator_noApis
# network_metrics$fl.rich <- fl.rich
# network_metrics$fl.rich.rare <- sRare
# network_metrics$fl.rich.rare.noApis <- sRare_noApis

write.csv(network_metrics, file = file.path("data", "network_metrics.csv"))










#### interaction rewiring and turnover ####
network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(webs[[1]]), "10.W" = as.matrix(webs[[2]]),
                                                         "52.B" = as.matrix(webs[[3]]), "52.W" = as.matrix(webs[[4]]),
                                                         "53N.B" = as.matrix(webs[[5]]), "53N.W" = as.matrix(webs[[6]]),
                                                         "53S.B" = as.matrix(webs[[7]]), "53S.W" = as.matrix(webs[[8]]),
                                                         "54S.B" = as.matrix(webs[[9]]), "54S.W" = as.matrix(webs[[10]]),
                                                         "57.B" = as.matrix(webs[[11]]), "57.W" = as.matrix(webs[[12]]),
                                                         "8.B" = as.matrix(webs[[13]]), "8.W" = as.matrix(webs[[14]]))), 
                                         index = "sorensen", partitioning="commondenom")
network_dissimilarity <- network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  #filter(block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
  #                        "54S.B-54S.W", "57.B-57.W", "8.B-8.W")) %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity")

network_dissimilarity$type <- factor(network_dissimilarity$type, levels=c("S", "WN", "ST", "OS"))


dissimilarity_plot <- network_dissimilarity %>%
  # filter(!type %in% c("S", "WN")) %>%
  filter(real_pairs == "real") %>%
  ggplot() +
  geom_boxplot(aes(type, dissimilarity, fill = type)) +
  #geom_jitter(aes(type, dissimilarity)) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  #scale_x_discrete(labels = c(expression(beta[WN]), expression(beta[ST]), expression(beta[OS]))) +
  # scale_x_discrete(labels = c(expression("Species turnover"), expression("Interaction rewiring"))) +
  ylim(c(0, 0.8)) +
  xlab("Dissimilarity component") +
  theme(legend.position = "none") +
  ylab(expression(paste("Dissimilarity between patch pairs"))) +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
dissimilarity_plot

# pdf(file = file.path("plots", "dissimilarity_plot.pdf"), width = 8, height = 9)
# dissimilarity_plot
# dev.off()

turnover.rewiring <- network_dissimilarity %>%
  filter(type %in% c("OS", "ST"))

m.dissimiliarty <- t.test(dissimilarity ~ type, data = turnover.rewiring, var.equal = TRUE)
m.dissimiliarty






