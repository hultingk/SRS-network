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
# # Calculate network nestedness - weighted
# net.metrics.weightnest <- lapply(webs, networklevel, index = 'weighted NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity <- lapply(webs, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density <- lapply(webs, networklevel, index = "linkage density")
# modularity 
net.metrics.modularity <- lapply(webs, computeModules, method = "Beckett")

# Make null models for all sites 
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 1000) # using the vaznull null - maintains connectance 
net.nulls.swap <- lapply(webs, nullmodel, method = "swap.web", N = 1000) # using the swap.web null - maintains row and column sums


# getting swap.web null for each metric
swap.nest <- net.null.networklevel(nulls = net.nulls.swap, metric = "NODF")
# swap.weightnest <- net.null.weightnest(net.nulls.swap)
swap.h2 <- net.null.networklevel(nulls = net.nulls.swap, metric = "H2")
swap.diversity <- net.null.networklevel(nulls = net.nulls.swap, metric = "Shannon diversity")
swap.links <- net.null.networklevel(nulls = net.nulls.swap, metric = "links per species")
swap.density <- net.null.networklevel(nulls = net.nulls.swap, metric = "linkage density")

# getting vaznull null for each metric
vaz.nest <- net.null.networklevel(nulls = net.nulls.vaz, metric = "NODF")
#vaz.weightnest <- net.null.networklevel(nulls = net.nulls.vaz, metric = "weighted NODF")
vaz.h2 <- net.null.networklevel(nulls = net.nulls.vaz, metric = "H2")
vaz.diversity <- net.null.networklevel(nulls = net.nulls.vaz, metric = "Shannon diversity")
vaz.links <- net.null.networklevel(nulls = net.nulls.vaz, metric = "links per species") # NEED TO ADD ARGUEMENT FOR UPPER
vaz.density <- net.null.networklevel(nulls = net.nulls.vaz, metric = "linkage density")
vaz.module <- net.null.computeModule(nulls = net.nulls.vaz)


# getting z score - swap.web
swap.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = swap.nest, metric = "NODF")
# swap.weightnest.zscore <- zscore_metric(obsval = net.metrics.weightnest,
#                                        nullval = swap.weightnest, metric = "weighted NODF")
swap.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = swap.h2, metric = "H2")
swap.diversity.zscore <- zscore_metric(obsval = net.metrics.diversity,
                                      nullval = swap.diversity, metric = "Shannon diversity")
swap.links.zscore <- zscore_metric(obsval = net.metrics.links,
                                  nullval = swap.links, metric = "links per species")
swap.density.zscore <- zscore_metric(obsval = net.metrics.density,
                                    nullval = swap.density, metric = "linkage density")

# getting z score - vaznull
vaz.nest.zscore <- zscore_metric(obsval = net.metrics.nest,
                                 nullval = vaz.nest, metric = "NODF")
# vaz.weightnest.zscore <- zscore_metric(obsval = net.metrics.weightnest,
#                                        nullval = vaz.weightnest, metric = "weighted NODF")
vaz.h2.zscore <- zscore_metric(obsval = net.metrics.h2,
                               nullval = vaz.h2, metric = "H2")
vaz.diversity.zscore <- zscore_metric(obsval = net.metrics.diversity,
                                      nullval = vaz.diversity, metric = "Shannon diversity")
vaz.links.zscore <- zscore_metric(obsval = net.metrics.links,
                                  nullval = vaz.links, metric = "links per species")
vaz.density.zscore <- zscore_metric(obsval = net.metrics.density,
                                    nullval = vaz.density, metric = "linkage density")

# # creating dataframes - swap.web


swap.nestedness <- as.data.frame(do.call('rbind', swap.nest.zscore) )
swap.nestedness <- swap.nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.NODF = NODF)

# swap.weightnestedness <- as.data.frame(do.call('rbind', swap.weightnest.zscore) )
# swap.weightnestedness <- swap.weightnestedness %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(swap.weightedNODF = `weighted NODF`)

swap.h2 <- as.data.frame(do.call('rbind', swap.h2.zscore) )
swap.h2 <- swap.h2 %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.h2 = H2)

swap.shannon <- as.data.frame(do.call('rbind', swap.diversity.zscore) )
swap.shannon <- swap.shannon %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.shannon = `Shannon diversity`)

swap.links <- as.data.frame(do.call('rbind', swap.links.zscore))
swap.links <- swap.links %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.links = `links per species`)

swap.density <- as.data.frame(do.call('rbind', swap.density.zscore))
swap.density <- swap.density %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.density = `linkage density`)


# creating dataframes - vaznull
vaz.nestedness <- as.data.frame(do.call('rbind', vaz.nest.zscore) )
vaz.nestedness <- vaz.nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF = NODF)

# vaz.weightnestedness <- as.data.frame(do.call('rbind', vaz.weightnest.zscore) )
# vaz.weightnestedness <- vaz.weightnestedness %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(vaz.weightedNODF = `weighted NODF`)

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


# creating dataframes for connectance
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
# # Calculate network nestedness - weighted
# net.metrics.weightnest_noApis <- lapply(webs_noApis, networklevel, index = 'weighted NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noApis <- lapply(webs_noApis, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity_noApis <- lapply(webs_noApis, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density_noApis <- lapply(webs_noApis, networklevel, index = "linkage density")


# Make null models for all sites 
net.nulls.vaz_noApis <- lapply(webs_noApis, nullmodel, method = "vaznull", N = 1000) # using the vaznull null - maintains connectance 
#net.nulls.swap_noApis <- lapply(webs_noApis, nullmodel, method = "swap.web", N = 1000) # using the swap.web null - maintains row and column sums


# getting swap.web null for each metric
swap.nest_noApis <- net.null.nest(net.nulls.swap_noApis)
# swap.weightnest_noApis <- net.null.weightnest(net.nulls.swap_noApis)
swap.h2_noApis <- net.null.h2(net.nulls.swap_noApis)
swap.diversity_noApis <- net.null.diversity(net.nulls.swap_noApis)
swap.links_noApis <- net.null.links(net.nulls.swap_noApis)
swap.density_noApis <- net.null.density(net.nulls.swap_noApis)

# getting vaznull null for each metric
vaz.nest_noApis <- net.null.nest(net.nulls.vaz_noApis)
#vaz.weightnest_noApis <- net.null.weightnest(net.nulls.vaz_noApis)
vaz.h2_noApis <- net.null.h2(net.nulls.vaz_noApis)
vaz.diversity_noApis <- net.null.diversity(net.nulls.vaz_noApis)
vaz.links_noApis <- net.null.links(net.nulls.vaz_noApis)
vaz.density_noApis <- net.null.density(net.nulls.vaz_noApis)


# getting z score - swap.web

swap.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis,
                                        nullval = swap.nest_noApis, metric = "NODF")
# swap.weightnest.zscore_noApis <- zscore_metric(obsval = net.metrics.weightnest_noApis,
#                                               nullval = swap.weightnest_noApis, metric = "weighted NODF")
swap.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis,
                                      nullval = swap.h2_noApis, metric = "H2")
swap.diversity.zscore_noApis <- zscore_metric(obsval = net.metrics.diversity_noApis,
                                             nullval = swap.diversity_noApis, metric = "Shannon diversity")
swap.links.zscore_noApis <- zscore_metric(obsval = net.metrics.links_noApis,
                                         nullval = swap.links_noApis, metric = "links per species")
swap.density.zscore_noApis <- zscore_metric(obsval = net.metrics.density_noApis,
                                           nullval = swap.density_noApis, metric = "linkage density")

# getting z score - vaznull
vaz.nest.zscore_noApis <- zscore_metric(obsval = net.metrics.nest_noApis, 
                                      nullval = vaz.nest_noApis, metric = "NODF")
# vaz.weightnest.zscore_noApis <- zscore_metric(obsval = net.metrics.weightnest_noApis, 
#                                                   nullval = vaz.weightnest_noApis, metric = "weighted NODF")
vaz.h2.zscore_noApis <- zscore_metric(obsval = net.metrics.h2_noApis, 
                                  nullval = vaz.h2_noApis, metric = "H2")
vaz.diversity.zscore_noApis <- zscore_metric(obsval = net.metrics.diversity_noApis, 
                                                nullval = vaz.diversity_noApis, metric = "Shannon diversity")
vaz.links.zscore_noApis <- zscore_metric(obsval = net.metrics.links_noApis, 
                                        nullval = vaz.links_noApis, metric = "links per species") # CAN'T USE - NO SD FOR NULL MODEL, USE REAL VALUES IN ANALYSIS
vaz.density.zscore_noApis <- zscore_metric(obsval = net.metrics.density_noApis, 
                                            nullval = vaz.density_noApis, metric = "linkage density") 

# creating dataframes - swap.web

swap.nestedness_noApis <- as.data.frame(do.call('rbind', swap.nest.zscore_noApis) )
swap.nestedness_noApis <- swap.nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.NODF_noApis = NODF)

# swap.weightnestedness_noApis <- as.data.frame(do.call('rbind', swap.weightnest.zscore_noApis) )
# swap.weightnestedness_noApis <- swap.weightnestedness_noApis %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(swap.weightedNODF_noApis = `weighted NODF`)

swap.h2_noApis <- as.data.frame(do.call('rbind', swap.h2.zscore_noApis) )
swap.h2_noApis <- swap.h2_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.h2_noApis = H2)

swap.shannon_noApis <- as.data.frame(do.call('rbind', swap.diversity.zscore_noApis) )
swap.shannon_noApis <- swap.shannon_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.shannon_noApis = `Shannon diversity`)

swap.links_noApis <- as.data.frame(do.call('rbind', swap.links.zscore_noApis))
swap.links_noApis <- swap.links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.links_noApis = `links per species`)

swap.density_noApis <- as.data.frame(do.call('rbind', swap.density.zscore_noApis))
swap.density_noApis <- swap.density_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(swap.density_noApis = `linkage density`)


# creating dataframes - vaznull
vaz.nestedness_noApis <- as.data.frame(do.call('rbind', vaz.nest.zscore_noApis) )
vaz.nestedness_noApis <- vaz.nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.NODF_noApis = NODF)

# vaz.weightnestedness_noApis <- as.data.frame(do.call('rbind', vaz.weightnest.zscore_noApis) )
# vaz.weightnestedness_noApis <- vaz.weightnestedness_noApis %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(vaz.weightedNODF_noApis = `weighted NODF`)

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

# creating connectance dataframes
# net.connect_noApis <- as.data.frame(do.call('rbind', net.metrics.connect_noApis) )
# net.connect_noApis <- net.connect_noApis %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(connectance_noApis = connectance)
# 
# net.weightconnect_noApis <- as.data.frame(do.call('rbind', net.metrics.weightconnect_noApis) )
# net.weightconnect_noApis <- net.weightconnect_noApis %>%
#   rownames_to_column(var = "unique_ID") %>%
#   separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
#   dplyr::rename(weightconnectance_noApis = `weighted connectance`)


#### exporting csv ####
# adding all together
network_metrics <- net.connect %>%
  left_join(net.weightconnect, by = c("block", "patch")) %>%
  left_join(swap.nestedness, by = c("block", "patch")) %>%
  # left_join(swap.weightnestedness, by = c("block", "patch")) %>%
  left_join(swap.h2, by = c("block", "patch")) %>%
  left_join(swap.shannon, by = c("block", "patch")) %>%
  left_join(swap.links, by = c("block", "patch")) %>%
  left_join(swap.density, by = c("block", "patch")) %>%
  left_join(vaz.nestedness, by = c("block", "patch")) %>%
  # left_join(vaz.weightnestedness, by = c("block", "patch")) %>%
  left_join(vaz.h2, by = c("block", "patch")) %>%
  left_join(vaz.shannon, by = c("block", "patch")) %>%
  left_join(vaz.links, by = c("block", "patch")) %>%
  left_join(vaz.density, by = c("block", "patch"))# %>%
  # left_join(net.connect_noApis, by = c("block", "patch")) %>%
  # left_join(net.weightconnect_noApis, by = c("block", "patch")) %>%
  # left_join(swap.nestedness_noApis, by = c("block", "patch")) %>%
  # left_join(swap.weightnestedness_noApis, by = c("block", "patch")) %>%
 #  left_join(swap.h2_noApis, by = c("block", "patch")) %>%
 #  left_join(swap.shannon_noApis, by = c("block", "patch")) %>%
 #  left_join(swap.links_noApis, by = c("block", "patch")) %>%
 #  left_join(swap.density_noApis, by = c("block", "patch")) %>%
 #  left_join(vaz.nestedness_noApis, by = c("block", "patch")) %>%
 # # left_join(vaz.weightnestedness_noApis, by = c("block", "patch")) %>%
 #  left_join(vaz.h2_noApis, by = c("block", "patch")) %>%
 #  left_join(vaz.shannon_noApis, by = c("block", "patch")) %>%
 #  left_join(vaz.links_noApis, by = c("block", "patch")) %>%
 #  left_join(vaz.density_noApis, by = c("block", "patch")) 


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






