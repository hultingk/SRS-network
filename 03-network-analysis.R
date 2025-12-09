# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite)


# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


# removing removing non-IDed species, creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 

pollinator <- pollinator %>%
  filter(pollinator_species != "Poecilognathus sulphureus")

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


# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity <- lapply(webs, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density <- lapply(webs, networklevel, index = "linkage density")


# Make null models for all sites using the vaznull null
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) 


# getting null for each metric
vaz.nest <- net.null.nest(net.nulls.vaz)
vaz.h2 <- net.null.h2(net.nulls.vaz)
vaz.diversity <- net.null.diversity(net.nulls.vaz)
vaz.links <- net.null.links(net.nulls.vaz)
vaz.density <- net.null.density(net.nulls.vaz)


# getting z score
vaz.nest.zscore <- nest.zscore(vaz.nest)
vaz.h2.zscore <- h2.zscore(vaz.h2)
vaz.diversity.zscore <- diversity.zscore(vaz.diversity)
vaz.links.zscore <- links.zscore(vaz.links) # CAN'T USE - NO SD FOR NULL MODEL, USE REAL VALUES IN ANALYSIS
vaz.density.zscore <- density.zscore(vaz.density) 

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

density <- as.data.frame(do.call('rbind', net.metrics.density))
density <- density %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column


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
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")
names(webs_noApis) <- webs.names


# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest_noApis <- lapply(webs_noApis, networklevel, index = 'NODF') 
# Calculate network specialization for all plant-pollinator sites
net.metrics.h2_noApis <- lapply(webs_noApis, networklevel, index = 'H2') 
# Calculate network interaction diversity for all plant-pollinator sites
net.metrics.diversity_noApis <- lapply(webs_noApis, networklevel, index = 'Shannon diversity') 
# Calculate links per pollinator species for all plant-pollinator sites
net.metrics.links_noApis <- lapply(webs_noApis, networklevel, index = 'links per species', level = "higher") 
# linkage density
net.metrics.density_noApis <- lapply(webs_noApis, networklevel, index = "linkage density")


# Make null models for all sites using the vaznull null
net.nulls.vaz_noApis <- lapply(webs_noApis, nullmodel, method = "vaznull", N = 500) 

# null metrics
vaz.nest_noApis <- net.null.nest(net.nulls.vaz_noApis)
vaz.h2_noApis <- net.null.h2(net.nulls.vaz_noApis)
vaz.diversity_noApis <- net.null.diversity(net.nulls.vaz_noApis)
vaz.links_noApis <- net.null.links(net.nulls.vaz_noApis)
vaz.density_noApis <- net.null.density(net.nulls.vaz_noApis)


# z score
vaz.nest.zscore_noApis <- nest.zscore(vaz.nest_noApis)
vaz.h2.zscore_noApis <- h2.zscore(vaz.h2_noApis)
vaz.diversity.zscore_noApis <- diversity.zscore(vaz.diversity_noApis)
vaz.links.zscore_noApis <- links.zscore(vaz.links_noApis) # CAN'T USE - NO SD FOR NULL MODEL, USE REAL VALUES IN ANALYSIS
vaz.density.zscore_noApis <- density.zscore(vaz.density_noApis) 

# creating dataframes
nestedness_noApis <- as.data.frame(do.call('rbind', vaz.nest.zscore_noApis) )
nestedness_noApis <- nestedness_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  dplyr::rename(NODF_noApis = NODF) %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into columns

h2_noApis <- as.data.frame(do.call('rbind', vaz.h2.zscore_noApis) )
h2_noApis <- h2_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  dplyr::rename(h2_noApis = H2) %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

shannon_noApis <- as.data.frame(do.call('rbind', vaz.diversity.zscore_noApis) )
shannon_noApis <- shannon_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  dplyr::rename(shannon_noApis = `Shannon diversity`) %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

links_noApis <- as.data.frame(do.call('rbind', net.metrics.links_noApis))
links_noApis <- links_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  dplyr::rename(links_noApis = `links per species`) %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

density_noApis <- as.data.frame(do.call('rbind', net.metrics.density_noApis))
density_noApis <- density_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  dplyr::rename(density_noApis = `linkage density`) %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column








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










#### Pollinator abundance ####
# total pollinator visition
abundance <- pollinator %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance = n) %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  dplyr::select(!c("n"))
# pollinator visitation excuding Apis mellifera
abundance_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID) %>%
  mutate(abundance_noApis = n) %>%
  separate(unique_ID, into = c("block", "patch")) %>%
  dplyr::select(!c("n"))


#### Flower and pollinator species diversity ####
## flower diversity
floral_wider <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_diversity <- diversity(floral_wider, "shannon") # calculating diversity of flowers that are interacted with

# floral diversity with no Apis
floral_wider_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
floral_div_noApis <- diversity(floral_wider_noApis, "shannon") # calculating diversity of flowers that are interacted with


# floral richness - rarefied
spAbund <- rowSums(floral_wider) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare <- rarefy(floral_wider, 160) # now use function rarefy

# floral richness - rarefied no apis
spAbund_noApis <- rowSums(floral_wider_noApis) # calculating minimum # of observation 
min(spAbund_noApis) # 7 is the fewest interactions observed per patch
sRare_noApis <- rarefy(floral_wider_noApis, 141) # now use function rarefy


# floral richness - not rarefied
floral_wider[floral_wider>0] <- 1 
fl.rich <- rowSums(floral_wider)



## pollinator diversity
pollinator_wider <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_diversity <- diversity(pollinator_wider, "shannon") 

# pollinator diversity no apis
pollinator_wider_noApis <- pollinator %>%
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID for diversity measure
pollinator_div_noApis <- diversity(pollinator_wider_noApis, "shannon") 



# rarified pollinator richness
spAbund <- rowSums(pollinator_wider) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare_pollinator <- rarefy(pollinator_wider, 160) # now use function rarefy

# rarified pollinator richness no apis
spAbund <- rowSums(pollinator_wider_noApis) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare_pollinator_noApis <- rarefy(pollinator_wider_noApis, 141) # now use function rarefy

#### exporting csv ####
# adding all together
network_metrics <- nestedness %>%
  left_join(h2, by = c("block", "patch")) %>%
  left_join(shannon, by = c("block", "patch")) %>%
  left_join(links, by = c("block", "patch")) %>%
  left_join(density, by = c("block", "patch")) %>%
  left_join(nestedness_noApis, by = c("block", "patch")) %>%
  left_join(h2_noApis, by = c("block", "patch")) %>%
  left_join(shannon_noApis, by = c("block", "patch")) %>%
  left_join(links_noApis, by = c("block", "patch")) %>%
  left_join(density_noApis, by = c("block", "patch")) %>%
  left_join(abundance, by = c("block", "patch")) %>%
  left_join(abundance_noApis, by = c("block", "patch"))
network_metrics$floral_diversity <- floral_diversity
network_metrics$floral_div_noApis <- floral_div_noApis
network_metrics$pollinator_diversity <- pollinator_diversity
network_metrics$pollinator_div_noApis <- pollinator_div_noApis
network_metrics$pollinator.rich.rare <- sRare_pollinator
network_metrics$pollinator.rich.rare.noApis <- sRare_pollinator_noApis
network_metrics$fl.rich <- fl.rich
network_metrics$fl.rich.rare <- sRare
network_metrics$fl.rich.rare.noApis <- sRare_noApis

write.csv(network_metrics, file = file.path("data", "network_metrics.csv"))










