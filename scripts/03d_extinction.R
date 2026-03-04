# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

set.seed(100)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))
load(file = file.path("data", "L2_nulls", "nulls.RData"))


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 


#### extinction: WITH ALL SPECIES ####
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


# all species - lower level
net.extinct.LL <- lapply(webs, second.extinct, participant = 'lower', method = "abun") 
net.robustness.LL <- lapply(net.extinct.LL, robustness)

# apply to null networks
vaz.extinct.LL <- net.null.extinct(net.nulls.vaz, participant = "lower", method = "abun")

# zscore 
vaz.extinct.zscore.LL <- zscore_metric(obsval = net.robustness.LL,
                                 nullval = vaz.extinct.LL)
# getting into dataframe format
net.extinct.LL <- as.data.frame(do.call('rbind', vaz.extinct.zscore.LL) )
net.extinct.LL <- net.extinct.LL %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(extinct.robustness.LL = V1)


# all species - higher level
net.extinct.HL <- lapply(webs, second.extinct, participant = 'higher', method = "abun") 
net.robustness.HL <- lapply(net.extinct.HL, robustness)

# apply to null networks
vaz.extinct.HL <- net.null.extinct(net.nulls.vaz, participant = "higher", method = "abun")

# zscore 
vaz.extinct.zscore.HL <- zscore_metric(obsval = net.robustness.HL,
                                       nullval = vaz.extinct.HL)
# getting into dataframe format
net.extinct.HL <- as.data.frame(do.call('rbind', vaz.extinct.zscore.HL) )
net.extinct.HL <- net.extinct.HL %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(extinct.robustness.HL = V1)





#### extinction: Excluding Apis mellifera and Poecilognathus sulphureus ####
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


# excluding apis and poe - lower level
net.extinct.LL_noPoe_Apis <- lapply(webs_noPoe_Apis, second.extinct, participant = 'lower', method = "abun") 
net.robustness.LL_noPoe_Apis <- lapply(net.extinct.LL_noPoe_Apis, robustness)

# apply to null networks
vaz.extinct.LL_noPoe_Apis <- net.null.extinct(net.nulls.vaz_noPoe_Apis, participant = 'lower', method = "abun")
# zscore 
vaz.extinct.zscore.LL_noPoe_Apis <- zscore_metric(obsval = net.robustness.LL_noPoe_Apis,
                                    nullval = vaz.extinct.LL_noPoe_Apis)
# getting into dataframe format
net.extinct.LL_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.extinct.zscore.LL_noPoe_Apis) )
net.extinct.LL_noPoe_Apis <- net.extinct.LL_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(extinct.robustness.LL_noPoe_Apis = V1)



# excluding apis and poe - higher level
net.extinct.HL_noPoe_Apis <- lapply(webs_noPoe_Apis, second.extinct, participant = 'higher', method = "abun") 
net.robustness.HL_noPoe_Apis <- lapply(net.extinct.HL_noPoe_Apis, robustness)

# apply to null networks
vaz.extinct.HL_noPoe_Apis <- net.null.extinct(net.nulls.vaz_noPoe_Apis, participant = 'higher', method = "abun")
# zscore 
vaz.extinct.zscore.HL_noPoe_Apis <- zscore_metric(obsval = net.robustness.HL_noPoe_Apis,
                                                  nullval = vaz.extinct.HL_noPoe_Apis)
# getting into dataframe format
net.extinct.HL_noPoe_Apis <- as.data.frame(do.call('rbind', vaz.extinct.zscore.HL_noPoe_Apis) )
net.extinct.HL_noPoe_Apis <- net.extinct.HL_noPoe_Apis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(extinct.robustness.HL_noPoe_Apis = V1)


#### all together ####
extinction_metric <- net.extinct.LL %>%
  left_join(net.extinct.LL_noPoe_Apis, by = c("block", "patch")) %>%
  left_join(net.extinct.HL, by = c("block", "patch")) %>%
  left_join(net.extinct.HL_noPoe_Apis, by = c("block", "patch"))
  

write.csv(extinction_metric, file = file.path("data", "L4_metrics", "extinction_metric.csv"))
