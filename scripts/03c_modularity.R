# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, parallel, data.table)

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

### All species ####
# modularity 
net.metrics.modularity <- lapply(webs, networklevel, index = 'modularity', level = "both")
# null model modularity  - takes a very long time
vaz.module <- net.null.networklevel(nulls = net.nulls.vaz, metric = "modularity", level = "both", web.names = web.names)
#save(vaz.module, file = file.path("data", "L3_modularity", "vaz.module.RData"))
# z score
vaz.module.zscore <- zscore_metric(obsval = net.metrics.modularity,
                                   nullval = vaz.module, metric = "modularity Q")
# make dataframe
vaz.module.df <- as.data.frame(do.call('rbind', vaz.module.zscore)) 
vaz.module.df <- vaz.module.df %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.module = `modularity Q`)



#### NO APIS ####
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

# modularity 
net.metrics.modularity_noApis <- lapply(webs_noApis, networklevel, index = 'modularity', level = "both") 
# null model modularity
vaz.module_noApis <- net.null.networklevel(nulls = net.nulls.vaz_noApis, metric = "modularity", level = "both")
#save(vaz.module_noApis, file = file.path("data", "L3_modularity","vaz.module_noApis.RData"))
# z score
vaz.module.zscore_noApis <- zscore_metric(obsval = net.metrics.modularity_noApis,
                                          nullval = vaz.module_noApis, metric = "modularity Q")
# creating dataframe
vaz.module_noApis <- as.data.frame(do.call('rbind', vaz.module.zscore_noApis)) 
vaz.module_noApis <- vaz.module_noApis %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.module_noApis = `modularity Q`)


#### NO Poecilognathus sulphureus ####
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

# modularity 
net.metrics.modularity_noPoe <- lapply(webs_noPoe, networklevel, index = 'modularity', level = "both") 
# null model modularity 
vaz.module_noPoe <- net.null.networklevel(nulls = net.nulls.vaz_noPoe, metric = "modularity", level = "both") 
#save(vaz.module_noPoe, file = file.path("data", "L3_modularity","vaz.module_noPoe.RData"))
# z score
vaz.module.zscore_noPoe <- zscore_metric(obsval = net.metrics.modularity_noPoe,
                                         nullval = vaz.module_noPoe, metric = "modularity Q")
# dataframe
vaz.module_noPoe <- as.data.frame(do.call('rbind', vaz.module.zscore_noPoe)) 
vaz.module_noPoe <- vaz.module_noPoe %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(vaz.module_noPoe = `modularity Q`)









