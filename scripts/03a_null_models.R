# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, parallel, data.table)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 



#### ALL SPECIES ####
# splitting into webs
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
# Make null models for all sites 
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) # using the vaznull null - maintains connectance


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

# Make null models for all sites 
net.nulls.vaz_noApis <- lapply(webs_noApis, nullmodel, method = "vaznull", N = 500) # using the vaznull null - maintains connectance 



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


# Make null models for all sites 
net.nulls.vaz_noPoe <- lapply(webs_noPoe, nullmodel, method = "vaznull", N = 500) # using the vaznull null - maintains connectance 



# exporting as R data objects
save(net.nulls.vaz, net.nulls.vaz_noApis, net.nulls.vaz_noPoe, file = file.path("data", "nulls.RData"))
