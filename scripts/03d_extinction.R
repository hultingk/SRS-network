# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table)

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


net.extinct <- lapply(webs, second.extinct, participant = 'lower') 
net.robustness <- lapply(net.extinct, robustness)

# apply to null networks
vaz.extinct <- net.null.extinct(net.nulls.vaz)
# zscore 
vaz.extinct.zscore <- zscore_metric(obsval = net.robustness,
                                 nullval = vaz.extinct)

net.extinct <- as.data.frame(do.call('rbind', vaz.extinct.zscore) )
net.extinct <- net.extinct %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch")) %>% # seperating unique ID into columns
  dplyr::rename(extinct.robustness = V1)


