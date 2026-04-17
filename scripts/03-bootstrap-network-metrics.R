#remotes::install_github("valentinitnelav/bootstrapnet")
librarian::shelf(bootstrapnet, tidyverse, plyr, vegan, bipartite, data.table, glmmTMB, ggeffects)
set.seed(100)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))

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

webs.matrix <- lapply(webs, as.matrix)

#### bootstrapped H2 ####
lst_h2 <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "H2",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 50,
                    n_cpu = 4)

gg_networklevel(lst_h2)
h2_boot <- lst_h2[["H2"]][["stats_df"]]
h2_boot %>%
  dplyr::count(spl_size) %>%
  arrange(desc(n)) # 160 interactions is smallest web

h2_boot <- h2_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))


#### bootstrapped NODF ####
lst_nodf <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "NODF",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 50,
                    n_cpu = 4)

gg_networklevel(lst_nodf)
nodf_boot <- lst_nodf[["NODF"]][["stats_df"]]

nodf_boot <- nodf_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))


#### bootstrapped links per species ####
lst_links <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "links per species",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 50,
                    n_cpu = 4)

gg_networklevel(lst_links)
links_boot <- lst_links[["links per species"]][["stats_df"]]

links_boot <- links_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))



#### bootstrapped niche overlap ####
lst_niche <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "niche overlap",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 50,
                    n_cpu = 4)

gg_networklevel(lst_niche)
HL_niche_boot <- lst_niche[["niche.overlap.HL"]][["stats_df"]]
LL_niche_boot <- lst_niche[["niche.overlap.LL"]][["stats_df"]]

HL_niche_boot <- HL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))

LL_niche_boot <- LL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))













