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

pollinator %>%
  dplyr::count(pollinator_species, flower_species)


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

# bootstrapped h2
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
# 160 interactions is smallest web
h2_boot <- lst_h2[["H2"]][["stats_df"]]
h2_boot %>%
  dplyr::count(spl_size) %>%
  arrange(desc(n))

h2_boot <- h2_boot %>%
  filter(spl_size == 158) %>%
  separate(web, into = c("block", "patch"))

m1 <- glmmTMB(mean ~ patch + (1|block), 
                            data = h2_boot)
summary(m1)
m1_predict <- ggpredict(m1, terms = c("patch"), back_transform = T)
# plotting
m1_predict_plot <- m1_predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = mean, color = patch), data = h2_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m1_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("H2"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m1_predict_plot



### bootstrapped taking dominant species out
pollinator_split_subset <- pollinator %>%
  dplyr::filter(!pollinator_species %in% c("Poecilognathus sulphureus", "Apis mellifera")) %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 

# getting each network into correct format
webs_subset <- pollinator_split_subset %>%
  lapply(prepare_matrix)

# adding patch names to each network
names(webs_subset) <- webs.names

webs.matrix_subset <- lapply(webs_subset, as.matrix)

lst_h2_subset <- webs.matrix_subset %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "H2",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 50,
                    n_cpu = 4)

gg_networklevel(lst_h2_subset)


h2_boot_subset <- lst_h2_subset[["H2"]][["stats_df"]]
h2_boot_subset %>%
  dplyr::count(spl_size) %>%
  arrange(desc(n))

h2_boot_subset <- h2_boot_subset %>%
  filter(spl_size == 103) %>%
  separate(web, into = c("block", "patch"))

m1b <- glmmTMB(mean ~ patch + (1|block), 
              data = h2_boot_subset)
summary(m1b)
m1b_predict <- ggpredict(m1b, terms = c("patch"), back_transform = T)
# plotting
m1b_predict_plot <- m1b_predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = mean, color = patch), data = h2_boot_subset, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m1b_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("H2"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m1b_predict_plot

