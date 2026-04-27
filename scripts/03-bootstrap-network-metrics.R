#remotes::install_github("valentinitnelav/bootstrapnet")
librarian::shelf(bootstrapnet, tidyverse, plyr, vegan, bipartite, data.table, glmmTMB, ggeffects, DHARMa)
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
# lst_h2 <- webs.matrix %>%
#   lapply(web_matrix_to_df) %>%
#   boot_networklevel(col_lower = "lower", # column name for plants
#                     col_higher = "higher", # column name for insects
#                     index = "H2",
#                     level = "both", 
#                     start = 90,
#                     step = 1,
#                     n_boot = 50,
#                     n_cpu = 4)
# 
# gg_networklevel(lst_h2)
# h2_boot <- lst_h2[["H2"]][["stats_df"]]
# h2_boot %>%
#   dplyr::count(spl_size) %>%
#   arrange(desc(n)) # 160 interactions is smallest web
# 
# h2_boot <- h2_boot %>%
#   filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
#   separate(web, into = c("block", "patch"))


#### bootstrapped NODF ####
lst_nodf <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "NODF",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_nodf, file = file.path("data", "L2_boot_metrics", "nodf.RData"))
#lst_nodf <- readRDS(file = file.path("data", "L2_boot_metrics", "nodf.RData"))

gg_networklevel(lst_nodf)
nodf_boot <- lst_nodf[["NODF"]][["stats_df"]]

nodf_boot <- nodf_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(nodf = mean) %>%
  dplyr::select(block, patch, spl_size, nodf)



#### bootstrapped links per species ####
lst_links <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "links per species",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_links, file = file.path("data", "L2_boot_metrics", "links.RData"))
#lst_links <- readRDS(file = file.path("data", "L2_boot_metrics", "links.RData"))


links_boot <- lst_links[["links per species"]][["stats_df"]]

links_boot <- links_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(links_per_sp = mean) %>%
  dplyr::select(block, patch, spl_size, links_per_sp)


#links_boot <- lst_links[["links per species"]][["lines_df"]]
# links_boot_mean <- links_boot_mean %>%
#   separate(web, into = c("block", "patch"), remove = F)
# links_boot <- links_boot %>%
#   separate(web, into = c("block", "patch"), remove = F)
# head(links_boot)


# m1 <- glmmTMB(value ~ patch + spl_size + (1|block) + (1|web), data = links_boot)
# summary(m1)
# plot(simulateResiduals(m1))
# 
# m1b_predict <- ggpredict(m1, terms = c("spl_size", "patch"), back_transform = T)
# 
# links_boot %>%
#   ggplot() +
#   geom_line(aes(spl_size, value, color = patch, group = interaction(web, simulation_id)), 
#             data = links_boot, linewidth = 0.1, alpha = 0.1) +
#   #geom_line(aes(spl_size, mean, color = patch, group = web), data = links_boot_mean) +
#   geom_line(aes(x, predicted, color = group), data = m1b_predict, linewidth = 2) +
#   geom_ribbon(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), data = m1b_predict, alpha = 0.2)


# # plotting
# m1b_predict_plot <- m1b_predict %>%
#   ggplot() +
#   #geom_jitter(aes(x = patch, y = value, color = patch), data = links_boot, size = 1, alpha = 0.2,
#   #            width = 0, height = 0) +
#   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
#                 data = m1b_predict, width = 0, linewidth = 2.5) +
#   geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
#   geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   scale_color_manual(values=c("#F5097C","#F7B3D4")) +
#   scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
#   xlab("Patch type") +
#   ylab(expression(paste("Links per species"))) +
#   theme_classic(base_size = 20) +
#   theme(legend.position = "none") 
# m1b_predict_plot




#### bootstrapped niche overlap ####
lst_niche <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "niche overlap",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_niche, file = file.path("data", "L2_boot_metrics", "niche.RData"))
#lst_niche <- readRDS(file = file.path("data", "L2_boot_metrics", "niche.RData"))

gg_networklevel(lst_niche)
HL_niche_boot <- lst_niche[["niche.overlap.HL"]][["stats_df"]]
LL_niche_boot <- lst_niche[["niche.overlap.LL"]][["stats_df"]]

HL_niche_boot <- HL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(HL_niche = mean) %>%
  dplyr::select(block, patch, spl_size, HL_niche)

LL_niche_boot <- LL_niche_boot %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch")) %>%
  dplyr::rename(LL_niche = mean) %>%
  dplyr::select(block, patch, spl_size, LL_niche)




#### bootstrapped robustness to extinction ####
# dividing into 7 lists because of computing power issues
webs.matrix1 <- webs.matrix[c(1, 2)]
webs.matrix2 <- webs.matrix[c(3, 4)]
webs.matrix3 <- webs.matrix[c(5, 6)]
webs.matrix4 <- webs.matrix[c(7, 8)]
webs.matrix5 <- webs.matrix[c(9, 10)]
webs.matrix6 <- webs.matrix[c(11, 12)]
webs.matrix7 <- webs.matrix[c(13, 14)]


lst_robustness1 <- webs.matrix1 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness1, file = file.path("data", "L2_boot_metrics", "robustness1.RData"))
#lst_robustness1 <- readRDS(file = file.path("data", "L2_boot_metrics", "robustness1.RData"))


lst_robustness2 <- webs.matrix2 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness2, file = file.path("data", "L2_boot_metrics", "robustness2.RData"))

lst_robustness3 <- webs.matrix3 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness3, file = file.path("data", "L2_boot_metrics", "robustness3.RData"))

lst_robustness4 <- webs.matrix4 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness4, file = file.path("data", "L2_boot_metrics", "robustness4.RData"))

lst_robustness5 <- webs.matrix5 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness5, file = file.path("data", "L2_boot_metrics", "robustness5.RData"))

lst_robustness6 <- webs.matrix6 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness6, file = file.path("data", "L2_boot_metrics", "robustness6.RData"))

lst_robustness7 <- webs.matrix7 %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 100,
                    n_cpu = 4,
                    weighted = F)
#saveRDS(lst_robustness7, file = file.path("data", "L2_boot_metrics", "robustness7.RData"))




gg_networklevel(lst_robustness1)
HL_robustness_boot1 <- lst_robustness1[["robustness.HL"]][["stats_df"]]
LL_robustness_boot1 <- lst_robustness1[["robustness.LL"]][["stats_df"]]

HL_robustness_boot1 <- HL_robustness_boot1 %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))

LL_robustness_boot1 <- LL_robustness_boot1 %>%
  filter(spl_size == 158) %>% # subsetting at 158 interactions to compare networks
  separate(web, into = c("block", "patch"))





#### All together ####
network_metrics_boot <- nodf_boot %>%
  left_join()



