# -------------------------------------- #
#### Analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr)


# load data
network_metrics <- read.csv(file = file.path("data", "network_metrics.csv"))
diversity_metrics <- read.csv(file = file.path("data", "diversity_metrics.csv"))
network_vaznull <- read.csv(file = file.path("data", "network_vaznull.csv"))

#### connectance ####
# no vaznull metric for connectance
# unconnected patches have higher connectance, even excluding Apis and Poe
m.connect <- glmmTMB(connectance ~ patch + (1|block),
                         data = network_vaznull)
summary(m.connect)
m.connect_noApis <- glmmTMB(connectance_noApis ~ patch + (1|block),
                         data = network_vaznull)
summary(m.connect_noApis)
m.connect_noPoe <- glmmTMB(connectance_noPoe ~ patch + (1|block),
                            data = network_vaznull)
summary(m.connect_noPoe)

#### weighted connectance ####
# no vaznull metric for weighted connectance
# no difference in weighted connectance including or excluding api and poe
m.weightconnect <- glmmTMB(weightconnectance ~ patch + (1|block),
                     data = network_vaznull)
summary(m.weightconnect)
m.weightconnect_noApis <- glmmTMB(weightconnectance_noApis ~ patch + (1|block),
                            data = network_vaznull)
summary(m.weightconnect_noApis)
m.weightconnect_noPoe <- glmmTMB(weightconnectance_noPoe ~ patch + (1|block),
                           data = network_vaznull)
summary(m.weightconnect_noPoe)



#### nestedness ####
m.nestedness.vaz <- glmmTMB(vaz.NODF ~ patch + (1|block), # significantly lower in connected patches 
                            data = network_vaznull)
summary(m.nestedness.vaz)

m.nestedness.vaz_noApis <- glmmTMB(vaz.NODF_noApis ~ patch + (1|block), # no difference
                                   data = network_vaznull)
summary(m.nestedness.vaz_noApis)

m.nestedness.vaz_noPoe <- glmmTMB(vaz.NODF_noPoe ~ patch + (1|block), # no difference
                                   data = network_vaznull)
summary(m.nestedness.vaz_noPoe)

# # plotting NO APIS
# # model predictions for plotting
# m.nestedness.noApis.df <- network_metrics %>%
#   dplyr::select(c("block", "patch", "NODF_noApis"))
# m.nestedness.noApis.df$nestedness_pred <- predict(m.nestedness.noApis, re.form = NA)
# m.nestedness.noApis.df$patch <- factor(m.nestedness.noApis.df$patch, levels = c("B", "W"))
# # plotting
# nestedness.pred <- m.nestedness.noApis.df %>%
#   ggplot() +
#   geom_line(aes(x = patch, y = NODF_noApis, group = block), linewidth = 2.5, color = "grey25", alpha = 0.85) +
#   geom_line(aes(x = patch, y = nestedness_pred, group = block), linewidth = 5.5, linetype = 3) +
#   geom_point(aes(x = patch, y = NODF_noApis, color = patch), size = 11, alpha = 0.95) + 
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   scale_color_manual(values=c("#506D8F","#E2A03C")) +
#   xlab("Patch Type") +
#   ylab(expression(paste("NODF (z-score)"))) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size = 30)) + # axis tick mark size
#   theme(axis.title = element_text(size = 34)) #+ # axis label size
# nestedness.pred

# pdf(file = file.path("plots", "nestedness_noApis.pdf"), width = 8, height = 10)
# nestedness.pred
# dev.off()



#### H2' ####
m.h2.vaz <- glmmTMB(vaz.h2 ~ patch + (1|block), # significantly lower in connected patches 
                            data = network_vaznull)
summary(m.h2.vaz)

m.h2.vaz_noApis <- glmmTMB(vaz.h2_noApis ~ patch + (1|block), # no difference
                                   data = network_vaznull)
summary(m.h2.vaz_noApis)

m.h2.vaz_noPoe <- glmmTMB(vaz.h2_noPoe ~ patch + (1|block), # no difference
                                  data = network_vaznull)
summary(m.h2.vaz_noPoe)


# 
# (2.1497/5.4302)*100 # 39.58786 less nested
# # plotting
# # model predictions for plotting
# m.h2.noApis.df <- network_metrics %>%
#   dplyr::select(c("block", "patch", "h2_noApis"))
# m.h2.noApis.df$h2.pred <- predict(m.h2.noApis, re.form = NA)
# m.h2.noApis.df$patch <- factor(m.h2.noApis.df$patch, levels = c("B", "W"))
# # plotting
# h2.pred <- m.h2.noApis.df %>%
#   ggplot() +
#   geom_line(aes(x = patch, y = h2_noApis, group = block), linewidth = 2.5, color = "lightslategrey", alpha = 0.5) +
#   geom_line(aes(x = patch, y = h2.pred, group = block), linewidth = 5, linetype = 1) +
#   geom_point(aes(x = patch, y = h2_noApis, color = patch), size = 11, alpha = 0.8) + 
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   scale_color_manual(values=c("#506D8F","#E2A03C")) +
#   xlab("Patch Type") +
#   ylab(expression(paste("H2' (z-score)"))) +
#   theme_classic() +
#   ylim(c(2, 17)) +
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size = 30)) + # axis tick mark size
#   theme(axis.title = element_text(size = 34)) #+ # axis label size
# h2.pred

# 
# pdf(file = file.path("plots", "h2_noApis.pdf"), width = 8, height = 10)
# h2.pred
# dev.off()




#### interaction diversity ####
m.shannon.vaz <- glmmTMB(vaz.shannon ~ patch + (1|block), # no difference
                    data = network_vaznull)
summary(m.shannon.vaz)

m.shannon.vaz_noApis <- glmmTMB(vaz.shannon_noApis ~ patch + (1|block), # no difference
                           data = network_vaznull)
summary(m.shannon.vaz_noApis)

m.shannon.vaz_noPoe <- glmmTMB(vaz.shannon_noPoe ~ patch + (1|block), # no difference
                          data = network_vaznull)
summary(m.shannon.vaz_noPoe)


#### interaction evenness ####
m.evenness.vaz <- glmmTMB(vaz.evenness ~ patch + (1|block), # no difference
                         data = network_vaznull)
summary(m.evenness.vaz)

m.evenness.vaz_noApis <- glmmTMB(vaz.evenness_noApis ~ patch + (1|block), # no difference
                                data = network_vaznull)
summary(m.evenness.vaz_noApis)

m.evenness.vaz_noPoe <- glmmTMB(vaz.evenness_noPoe ~ patch + (1|block), # no difference
                               data = network_vaznull)
summary(m.evenness.vaz_noPoe)


#### pollinator: links per species ####
m.pol.links <- glmmTMB(pol.links ~ patch + (1|block), # lower in unconnected patches
                         data = network_vaznull)
summary(m.pol.links)

m.pol.links_noApis <- glmmTMB(pol.links_noApis ~ patch + (1|block), # lower in unconnected patches
                                data = network_vaznull)
summary(m.pol.links_noApis)

m.pol.links_noPoe <- glmmTMB(pol.links_noPoe ~ patch + (1|block), # lower in unconnected patches
                               data = network_vaznull)
summary(m.pol.links_noPoe)


#### plant: links per species ####
m.plant.links <- glmmTMB(plant.links ~ patch + (1|block), # lower in unconnected patches
                       data = network_vaznull)
summary(m.plant.links)

m.plant.links_noApis <- glmmTMB(plant.links_noApis ~ patch + (1|block), # lower in unconnected patches
                              data = network_vaznull)
summary(m.plant.links_noApis)

m.plant.links_noPoe <- glmmTMB(plant.links_noPoe ~ patch + (1|block), # lower in unconnected patches
                             data = network_vaznull)
summary(m.plant.links_noPoe)





#### DIVERSITY and ABUNDANCE ####
#### abundance ####
# abundance all interactions
m.abundance <- glmmTMB(abundance ~ patch + (1|block), # no difference
                       data = diversity_metrics, 
                       family = "nbinom2")
summary(m.abundance)

# abundance no apis
m.abundance_noApis <- glmmTMB(abundance_noApis ~ patch + (1|block), # significantly lower in unconnected
                       data = diversity_metrics, 
                       family = "nbinom2")
summary(m.abundance_noApis)

# abundance no poe
m.abundance_noPoe <- glmmTMB(abundance_noPoe ~ patch + (1|block), # no difference
                              data = diversity_metrics, 
                              family = "nbinom2")
summary(m.abundance_noPoe)

# abundance no poe or apis
m.abundance_noApis_Poe <- glmmTMB(abundance_noApis_Poe ~ patch + (1|block), # significantly lower in unconnected
                             data = diversity_metrics, 
                             family = "nbinom2")
summary(m.abundance_noApis_Poe)

# abundance butterflies
m.abundance_lep <- glmmTMB(abundance_lep ~ patch + (1|block), # marginally significantly lower in unconnected
                        data = diversity_metrics, 
                        family = "nbinom2")
summary(m.abundance_lep)

# abundance bee
m.abundance_bee <- glmmTMB(abundance_bee ~ patch + (1|block), # no difference
                           data = diversity_metrics, 
                           family = "nbinom2")
summary(m.abundance_bee)

# abundance fly
m.abundance_fly <- glmmTMB(abundance_fly ~ patch + (1|block), # no difference
                           data = diversity_metrics, 
                           family = "nbinom2")
summary(m.abundance_fly)


##### floral diversity #####
### ALTOGETHER - pollinators are interacting with a lower floral diversity in unconnected patches, but this effect goes away when taking out the most dominant fly species 
# shannon all species
m.floral_1 <- glmmTMB(floral_1 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_1)

# simpsons all species
m.floral_2 <- glmmTMB(floral_2 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_2)

# shannon no apis
m.floral_noApis_1 <- glmmTMB(floral_noApis_1 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_noApis_1)

# simpsons no apis
m.floral_noApis_2 <- glmmTMB(floral_noApis_2 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_noApis_2)

# shannon no poe
m.floral_noPoe_1 <- glmmTMB(floral_noPoe_1 ~ patch + (1|block), # no difference
                             data = diversity_metrics)
summary(m.floral_noPoe_1)

# simpsons no poe
m.floral_noPoe_2 <- glmmTMB(floral_noPoe_2 ~ patch + (1|block), # no difference
                             data = diversity_metrics)
summary(m.floral_noPoe_2)

# shannon no poe and apis
m.floral_noPoe_Apis_1 <- glmmTMB(floral_noPoe_Apis_1 ~ patch + (1|block), # no difference
                            data = diversity_metrics)
summary(m.floral_noPoe_Apis_1)

# simpsons no poe and apis
m.floral_noPoe_Apis_2 <- glmmTMB(floral_noPoe_Apis_2 ~ patch + (1|block), # no difference
                            data = diversity_metrics)
summary(m.floral_noPoe_Apis_2)





##### pollinator diversity #####
# lower pollinator diversity in unconnected patches no matter what metric
# shannon all species
m.pollinator_1 <- glmmTMB(pollinator_1 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.pollinator_1)

# simpsons all species
m.pollinator_2 <- glmmTMB(pollinator_2 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.pollinator_2)

# shannon no apis
m.pollinator_noApis_1 <- glmmTMB(pollinator_noApis_1 ~ patch + (1|block), # significantly lower in unconnected
                             data = diversity_metrics)
summary(m.pollinator_noApis_1)

# simpsons no apis
m.pollinator_noApis_2 <- glmmTMB(pollinator_noApis_2 ~ patch + (1|block), # significantly lower in unconnected
                             data = diversity_metrics)
summary(m.pollinator_noApis_2)

# shannon no poe
m.pollinator_noPoe_1 <- glmmTMB(pollinator_noPoe_1 ~ patch + (1|block), # significantly lower in unconnected
                            data = diversity_metrics)
summary(m.pollinator_noPoe_1)

# simpsons no poe
m.pollinator_noPoe_2 <- glmmTMB(pollinator_noPoe_2 ~ patch + (1|block), # significantly lower in unconnected
                            data = diversity_metrics)
summary(m.pollinator_noPoe_2)

# shannon no poe and apis
m.pollinator_noPoe_Apis_1 <- glmmTMB(pollinator_noPoe_Apis_1 ~ patch + (1|block), # significantly lower in unconnected
                                 data = diversity_metrics)
summary(m.pollinator_noPoe_Apis_1)
-2.3491/26.4289 * 100 # 8.9% decrease

# simpsons no poe and apis
m.pollinator_noPoe_Apis_2 <- glmmTMB(pollinator_noPoe_Apis_2 ~ patch + (1|block), # significantly lower in unconnected
                                 data = diversity_metrics)
summary(m.pollinator_noPoe_Apis_2)




##### order diversity ####
m.lep_1 <- glmmTMB(lep_1 ~ patch + (1|block), # no difference
                   data = diversity_metrics)
summary(m.lep_1)

m.lep_2 <- glmmTMB(lep_2 ~ patch + (1|block), # no difference
                   data = diversity_metrics)
summary(m.lep_2)

m.bee_1 <- glmmTMB(bee_1 ~ patch + (1|block), # marginally significantly lower in unconnected patches
                   data = diversity_metrics)
summary(m.bee_1)

m.bee_2 <- glmmTMB(bee_2 ~ patch + (1|block), # no difference
                   data = diversity_metrics)
summary(m.bee_2)

m.fly_1 <- glmmTMB(fly_1 ~ patch + (1|block), # no difference
                   data = diversity_metrics)
summary(m.fly_1)

m.fly_2 <- glmmTMB(fly_2 ~ patch + (1|block), # no difference
                   data = diversity_metrics)
summary(m.fly_2)

