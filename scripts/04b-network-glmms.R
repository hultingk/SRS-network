# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr)


# load data
network_vaznull <- read.csv(file = file.path("data", "L4_metrics", "network_vaznull.csv"))
modularity <- read.csv(file = file.path("data", "L4_metrics", "modularity.csv"))

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
# no difference in weighted connectance including or excluding apis and poe
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



#### modularity ####
m.module <- glmmTMB(vaz.module ~ patch + (1|block),
                    data = modularity)
summary(m.module)

m.module_noApis <- glmmTMB(vaz.module_noApis ~ patch + (1|block),
                           data = modularity)
summary(m.module_noApis)

m.module_noPoe <- glmmTMB(vaz.module_noPoe ~ patch + (1|block),
                          data = modularity)
summary(m.module_noPoe)