# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr)


# load data
network_vaznull <- read.csv(file = file.path("data", "L4_metrics", "network_vaznull.csv"))

# #### connectance ####
# # no vaznull metric for connectance
# # unconnected patches have higher connectance, even excluding Apis and Poe
# m.connect <- glmmTMB(connectance ~ patch + (1|block),
#                      data = network_vaznull)
# summary(m.connect)
# 
# m.connect_noPoe_Apis <- glmmTMB(connectance_noPoe_Apis ~ patch + (1|block),
#                            data = network_vaznull)
# summary(m.connect_noPoe_Apis)
# 
# 
# 
# #### weighted connectance ####
# # no vaznull metric for weighted connectance
# # no difference in weighted connectance including or excluding apis and poe, except when both are excluded
# m.weightconnect <- glmmTMB(weightconnectance ~ patch + (1|block),
#                            data = network_vaznull)
# summary(m.weightconnect)
# 
# m.weightconnect_noPoe_Apis <- glmmTMB(weightconnectance_noPoe_Apis ~ patch + (1|block),
#                                  data = network_vaznull)
# summary(m.weightconnect_noPoe_Apis)



#### nestedness ####
m.nestedness.vaz <- glmmTMB(vaz.NODF ~ patch + (1|block), # significantly lower in connected patches 
                            data = network_vaznull)
summary(m.nestedness.vaz)

m.nestedness.vaz_noPoe_Apis <- glmmTMB(vaz.NODF_noPoe_Apis ~ patch + (1|block), # no difference
                                  data = network_vaznull)
summary(m.nestedness.vaz_noPoe_Apis)

# plotting
# model predictions for plotting - all species
m.nest_predict <- ggpredict(m.nestedness.vaz, terms = c("patch"), back_transform = T)
# plotting
m.nest_predict_plot <- m.nest_predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = vaz.NODF, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.nest_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("NODF (z-score)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.nest_predict_plot

# excluding 2 dominant
# model predictions for plotting 
m.nest_predict_noPoe_Apis <- ggpredict(m.nestedness.vaz_noPoe_Apis, terms = c("patch"), back_transform = T)
# plotting
m.nest_predict_plot_noPoe_Apis <- m.nest_predict_noPoe_Apis %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = vaz.NODF_noPoe_Apis, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.nest_predict_noPoe_Apis, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(atop("NODF (z-score)", paste("excluding two dominant pollinators")))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.nest_predict_plot_noPoe_Apis

# all together
nest.plot <- cowplot::plot_grid(m.nest_predict_plot, m.nest_predict_plot_noPoe_Apis, rel_widths = c(1, 1.1))
nest.plot

# pdf(file = file.path("plots", "nest.plot.pdf"), width = 13, height = 6.5)
# nest.plot
# dev.off()



#### H2' ####
m.h2.vaz <- glmmTMB(vaz.h2 ~ patch + (1|block), # significantly lower in connected patches 
                    data = network_vaznull)
summary(m.h2.vaz)

m.h2.vaz_noPoe_Apis <- glmmTMB(vaz.h2_noPoe_Apis ~ patch + (1|block), # no difference
                          data = network_vaznull)
summary(m.h2.vaz_noPoe_Apis)


# model predictions for plotting - all species
m.h2_predict <- ggpredict(m.h2.vaz, terms = c("patch"), back_transform = T)
# plotting
m.h2_predict_plot <- m.h2_predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = vaz.h2, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.h2_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("H2' (z-score)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.h2_predict_plot


# excluding 2 dominant
# model predictions for plotting 
m.h2_predict_noPoe_Apis <- ggpredict(m.h2.vaz_noPoe_Apis, terms = c("patch"), back_transform = T)
# plotting
m.h2_predict_plot_noPoe_Apis <- m.h2_predict_noPoe_Apis %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = vaz.h2_noPoe_Apis, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.h2_predict_noPoe_Apis, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(atop("H2' (z-score)", paste("excluding two dominant pollinators")))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.h2_predict_plot_noPoe_Apis

# all together
h2.plot <- cowplot::plot_grid(m.h2_predict_plot, m.h2_predict_plot_noPoe_Apis, rel_widths = c(1, 1.1))
h2.plot

# pdf(file = file.path("plots", "h2.plot.pdf"), width = 13, height = 6.5)
# h2.plot
# dev.off()








### links per species ###
m.links <- glmmTMB(links.per.sp ~ patch + (1|block), 
                   data = network_vaznull)
summary(m.links)

m.links_noPoe_Apis <- glmmTMB(links.per.sp_noPoe_Apis ~ patch + (1|block), 
                   data = network_vaznull)
summary(m.links_noPoe_Apis)

# links per species plotting - all species
# model predictions for plotting
m.links_predict <- ggpredict(m.links, terms = c("patch"), back_transform = T)
# plotting
m.links_predict_plot <- m.links_predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = links.per.sp, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.links_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("Links per species"))) +
  theme_classic(base_size = 20) +
  ylim(1.2, 1.7) +
  theme(legend.position = "none") 
m.links_predict_plot


# links per species plotting - no apis or poe
# model predictions for plotting
m.links_predict_noPoe_Apis <- ggpredict(m.links_noPoe_Apis, terms = c("patch"), back_transform = T)
# plotting
m.links_predict_plot_noPoe_Apis <- m.links_predict_noPoe_Apis %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = links.per.sp_noPoe_Apis, color = patch), data = network_vaznull, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.links_predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(atop("Links per species", paste("excluding two dominant pollinators")))) +
  theme_classic(base_size = 20) +
  ylim(1.2, 1.7) +
  theme(legend.position = "none") 
m.links_predict_plot_noPoe_Apis


# all together
links_predict_plot <- cowplot::plot_grid(m.links_predict_plot, m.links_predict_plot_noPoe_Apis, rel_widths = c(1, 1.1))
links_predict_plot

# pdf(file = file.path("plots", "links_predict_plot.pdf"), width = 13, height = 6.5)
# links_predict_plot
# dev.off()



#### linkage density ####
m.density <- glmmTMB(vaz.density ~ patch + (1|block), 
                   data = network_vaznull)
summary(m.density)

m.density_noPoe_Apis <- glmmTMB(vaz.density_noPoe_Apis ~ patch + (1|block), 
                              data = network_vaznull)
summary(m.density_noPoe_Apis)


#### web asymmetry ####
m.asymmetry <- glmmTMB(web.asymmetry ~ patch + (1|block), 
                     data = network_vaznull)
summary(m.asymmetry)

m.asymmetry_noPoe_Apis <- glmmTMB(web.asymmetry_noPoe_Apis ~ patch + (1|block), 
                                data = network_vaznull)
summary(m.asymmetry_noPoe_Apis)





# #### pollinator: mean links ####
# m.pol.links <- glmmTMB(pol.links ~ patch + (1|block), # lower in unconnected patches
#                        data = network_vaznull)
# summary(m.pol.links)
# 
# m.pol.links_noPoe_Apis <- glmmTMB(pol.links_noPoe_Apis ~ patch + (1|block), 
#                              data = network_vaznull)
# summary(m.pol.links_noPoe_Apis)
# 
# # # plotting
# # # links per species pollinator all species
# # # model predictions for plotting
# # m.pol.links_predict<- ggpredict(m.pol.links, terms = c("patch"), back_transform = T)
# # # plotting
# # pol.links.pred <- m.pol.links_predict %>%
# #   ggplot() +
# #   geom_jitter(aes(x = patch, y = pol.links, color = patch), data = network_vaznull, size = 5, alpha = 0.55,
# #               width = 0.08, height = 0) +
# #   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
# #                 data = m.pol.links_predict, width = 0, linewidth = 2.5) +
# #   geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
# #   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
# #   scale_color_manual(values=c("#F5097C","#F7B3D4")) +
# #   scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
# #   xlab("Patch type") +
# #   ylab(expression(paste("Links per pollinator species"))) +
# #   theme_classic(base_size = 20) +
# #   ylim(1.2, 1.75) +
# #   theme(legend.position = "none") 
# # 
# # # links per species pollinator excluding Apis and Poe
# # # model predictions for plotting
# # m.pol.links_noPoe_Apis_predict<- ggpredict(m.pol.links_noPoe_Apis, terms = c("patch"), back_transform = T)
# # # plotting
# # pol.links_noPoe_Apis.pred <- m.pol.links_noPoe_Apis_predict %>%
# #   ggplot() +
# #   geom_jitter(aes(x = patch, y = pol.links_noPoe_Apis, color = patch), data = network_vaznull, size = 5, alpha = 0.55,
# #               width = 0.08, height = 0) +
# #   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
# #                 data = m.pol.links_noPoe_Apis_predict, width = 0, linewidth = 2.5) +
# #   geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
# #   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
# #   scale_color_manual(values=c("#F5097C","#F7B3D4")) +
# #   scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
# #   xlab("Patch type") +
# #   ylab(expression(paste("Links per pollinator species (exclusing 2 dominant pollinators)"))) +
# #   theme_classic(base_size = 20) +
# #   ylim(1.2, 1.75) +
# #   theme(legend.position = "none") 
# # 
# # pol.links_all_plot <- cowplot::plot_grid(pol.links.pred, pol.links_noPoe_Apis.pred)
# # pol.links_all_plot
# 
# 
# 
# 
# 
# #### plant: links per species ####
# m.plant.links <- glmmTMB(plant.links ~ patch + (1|block), # lower in unconnected patches
#                          data = network_vaznull)
# summary(m.plant.links)
# 
# m.plant.links_noPoe_Apis <- glmmTMB(plant.links_noPoe_Apis ~ patch + (1|block), # higher in unconnected patches
#                                data = network_vaznull)
# summary(m.plant.links_noPoe_Apis)
# 
# 
# # plotting
# # links per species plant all species
# # model predictions for plotting
# m.plant.links_predict<- ggpredict(m.plant.links, terms = c("patch"), back_transform = T)
# # plotting
# plant.links.pred <- m.plant.links_predict %>%
#   ggplot() +
#   geom_jitter(aes(x = patch, y = plant.links, color = patch), data = network_vaznull, size = 5, alpha = 0.55,
#               width = 0.08, height = 0) +
#   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
#                 data = m.plant.links_predict, width = 0, linewidth = 2.5) +
#   geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   scale_color_manual(values=c("#F5097C","#F7B3D4")) +
#   scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
#   xlab("Patch type") +
#   ylab(expression(paste("Mean number of links (z-score)"))) +
#   theme_classic(base_size = 20) +
#   theme(legend.position = "none") 
# plant.links.pred
# 
# # links per plant species excluding dominant pollinators
# # model predictions for plotting
# m.plant.links_predict_noPoe_Apis <- ggpredict(m.plant.links_noPoe_Apis, terms = c("patch"), back_transform = T)
# # plotting
# plant.links.pred_noPoe_Apis <- m.plant.links_predict_noPoe_Apis %>%
#   ggplot() +
#   geom_jitter(aes(x = patch, y = plant.links_noPoe_Apis, color = patch), data = network_vaznull, size = 5, alpha = 0.55,
#               width = 0.08, height = 0) +
#   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
#                 data = m.plant.links_predict_noPoe_Apis, width = 0, linewidth = 2.5) +
#   geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   scale_color_manual(values=c("#F5097C","#F7B3D4")) +
#   scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
#   xlab("Patch type") +
#   ylab(expression(paste("Links per plant species (z-score) excluding 2 dominant pollinators"))) +
#   theme_classic(base_size = 20) +
#   theme(legend.position = "none") 
# plant.links.pred_noPoe_Apis





