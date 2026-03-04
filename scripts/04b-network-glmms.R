# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr)


# load data
network_vaznull <- read.csv(file = file.path("data", "L4_metrics", "network_vaznull.csv"))
modularity <- read.csv(file = file.path("data", "L4_metrics", "modularity.csv"))
extinction <- read.csv(file = file.path("data", "L4_metrics", "extinction_metric.csv"))

#### connectance ####
# no vaznull metric for connectance
# unconnected patches have higher connectance, even excluding Apis and Poe
m.connect <- glmmTMB(connectance ~ patch + (1|block),
                     data = network_vaznull)
summary(m.connect)

m.connect_noPoe_Apis <- glmmTMB(connectance_noPoe_Apis ~ patch + (1|block),
                           data = network_vaznull)
summary(m.connect_noPoe_Apis)



#### weighted connectance ####
# no vaznull metric for weighted connectance
# no difference in weighted connectance including or excluding apis and poe, except when both are excluded
m.weightconnect <- glmmTMB(weightconnectance ~ patch + (1|block),
                           data = network_vaznull)
summary(m.weightconnect)

m.weightconnect_noPoe_Apis <- glmmTMB(weightconnectance_noPoe_Apis ~ patch + (1|block),
                                 data = network_vaznull)
summary(m.weightconnect_noPoe_Apis)



#### nestedness ####
m.nestedness.vaz <- glmmTMB(vaz.NODF ~ patch + (1|block), # significantly lower in connected patches 
                            data = network_vaznull)
summary(m.nestedness.vaz)

m.nestedness.vaz_noPoe_Apis <- glmmTMB(vaz.NODF_noPoe_Apis ~ patch + (1|block), # no difference
                                  data = network_vaznull)
summary(m.nestedness.vaz_noPoe_Apis)

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

m.h2.vaz_noPoe_Apis <- glmmTMB(vaz.h2_noPoe_Apis ~ patch + (1|block), # no difference
                          data = network_vaznull)
summary(m.h2.vaz_noPoe_Apis)


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
  geom_jitter(aes(x = patch, y = links.per.sp, color = patch), data = network_vaznull, size = 5, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.links_predict, width = 0, linewidth = 2.5) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("Links per species"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.links_predict_plot


# model predictions for plotting
m.links.df <- network_vaznull %>%
  dplyr::select(c("block", "patch", "links.per.sp"))
m.links.df$links.pred <- predict(m.links, re.form = NA)
m.links.df$patch <- factor(m.links.df$patch, levels = c("B", "W"))
# plotting
links.pred <- m.links.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = links.per.sp, group = block), linewidth = 2.5, color = "gray40", alpha = 0.5) +
  geom_line(aes(x = patch, y = links.pred, group = block), linewidth = 5, linetype = 1) +
  geom_point(aes(x = patch, y = links.per.sp, color = patch), size = 11, alpha = 0.8) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch Type") +
  ylab(expression(paste("Links per species"))) +
  theme_classic() +
  ylim(1.2, 1.8) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
links.pred


# links per species plotting - exclusing 2 dominant species
# model predictions for plotting
m.links.df_noPoe_Apis <- network_vaznull %>%
  dplyr::select(c("block", "patch", "links.per.sp_noPoe_Apis"))
m.links.df_noPoe_Apis$links.pred_noPoe_Apis <- predict(m.links_noPoe_Apis, re.form = NA)
m.links.df_noPoe_Apis$patch <- factor(m.links.df_noPoe_Apis$patch, levels = c("B", "W"))
# plotting
links.pred_noPoe_Apis <- m.links.df_noPoe_Apis%>%
  ggplot() +
  geom_line(aes(x = patch, y = links.per.sp_noPoe_Apis, group = block), linewidth = 2.5, color = "gray40", alpha = 0.5) +
  geom_line(aes(x = patch, y = links.pred_noPoe_Apis, group = block), linewidth = 5, linetype = 1) +
  geom_point(aes(x = patch, y = links.per.sp_noPoe_Apis, color = patch), size = 11, alpha = 0.8) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch Type") +
  ylab(expression(paste("Links per species"))) +
  theme_classic() +
  ylim(1.2, 1.8) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
links.pred_noPoe_Apis


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




#### robustness to extinction - Lower level ####
m.extinct.LL <- glmmTMB(extinct.robustness.LL ~ patch + (1|block),
                     data = extinction)
summary(m.extinct.LL)

m.extinct.LL_noPoe_Apis <- glmmTMB(extinct.robustness.LL_noPoe_Apis ~ patch + (1|block),
                     data = extinction)
summary(m.extinct.LL_noPoe_Apis)

# all species - plotting robustness to secondary extinction - robustness to plant extinctions 
# model predictions for plotting
m.extinct.LL.predict <- ggpredict(m.extinct.LL, terms = c("patch"), back_transform = T)
# plotting
m.extinct.LL_plot <- m.extinct.LL.predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = extinct.robustness.LL, color = patch), data = extinction, size = 5, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.extinct.LL.predict, width = 0, linewidth = 2.5) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("Robustness to plant extinctions (z-score)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.extinct.LL_plot

# excluding 2 dominant pollinators - plotting robustness to secondary extinction - robustness to plant extinctions 
# model predictions for plotting
m.extinct.LL.predict_noPoe_Apis <- ggpredict(m.extinct.LL_noPoe_Apis, terms = c("patch"), back_transform = T)
# plotting
m.extinct.LL_plot_noPoe_Apis <- m.extinct.LL.predict_noPoe_Apis %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = extinct.robustness.LL_noPoe_Apis, color = patch), data = extinction, size = 5, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.extinct.LL.predict_noPoe_Apis, width = 0, linewidth = 2.5) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("Robustness to plant extinctions excluding 2 dominant pollinators (z-score)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.extinct.LL_plot_noPoe_Apis


#### robustness to extinction - higher level ####
m.extinct.HL <- glmmTMB(extinct.robustness.HL ~ patch + (1|block),
                        data = extinction)
summary(m.extinct.HL)

m.extinct.HL_noPoe_Apis <- glmmTMB(extinct.robustness.HL_noPoe_Apis ~ patch + (1|block),
                                   data = extinction)
summary(m.extinct.HL_noPoe_Apis)


# all species - plotting robustness to secondary extinction - robustness to pollinator extinctions 
# model predictions for plotting
m.extinct.HL.predict <- ggpredict(m.extinct.HL, terms = c("patch"), back_transform = T)
# plotting
m.extinct.HL_plot <- m.extinct.HL.predict %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = extinct.robustness.HL, color = patch), data = extinction, size = 5, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.extinct.HL.predict, width = 0, linewidth = 2.5) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 6, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(paste("Robustness to plant extinctions (z-score)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") 
m.extinct.HL_plot
