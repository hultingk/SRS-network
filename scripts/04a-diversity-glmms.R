# -------------------------------------- #
#### Analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries 
librarian::shelf(tidyverse, glmmTMB, DHARMa, performance, ggeffects, ggpubr)


# load data
#network_metrics <- read.csv(file = file.path("data", "network_metrics.csv"))
diversity_metrics <- read.csv(file = file.path("data", "L4_metrics", "diversity_metrics.csv"))


#### DIVERSITY and ABUNDANCE ####
#### abundance ####
# abundance all interactions
m.abundance <- glmmTMB(abundance ~ patch + (1|block), # no difference
                       data = diversity_metrics, 
                       family = "nbinom2")
summary(m.abundance)

# abundance no poe or apis
m.abundance_noApis_Poe <- glmmTMB(abundance_noApis_Poe ~ patch + (1|block), # significantly lower in unconnected
                             data = diversity_metrics, 
                             family = "nbinom2")
summary(m.abundance_noApis_Poe)




##### floral diversity #####
### ALTOGETHER - pollinators are interacting with a lower floral diversity in unconnected patches, but this effect goes away when taking out the most dominant fly species 
# shannon all species
m.floral_1 <- glmmTMB(floral_1 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_1)
-4.486/11.763 * 100 # 38.13653% decrease

# simpsons all species
m.floral_2 <- glmmTMB(floral_2 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.floral_2)


# shannon no poe and apis
m.floral_noPoe_Apis_1 <- glmmTMB(floral_noPoe_Apis_1 ~ patch + (1|block), # no difference
                            data = diversity_metrics)
summary(m.floral_noPoe_Apis_1)
-1.400/10.106 * 100 # 13.85316% decrease, but not significant

# simpsons no poe and apis
m.floral_noPoe_Apis_2 <- glmmTMB(floral_noPoe_Apis_2 ~ patch + (1|block), # no difference
                            data = diversity_metrics)
summary(m.floral_noPoe_Apis_2)


##### floral diversity plotting: shannon all species and shannon no Apis or Poe ####
# Shannon all species
# model predictions for plotting
m.floral_1.df <- ggpredict(m.floral_1, terms = c("patch"), back_transform = TRUE)
# plotting
floral_1.pred <- m.floral_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = floral_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.floral_1.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression("Floral diversity (Hill-Shannon)")) +
  theme_classic(base_size = 20) +
  ylim(2, 16.5) +
  theme(legend.position = "none") 
floral_1.pred

# Shannon excluding Apis and Poe interactions
# model predictions for plotting
m.floral_noPoe_Apis_1.df <- ggpredict(m.floral_noPoe_Apis_1, terms = c("patch"), back_transform = TRUE)
# plotting
floral_1.pred_noPoe_Apis <- m.floral_noPoe_Apis_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = floral_noPoe_Apis_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.floral_noPoe_Apis_1.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(atop("Floral diversity (Hill-Shannon)", paste("excluding two dominant pollinators")))) +
  theme_classic(base_size = 20) +
  ylim(2, 16.5) +
  theme(legend.position = "none") 
floral_1.pred_noPoe_Apis


# all together
floral_div_plot <- cowplot::plot_grid(floral_1.pred, floral_1.pred_noPoe_Apis)
floral_div_plot

# exporting
# pdf(file = file.path("plots", "floral_div_plot.pdf"), width = 13, height = 6.5)
# floral_div_plot
# dev.off()



##### pollinator diversity #####
# lower pollinator diversity in unconnected patches no matter what metric
# shannon all species
m.pollinator_1 <- glmmTMB(pollinator_1 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.pollinator_1)
-11.290/27.638 * 100 # 40.84955% decrease

# simpsons all species
m.pollinator_2 <- glmmTMB(pollinator_2 ~ patch + (1|block), # significantly lower in unconnected
                      data = diversity_metrics)
summary(m.pollinator_2)

# shannon no poe and apis
m.pollinator_noPoe_Apis_1 <- glmmTMB(pollinator_noPoe_Apis_1 ~ patch + (1|block), # significantly lower in unconnected
                                 data = diversity_metrics)
summary(m.pollinator_noPoe_Apis_1)
-2.2717/23.9314 * 100 # 9.49255% decrease

# simpsons no poe and apis
m.pollinator_noPoe_Apis_2 <- glmmTMB(pollinator_noPoe_Apis_2 ~ patch + (1|block), # significantly lower in unconnected
                                 data = diversity_metrics)
summary(m.pollinator_noPoe_Apis_2)
-1.5995/18.3849 * 100 # 8.700075% decrease


##### pollinator diversity plotting: shannon all species and shannon no Apis or Poe ####
# Shannon all species
# model predictions for plotting
# Shannon all species
# model predictions for plotting
m.pollinator_1.df <- ggpredict(m.pollinator_1, terms = c("patch"), back_transform = TRUE)
# plotting
pollinator_1.pred <- m.pollinator_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = pollinator_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.pollinator_1.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression("Pollinator diversity (Hill-Shannon)")) +
  theme_classic(base_size = 20) +
  ylim(10, 35) +
  theme(legend.position = "none") 
pollinator_1.pred

# shannon excluding 2 dominant pollinators
# model predictions for plotting
m.pollinator_noPoe_Apis_1.df <- ggpredict(m.pollinator_noPoe_Apis_1, terms = c("patch"), back_transform = TRUE)
# plotting
pollinator_1.pred_noPoe_Apis <- m.pollinator_noPoe_Apis_1.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = pollinator_noPoe_Apis_1, color = patch), data = diversity_metrics, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.pollinator_noPoe_Apis_1.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  xlab("Patch type") +
  ylab(expression(atop("Pollinator diversity (Hill-Shannon)", paste("excluding two dominant pollinators")))) +
  theme_classic(base_size = 20) +
  ylim(10, 35) +
  theme(legend.position = "none") 
pollinator_1.pred_noPoe_Apis


pollinator_div_plot <- cowplot::plot_grid(pollinator_1.pred, pollinator_1.pred_noPoe_Apis)
pollinator_div_plot


# exporting
# pdf(file = file.path("plots", "pollinator_div_plot.pdf"), width = 13, height = 6.5)
# pollinator_div_plot
# dev.off()



# ##### order diversity ####
# m.lep_1 <- glmmTMB(lep_1 ~ patch + (1|block), # no difference
#                    data = diversity_metrics)
# summary(m.lep_1)
# 
# m.lep_2 <- glmmTMB(lep_2 ~ patch + (1|block), # no difference
#                    data = diversity_metrics)
# summary(m.lep_2)
# 
# m.bee_1 <- glmmTMB(bee_1 ~ patch + (1|block), # marginally significantly lower in unconnected patches
#                    data = diversity_metrics)
# summary(m.bee_1)
# 
# m.bee_2 <- glmmTMB(bee_2 ~ patch + (1|block), # no difference
#                    data = diversity_metrics)
# summary(m.bee_2)
# 
# m.fly_1 <- glmmTMB(fly_1 ~ patch + (1|block), # no difference
#                    data = diversity_metrics)
# summary(m.fly_1)
# 
# m.fly_2 <- glmmTMB(fly_2 ~ patch + (1|block), # no difference
#                    data = diversity_metrics)
# summary(m.fly_2)
# 
