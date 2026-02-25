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

# # abundance butterflies
# m.abundance_lep <- glmmTMB(abundance_lep ~ patch + (1|block), # marginally significantly lower in unconnected
#                         data = diversity_metrics, 
#                         family = "nbinom2")
# summary(m.abundance_lep)
# 
# # abundance bee
# m.abundance_bee <- glmmTMB(abundance_bee ~ patch + (1|block), # no difference
#                            data = diversity_metrics, 
#                            family = "nbinom2")
# summary(m.abundance_bee)
# 
# # abundance fly
# m.abundance_fly <- glmmTMB(abundance_fly ~ patch + (1|block), # no difference
#                            data = diversity_metrics, 
#                            family = "nbinom2")
# summary(m.abundance_fly)


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


##### floral diversity plotting: shannon all species and shannon no Apis or Poe ####
# Shannon all species
# model predictions for plotting
m.floral_1.df <- diversity_metrics %>%
  dplyr::select(c("block", "patch", "floral_1"))
m.floral_1.df$floral_1_predict <- predict(m.floral_1, re.form = NA)
m.floral_1.df$patch <- factor(m.floral_1.df$patch, levels = c("B", "W"))
# plotting
floral_1.pred <- m.floral_1.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = floral_1, group = block), linewidth = 1.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = floral_1_predict, group = block), linewidth = 3.5, linetype = 1) +
  geom_point(aes(x = patch, y = floral_1, color = patch), size = 8, alpha = 0.55) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#DC267F","#FFB000")) +
  xlab("Patch type") +
  ylab(expression(paste("Floral diversity (Hill-Shannon)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") #+
  #theme(axis.text = element_text(size = 30)) + # axis tick mark size
  #theme(axis.title = element_text(size = 34)) #+ # axis label size
floral_1.pred

# Shannon excluding Apis and Poe interactions
# model predictions for plotting
m.floral_noPoe_Apis_1.df <- diversity_metrics %>%
  dplyr::select(c("block", "patch", "floral_noPoe_Apis_1"))
m.floral_noPoe_Apis_1.df$floral_noPoe_Apis_1_predict <- predict(m.floral_noPoe_Apis_1, re.form = NA)
m.floral_noPoe_Apis_1.df$patch <- factor(m.floral_noPoe_Apis_1.df$patch, levels = c("B", "W"))
# plotting
floral_noPoe_Apis_1.pred <- m.floral_noPoe_Apis_1.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = floral_noPoe_Apis_1, group = block), linewidth = 1.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = floral_noPoe_Apis_1_predict, group = block), linewidth = 3.5, linetype = 2) +
  geom_point(aes(x = patch, y = floral_noPoe_Apis_1, color = patch), size = 8, alpha = 0.55) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#DC267F","#FFB000")) +
  xlab("Patch type") +
  ylab(expression(atop("Floral diversity (Hill-Shannon)", paste("excluding 2 dominant pollinators")))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") #+
  #theme(axis.text = element_text(size = 30)) + # axis tick mark size
  #theme(axis.title = element_text(size = 34)) #+ # axis label size
floral_noPoe_Apis_1.pred

floral_div_plot <- cowplot::plot_grid(floral_1.pred, floral_noPoe_Apis_1.pred)
floral_div_plot



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
-2.2717/23.9314 * 100 # 9.49255% decrease

# simpsons no poe and apis
m.pollinator_noPoe_Apis_2 <- glmmTMB(pollinator_noPoe_Apis_2 ~ patch + (1|block), # significantly lower in unconnected
                                 data = diversity_metrics)
summary(m.pollinator_noPoe_Apis_2)


##### pollinator diversity plotting: shannon all species and shannon no Apis or Poe ####
# Shannon all species
# model predictions for plotting
m.pollinator_1.df <- diversity_metrics %>%
  dplyr::select(c("block", "patch", "pollinator_1"))
m.pollinator_1.df$pollinator_1_predict <- predict(m.pollinator_1, re.form = NA)
m.pollinator_1.df$patch <- factor(m.pollinator_1.df$patch, levels = c("B", "W"))
# plotting
pollinator_1.pred <- m.pollinator_1.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = pollinator_1, group = block), linewidth = 1.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = pollinator_1_predict, group = block), linewidth = 3.5, linetype = 1) +
  geom_point(aes(x = patch, y = pollinator_1, color = patch), size = 8, alpha = 0.55) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#DC267F","#FFB000")) +
  xlab("Patch type") +
  ylab(expression(paste("Pollinator diversity (Hill-Shannon)"))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") #+
#theme(axis.text = element_text(size = 30)) + # axis tick mark size
#theme(axis.title = element_text(size = 34)) #+ # axis label size
pollinator_1.pred

# Shannon excluding Apis and Poe interactions
# model predictions for plotting
m.pollinator_noPoe_Apis_1.df <- diversity_metrics %>%
  dplyr::select(c("block", "patch", "pollinator_noPoe_Apis_1"))
m.pollinator_noPoe_Apis_1.df$pollinator_noPoe_Apis_1_predict <- predict(m.pollinator_noPoe_Apis_1, re.form = NA)
m.pollinator_noPoe_Apis_1.df$patch <- factor(m.pollinator_noPoe_Apis_1.df$patch, levels = c("B", "W"))
# plotting
pollinator_noPoe_Apis_1.pred <- m.pollinator_noPoe_Apis_1.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = pollinator_noPoe_Apis_1, group = block), linewidth = 1.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = pollinator_noPoe_Apis_1_predict, group = block), linewidth = 3.5, linetype = 1) +
  geom_point(aes(x = patch, y = pollinator_noPoe_Apis_1, color = patch), size = 8, alpha = 0.55) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#DC267F","#FFB000")) +
  xlab("Patch type") +
  ylab(expression(atop("Pollinator diversity (Hill-Shannon)", paste("excluding 2 dominant pollinators")))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") #+
#theme(axis.text = element_text(size = 30)) + # axis tick mark size
#theme(axis.title = element_text(size = 34)) #+ # axis label size
pollinator_noPoe_Apis_1.pred

pollinator_div_plot <- cowplot::plot_grid(pollinator_1.pred, pollinator_noPoe_Apis_1.pred)
pollinator_div_plot






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
