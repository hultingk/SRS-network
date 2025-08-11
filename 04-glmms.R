# -------------------------------------- #
#### Analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries 
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggeffects)
library(ggpubr)
library(tidyverse)

# load data
network_metrics <- read.csv(file = file.path("data", "network_metrics.csv"))


#### nestedness ####
m.nestedness <- glmmTMB(NODF ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.nestedness)
plot(simulateResiduals(m.nestedness))
check_model(m.nestedness)

m.nestedness.noApis <- glmmTMB(NODF_noApis ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.nestedness.noApis)

# plotting NO APIS
# model predictions for plotting
m.nestedness.noApis.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "NODF_noApis"))
m.nestedness.noApis.df$nestedness_pred <- predict(m.nestedness.noApis, re.form = NA)
m.nestedness.noApis.df$patch <- factor(m.nestedness.noApis.df$patch, levels = c("B", "W"))
# plotting
nestedness.pred <- m.nestedness.noApis.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = NODF_noApis, group = block), linewidth = 2.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = nestedness_pred, group = block), linewidth = 5.5, linetype = 3) +
  geom_point(aes(x = patch, y = NODF_noApis, color = patch), size = 11, alpha = 0.95) + 
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("NODF (z-score)"))) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
nestedness.pred

pdf(file = file.path("plots", "nestedness_noApis.pdf"), width = 8, height = 10)
nestedness.pred
dev.off()




#### H2' ####
m.h2 <- glmmTMB(H2 ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.h2)
plot(simulateResiduals(m.h2))
check_model(m.h2)

m.h2.noApis <- glmmTMB(h2_noApis ~ patch + (1|block),
                data = network_metrics,
                family = "gaussian")
summary(m.h2.noApis)
(2.1497/5.4302)*100 # 39.58786 less nested
# plotting
# model predictions for plotting
m.h2.noApis.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "h2_noApis"))
m.h2.noApis.df$h2.pred <- predict(m.h2.noApis, re.form = NA)
m.h2.noApis.df$patch <- factor(m.h2.noApis.df$patch, levels = c("B", "W"))
# plotting
h2.pred <- m.h2.noApis.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = h2_noApis, group = block), linewidth = 2.5, color = "lightslategrey", alpha = 0.5) +
  geom_line(aes(x = patch, y = h2.pred, group = block), linewidth = 5, linetype = 1) +
  geom_point(aes(x = patch, y = h2_noApis, color = patch), size = 11, alpha = 0.8) + 
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("H2' (z-score)"))) +
  theme_classic() +
  ylim(c(2, 17)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
h2.pred


pdf(file = file.path("plots", "h2_noApis.pdf"), width = 8, height = 10)
h2.pred
dev.off()




#### interaction diversity ####
m.interaction.div <- glmmTMB(Shannon.diversity ~ patch + (1|block),
                data = network_metrics,
                family = "gaussian")
summary(m.interaction.div)
plot(simulateResiduals(m.interaction.div))
check_model(m.interaction.div)

m.interaction.div.noApis <- glmmTMB(shannon_noApis ~ patch + (1|block),
                             data = network_metrics,
                             family = "gaussian")
summary(m.interaction.div.noApis)
(-5.633/-4.397)*100
#### links per species ####
m.links <- glmmTMB(links.per.species ~ patch + (1|block),
                             data = network_metrics,
                             family = "gaussian")
summary(m.links)
plot(simulateResiduals(m.links))
check_model(m.links)

m.links.noApis <- glmmTMB(links_noApis ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.links.noApis)


#### linkage density ####
m.density <- glmmTMB(linkage.density ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.density)
plot(simulateResiduals(m.density))
check_model(m.density)

m.density.noApis <- glmmTMB(density_noApis ~ patch + (1|block),
                     data = network_metrics,
                     family = "gaussian")
summary(m.density.noApis)
(-0.8167/4.8001)*100 #17% decrease in linkage density
# plotting
# model predictions for plotting
m.density.noApis.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "density_noApis"))
m.density.noApis.df$linkage.density.pred <- predict(m.density.noApis, re.form = NA)
m.density.noApis.df$patch <- factor(m.density.noApis.df$patch, levels = c("B", "W"))
# plotting
linkage.density.pred <- m.density.noApis.df %>%
  ggplot() +
  geom_line(aes(x = patch, y = density_noApis, group = block), linewidth = 2.5, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = linkage.density.pred, group = block), linewidth = 5.5) +
  geom_point(aes(x = patch, y = density_noApis, color = patch), size = 11, alpha = 0.95) + 
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("Linkage density (z-score)"))) +
  theme_classic() +
  ylim(c(2.5, 5.5)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
linkage.density.pred

pdf(file = file.path("plots", "linkage.density_noApis.pdf"), width = 8, height = 10)
linkage.density.pred
dev.off()




#### flower diversity ####
# with Apis
m.floral.div <- glmmTMB(floral_diversity ~ patch + (1|block), # diversity
                   data = network_metrics,
                   family = "gaussian")
m.floral.div <- glmmTMB(fl.rich ~ patch + (1|block), # richness
                        data = network_metrics,
                        family = "nbinom2")
summary(m.floral.div)
plot(simulateResiduals(m.floral.div))
check_model(m.floral.div)

exp(3.21803)
exp(3.21803-0.20816)

# no Apis
m.floral.div.noApis <- glmmTMB(floral_div_noApis ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
m.floral.div.noApis <- glmmTMB(fl.rich.no_Apis ~ patch + (1|block), # richness
                        data = network_metrics,
                        family = "nbinom2")
summary(m.floral.div.noApis)
plot(simulateResiduals(m.floral.div.noApis))
check_model(m.floral.div.noApis)
(-3.698/21.616)*100
(-0.36451/2.51352)*100

exp(3.21770)
exp(3.21770-0.21473)

# plotting
# model predictions for plotting
m.floral.div.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "fl.rich.no_Apis"))
m.floral.div.df$floral.div.pred <- predict(m.floral.div.noApis, re.form = NA)
m.floral.div.df$patch <- factor(m.floral.div.df$patch, levels = c("B", "W"))
# plotting
floral.div.pred <- m.floral.div.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = fl.rich.no_Apis, color = patch), size = 9, alpha = 0.95) + 
  geom_line(aes(x = patch, y = fl.rich.no_Apis, group = block), linewidth = 2, color = "grey25", alpha = 0.85) +
  geom_line(aes(x = patch, y = exp(floral.div.pred), group = block), linewidth = 5.5) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("Floral richness"))) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
floral.div.pred

pdf(file = file.path("plots", "floral_diversity.pdf"), width = 8, height = 10)
floral.div.pred
dev.off()


# plotting
# model predictions for plotting
m.floral.div.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "floral_diversity"))
m.floral.div.df$floral.div.pred <- predict(m.floral.div, re.form = NA)
m.floral.div.df$patch <- factor(m.floral.div.df$patch, levels = c("B", "W"))
# plotting
floral.div.pred <- m.floral.div.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = floral_diversity, color = patch), size = 9, alpha = 0.7) + 
  geom_line(aes(x = patch, y = floral_diversity, group = block), linewidth = 2, color = "black", alpha = 0.2) +
  geom_line(aes(x = patch, y = floral.div.pred, group = block), linewidth = 4) +
  scale_x_discrete(labels = c('Connected', 'Winged')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("Floral diversity (shannon)"))) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
floral.div.pred

pdf(file = file.path("plots", "floral_diversity.pdf"), width = 8, height = 10)
floral.div.pred
dev.off()












#### pollinator diversity ####
m.pollinator.div <- glmmTMB(pollinator_diversity ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.pollinator.div)
plot(simulateResiduals(m.pollinator.div))
check_model(m.pollinator.div)

# pollinator richness
m.pollinator.rich <- glmmTMB(pollinator.rich ~ patch + (1|block),
                            data = network_metrics,
                            family = "poisson")
summary(m.pollinator.rich)

# no Apis
m.pollinator.div.noApis <- glmmTMB(pollinator_div_noApis ~ patch + (1|block),
                            data = network_metrics,
                            family = "gaussian")
summary(m.pollinator.div.noApis)
plot(simulateResiduals(m.pollinator.div.noApis))
check_model(m.pollinator.div.noApis)


# plotting
# model predictions for plotting
m.pollinator.div.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "pollinator_div_noApis"))
m.pollinator.div.df$pollinator_diversity.pred <- predict(m.pollinator.div, re.form = NA)
m.pollinator.div.df$patch <- factor(m.pollinator.div.df$patch, levels = c("B", "W"))
# plotting
pollinator_diversity.pred <- m.pollinator.div.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = pollinator_div_noApis, color = patch), size = 9, alpha = 0.7) + 
  geom_line(aes(x = patch, y = pollinator_div_noApis, group = block), linewidth = 2, color = "black", alpha = 0.2) +
  geom_line(aes(x = patch, y = pollinator_diversity.pred, group = block), linewidth = 4) +
  scale_x_discrete(labels = c('Connected', 'Winged')) +
  scale_color_manual(values=c("#506D8F","#E2A03C")) +
  xlab("Patch Type") +
  ylab(expression(paste("Pollinator diversity (shannon)"))) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 30)) + # axis tick mark size
  theme(axis.title = element_text(size = 34)) #+ # axis label size
pollinator_diversity.pred

pdf(file = file.path("plots", "pollinator_diversity.pdf"), width = 8, height = 10)
pollinator_diversity.pred
dev.off()










