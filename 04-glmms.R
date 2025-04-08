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

m.nestedness <- glmmTMB(NODF_noApis ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.nestedness)

# plotting
# model predictions for plotting
m.nestedness.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "NODF"))
m.nestedness.df$nestedness_pred <- predict(m.nestedness, re.form = NA)
m.nestedness.df$patch <- factor(m.nestedness.df$patch, levels = c("B", "W"))
# plotting
nestedness.pred <- m.nestedness.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = NODF, color = patch), size = 9, alpha = 0.7) + 
  geom_line(aes(x = patch, y = NODF, group = block), linewidth = 2, color = "black", alpha = 0.2) +
  geom_line(aes(x = patch, y = nestedness_pred, group = block), linewidth = 4, linetype = 2) +
  scale_x_discrete(labels = c('Connected', 'Winged')) +
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

m.h2 <- glmmTMB(h2_noApis ~ patch + (1|block),
                data = network_metrics,
                family = "gaussian")
summary(m.h2)

# plotting
# model predictions for plotting
m.h2.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "h2_noApis"))
m.h2.df$h2.pred <- predict(m.h2, re.form = NA)
m.h2.df$patch <- factor(m.h2.df$patch, levels = c("B", "W"))
# plotting
h2.pred <- m.h2.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = h2_noApis, color = patch), size = 9, alpha = 0.7) + 
  geom_line(aes(x = patch, y = h2_noApis, group = block), linewidth = 2, color = "black", alpha = 0.2) +
  geom_line(aes(x = patch, y = h2.pred, group = block), linewidth = 4, linetype = 2) +
  scale_x_discrete(labels = c('Connected', 'Winged')) +
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

m.interaction.div <- glmmTMB(shannon_noApis ~ patch + (1|block),
                             data = network_metrics,
                             family = "gaussian")
summary(m.interaction.div)

#### links per species ####
m.links <- glmmTMB(links.per.species ~ patch + (1|block),
                             data = network_metrics,
                             family = "gaussian")
summary(m.links)
plot(simulateResiduals(m.links))
check_model(m.links)

m.links <- glmmTMB(links_noApis ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.links)


#### linkage density ####
m.density <- glmmTMB(linkage.density ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.density)
plot(simulateResiduals(m.density))
check_model(m.density)

m.density <- glmmTMB(density_noApis ~ patch + (1|block),
                     data = network_metrics,
                     family = "gaussian")
summary(m.density)

# plotting
# model predictions for plotting
m.density.df <- network_metrics %>%
  dplyr::select(c("block", "patch", "density_noApis"))
m.density.df$linkage.density.pred <- predict(m.density, re.form = NA)
m.density.df$patch <- factor(m.density.df$patch, levels = c("B", "W"))
# plotting
linkage.density.pred <- m.density.df %>%
  ggplot() +
  geom_point(aes(x = patch, y = density_noApis, color = patch), size = 9, alpha = 0.7) + 
  geom_line(aes(x = patch, y = density_noApis, group = block), linewidth = 2, color = "black", alpha = 0.2) +
  geom_line(aes(x = patch, y = linkage.density.pred, group = block), linewidth = 4) +
  scale_x_discrete(labels = c('Connected', 'Winged')) +
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
m.floral.div <- glmmTMB(floral_diversity ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.floral.div)
plot(simulateResiduals(m.floral.div))
check_model(m.floral.div)

# no Apis
m.floral.div <- glmmTMB(floral_div_noApis ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.floral.div)
plot(simulateResiduals(m.floral.div))
check_model(m.floral.div)

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

# no Apis
m.pollinator.div <- glmmTMB(pollinator_div_noApis ~ patch + (1|block),
                            data = network_metrics,
                            family = "gaussian")
summary(m.pollinator.div)
plot(simulateResiduals(m.pollinator.div))
check_model(m.pollinator.div)

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










