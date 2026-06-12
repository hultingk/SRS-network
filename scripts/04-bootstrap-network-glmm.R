librarian::shelf(glmmTMB, DHARMa, tidyverse, performance, ggeffects)

network_metrics_boot <- read.csv(file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"))

# nestedness
m_nodf <- glmmTMB(nodf ~ patch + (1|block),
                  data = network_metrics_boot)
summary(m_nodf)
# nestedness plot
# plotting nodf per species
m.nodf.df <- ggpredict(m_nodf, terms = c("patch"), back_transform = TRUE)
# plotting
nodf.pred <- m.nodf.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = nodf, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.nodf.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 28) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#1E395FFF","#5B859EFF")) +
  scale_fill_manual(values=c("#1E395FFF","#5B859EFF")) +
  xlab("Patch type") +
  ylab(expression("NODF")) +
  theme(legend.position = "none") 
nodf.pred

# pdf(file = file.path("plots", "nodf.pred.pdf"), width = 7, height = 6.5)
# nodf.pred
# dev.off()


# links per species 
m_links <- glmmTMB(links_per_sp ~ patch + (1|block),
                   data = network_metrics_boot)
summary(m_links)

# plotting links per species
m.links.df <- ggpredict(m_links, terms = c("patch"), back_transform = TRUE)
# plotting
links.pred <- m.links.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = links_per_sp, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.links.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 28) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#1E395FFF","#5B859EFF")) +
  scale_fill_manual(values=c("#1E395FFF","#5B859EFF")) +
  xlab("Patch type") +
  ylab(expression("Links per species")) +
  theme(legend.position = "none") 
links.pred

# pdf(file = file.path("plots", "links.pred.pdf"), width = 7, height = 6.5)
# links.pred
# dev.off()



# pollinator niche overlap
m_HL_niche <- glmmTMB(HL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_HL_niche)

# plant niche overlap
m_LL_niche <- glmmTMB(LL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_LL_niche)

# web asymmetry
m_asymmetry <- glmmTMB(asymmetry ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_asymmetry)

# connectance
m_connectance <- glmmTMB(connectance ~ patch + (1|block),
                       data = network_metrics_boot)
summary(m_connectance)

# linkage density
m_density <- glmmTMB(linkage_density ~ patch + (1|block),
                         data = network_metrics_boot)
summary(m_density)

# plotting linkage density
m.linkage.df <- ggpredict(m_density, terms = c("patch"), back_transform = TRUE)
# plotting
linkage.pred <- m.linkage.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = linkage_density, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.linkage.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 28) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#1E395FFF","#5B859EFF")) +
  scale_fill_manual(values=c("#1E395FFF","#5B859EFF")) +
  xlab("Patch type") +
  ylab(expression("Linkage density")) +
  theme(legend.position = "none") 
linkage.pred

# pdf(file = file.path("plots", "linkage.pred.pdf"), width = 7, height = 6.5)
# linkage.pred
# dev.off()



# HL robustness 
m_HL_robustness <- glmmTMB(HL_robustness ~ patch + (1|block),
                     data = network_metrics_boot)
summary(m_HL_robustness)

# LL robustness 
m_LL_robustness <- glmmTMB(LL_robustness ~ patch + (1|block),
                           data = network_metrics_boot)
summary(m_LL_robustness)

# plotting robustness to plant extinctions
# plotting LL_robustness
m.LL_robustness.df <- ggpredict(m_LL_robustness, terms = c("patch"), back_transform = TRUE)
# plotting
LL_robustness.pred <- m.LL_robustness.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = LL_robustness, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.LL_robustness.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Connected', 'Unconnected')) +
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#1E395FFF","#5B859EFF")) +
  scale_fill_manual(values=c("#1E395FFF","#5B859EFF")) +
  xlab("Patch type") +
  ylab(expression("Robustness to plant extinction")) +
  theme(legend.position = "none") 
LL_robustness.pred

# pdf(file = file.path("plots", "LL_robustness.pred.pdf"), width = 7, height = 6.5)
# LL_robustness.pred
# dev.off()



