librarian::shelf(glmmTMB, DHARMa, tidyverse, performance, ggeffects)

network_metrics_boot <- read.csv(file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"))


# H2
m_h2 <- glmmTMB(h2 ~ patch + (1|block),
                data = network_metrics_boot)
summary(m_h2)

# plotting h2 per species
m.h2.df <- ggpredict(m_h2, terms = c("patch"), back_transform = TRUE)
# plotting
h2.pred <- m.h2.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = h2, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.h2.df, width = 0, linewidth = 2.5) +
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
  ylab(expression("Network specialization (H2')")) +
  theme(legend.position = "none") 
h2.pred

# pdf(file = file.path("plots", "h2.pred.pdf"), width = 7, height = 6.5)
# h2.pred
# dev.off()




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



