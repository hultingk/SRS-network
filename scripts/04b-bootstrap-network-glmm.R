librarian::shelf(glmmTMB, DHARMa, tidyverse, ggeffects)

network_metrics_boot <- read.csv(file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"))


#### H2 ####
m_h2 <- glmmTMB(h2 ~ patch + (1|block),
                data = network_metrics_boot)
summary(m_h2)
plot(simulateResiduals(m_h2))

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
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Network specialization (H2')")) +
  theme(legend.position = "none") 
h2.pred

# pdf(file = file.path("plots", "network-plots", "h2.pred.pdf"), width = 7, height = 6)
# h2.pred
# dev.off()




#### nestedness ####
m_nodf <- glmmTMB(nodf ~ patch + (1|block),
                  data = network_metrics_boot)
summary(m_nodf)
plot(simulateResiduals(m_nodf))

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
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Nestedness (NODF)")) +
  theme(legend.position = "none") 
nodf.pred

# pdf(file = file.path("plots", "network-plots", "nodf.pred.pdf"), width = 7, height = 6)
# nodf.pred
# dev.off()




#### links per species ####
m_links <- glmmTMB(links_per_sp ~ patch + (1|block),
                   data = network_metrics_boot)
summary(m_links)
plot(simulateResiduals(m_links))

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
  theme_classic(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Links per species")) +
  theme(legend.position = "none") 
links.pred

# pdf(file = file.path("plots", "network-plots", "links.pred.pdf"), width = 7, height = 6)
# links.pred
# dev.off()



#### pollinator niche overlap ####
m_HL_niche <- glmmTMB(HL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_HL_niche)
plot(simulateResiduals(m_HL_niche))

# plotting links per species
m.HL_niche.df <- ggpredict(m_HL_niche, terms = c("patch"), back_transform = TRUE)
# plotting
HL_niche.pred <- m.HL_niche.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = HL_niche, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.HL_niche.df, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
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
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Pollinator niche overlap")) +
  theme(legend.position = "none") 
HL_niche.pred

# pdf(file = file.path("plots", "network-plots", "HL_niche.pred.pdf"), width = 7, height = 6)
# HL_niche.pred
# dev.off()



#### plant niche overlap ####
m_LL_niche <- glmmTMB(LL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_LL_niche)
plot(simulateResiduals(m_LL_niche))

# plotting links per species
m.LL_niche.df <- ggpredict(m_LL_niche, terms = c("patch"), back_transform = TRUE)
# plotting
LL_niche.pred <- m.LL_niche.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = LL_niche, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.LL_niche.df, width = 0, linewidth = 2.5) +
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
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Plant niche overlap")) +
  theme(legend.position = "none") 
LL_niche.pred

# pdf(file = file.path("plots", "LL_niche.pred.pdf"), width = 7, height = 6.5)
# LL_niche.pred
# dev.off()



# ##### HL robustness ####
# m_HL_robustness <- glmmTMB(HL_robustness ~ patch + (1|block),
#                      data = network_metrics_boot)
# summary(m_HL_robustness)
# plot(simulateResiduals(m_HL_robustness))
# 
# # plotting links per species
# m.HL_robustness.df <- ggpredict(m_HL_robustness, terms = c("patch"), back_transform = TRUE)
# # plotting
# HL_robustness.pred <- m.HL_robustness.df %>%
#   ggplot() +
#   geom_jitter(aes(x = patch, y = HL_robustness, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
#               width = 0.08, height = 0) +
#   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
#                 data = m.HL_robustness.df, width = 0, linewidth = 2.5) +
#   geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 2) +
#   geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   theme_classic(base_size = 26) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
#         # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_line(color = "black", linewidth = 0.7),
#         strip.text.x = element_text(hjust = -0.05),
#         panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
#         plot.background = element_rect(fill = "transparent", color = NA)) +
#   scale_color_manual(values=c("#5B859EFF","#DCB254")) +
#   scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
#   xlab("Patch type") +
#   ylab(expression("Pollinator robustness to plant extinctions")) +
#   theme(legend.position = "none") 
# HL_robustness.pred
# 
# pdf(file = file.path("plots", "network-plots", "HL_robustness.pred.pdf"), width = 7, height = 6)
# HL_robustness.pred
# dev.off()
# 
# 
# 
# 
# #### LL robustness #####
# m_LL_robustness <- glmmTMB(LL_robustness ~ patch + (1|block),
#                            data = network_metrics_boot)
# summary(m_LL_robustness)
# plot(simulateResiduals(m_LL_robustness))
# 
# # plotting robustness to plant extinctions
# # plotting LL_robustness
# m.LL_robustness.df <- ggpredict(m_LL_robustness, terms = c("patch"), back_transform = TRUE)
# # plotting
# LL_robustness.pred <- m.LL_robustness.df %>%
#   ggplot() +
#   geom_jitter(aes(x = patch, y = LL_robustness, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
#               width = 0.08, height = 0) +
#   geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
#                 data = m.LL_robustness.df, width = 0, linewidth = 2.5) +
#   geom_line(aes(x = x, y = predicted, group = group), linewidth = 2, linetype = 1) +
#   geom_point(aes(x = x, y = predicted, fill = x), size = 8, colour="black", pch=21, stroke = 2) +
#   scale_x_discrete(labels = c('Connected', 'Unconnected')) +
#   theme_classic(base_size = 26) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
#         # panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_line(color = "black", linewidth = 0.7),
#         strip.text.x = element_text(hjust = -0.05),
#         panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
#         plot.background = element_rect(fill = "transparent", color = NA)) +
#   scale_color_manual(values=c("#5B859EFF","#DCB254")) +
#   scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
#   xlab("Patch type") +
#   ylab(expression("Robustness to plant extinction")) +
#   theme(legend.position = "none") 
# LL_robustness.pred
# 
# # pdf(file = file.path("plots", "LL_robustness.pred.pdf"), width = 7, height = 6.5)
# # LL_robustness.pred
# # dev.off()



#### shannon diversity ####
m_shannon <- glmmTMB(shannon ~ patch + (1|block),
                     data = network_metrics_boot)
summary(m_shannon)
plot(simulateResiduals(m_shannon))


# plotting shannon diversity of interactions
m.shannon.df <- ggpredict(m_shannon, terms = c("patch"), back_transform = TRUE)
# plotting
shannon.pred <- m.shannon.df %>%
  ggplot() +
  geom_jitter(aes(x = patch, y = shannon, color = patch), data = network_metrics_boot, size = 6, alpha = 0.55,
              width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                data = m.shannon.df, width = 0, linewidth = 2.5) +
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
  scale_color_manual(values=c("#5B859EFF","#DCB254")) +
  scale_fill_manual(values=c("#5B859EFF","#DCB254")) +
  xlab("Patch type") +
  ylab(expression("Shannon diversity of interactions")) +
  theme(legend.position = "none") 
shannon.pred

# pdf(file = file.path("plots", "network-plots", "shannon.pred.pdf"), width = 7, height = 6)
# shannon.pred
# dev.off()


# all together
network_plot <- cowplot::plot_grid(shannon.pred, h2_pred, nodf.pred,
                                     labels = c("A)", "B)", "C)", "D)"), label_size = 26, label_x = 0.19, label_y = 0.95)
network_plot

pdf(file = file.path("plots", "network-plots", "network_plot.pdf"), width = 13, height = 11)
network_plot
dev.off()



