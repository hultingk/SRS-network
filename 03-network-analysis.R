# -------------------------------------- #
#### Network analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libaries
library(tidyverse)
library(bipartite)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggeffects)
library(ggpubr)


# loading data 
pollinator <- read.csv(file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))

pollinator %>%
  count(pollinator_species)

pollinator_wider <- pollinator %>%
  filter(!pollinator_species == "") %>%
  filter(patch == "B") %>%
  count(pollinator_species, flower_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  column_to_rownames(var="flower_species")
plotweb(pollinator_wider)
#### Interaction abundance ####
# calculate abundance of interactions per patch and sampling round
abundance <- pollinator %>%
  count(block, patch, flower_species)

# how does connectivity affect the abundance of plant-pollinator interactions? 
m1 <- glmmTMB(n ~ patch + (1|block), # block is random effect
              data = abundance,
              family = nbinom2) # negative binomial family for overdispersed count data
summary(m1) # model output
plot(simulateResiduals(m1)) # looks good
check_overdispersion(m1) # no overdispersion for nbinom, overdispersed when family = poisson


# model predictions for plotting
m1.df <- abundance %>%
  dplyr::select(c("block", "patch", "n"))
m1.df$abund_pred <- exp(predict(m1, re.form = NA))
m1.df$patch <- factor(m1.df$patch, levels = c("B", "W"))
# plotting
interaction.abun.pred <- m1.df %>%
    ggplot() +
    geom_point(aes(x = patch, y = n, color = patch), size = 7, alpha = 0.7) + 
    geom_line(aes(x = patch, y = n, group = block), size = 1.5, color = "black", alpha = 0.2) +
    geom_line(aes(x = patch, y = abund_pred, group = block), size = 2) +
    scale_x_discrete(labels = c('Connected', 'Winged')) +
    scale_color_manual(values=c("#506D8F","#E2A03C")) +
    xlab("Patch Type") +
    ylab(expression(paste("Interaction Abundance"))) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 26)) + # axis tick mark size
    theme(axis.title = element_text(size = 30)) #+ # axis label size
interaction.abun.pred
pdf(file = file.path("plots", "interaction-abund.pdf"), width = 10, height = 8)
interaction.abun.pred
dev.off()


#### Floral diversity ####
floral_wider <- pollinator %>% # converting floral interaction data into a community matrix
  filter(!pollinator_species %in% c("Apis mellifera")) %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID in order to rarefy

fl.diversity <- diversity(floral_wider, "simpson") # calculating diversity of flowers that are interacted with
summarize_pollinator <- pollinator %>%
  count(block, patch)
summarize_pollinator$fl.diversity <- fl.diversity

# how does connectivity affect the diversity of flowers that pollinators forage on?
m2 <- glmmTMB(fl.diversity ~ patch + (1|block), # block is random effect
              data = summarize_pollinator, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(m2) # model summary
plot(simulateResiduals(m2)) # looks good


# plotting richness model predictions
m2.predict <- ggpredict(m2, terms=c("patch [all]"), back.transform = T, allow.new.levels=TRUE)
m2.stat.test <- tibble::tribble( # creating tibble of p-values for plotting
  ~group1, ~group2,   ~p.adj,
  "B",     "W", "0.01"
)
m2_plot <- m2.predict %>%
  ggplot(aes(x = x, y = predicted)) + # plotting predicted values
  geom_jitter(aes(x = patch, y = fl.diversity), data = summarize_pollinator, alpha = 0.2, width = 0.1, height = 0.1, size = 5)+ # adding raw data points
  geom_point(size = 5)+  # size of center prediction point
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, linewidth = 2.5) + # adding error bars of predicted confidence intervals
  theme_classic() + # choosing simplest ggplot theme
  stat_pvalue_manual( # adding p-values to plot
    m2.stat.test, 
    size = 9,# size of letters
    bracket.size = 1.5,
    y.position = 1.3, step.increase = 0.12, # position of brakets
    label = "p.adj"
  ) +
  labs(title = NULL, # no title
       x = "Patch Type", # changing axis labels
       y = "Floral Foraging Richness") +
  scale_x_discrete(labels= c("Connected", "Winged")) + # replacing "B" and "W" with connected and winged
  theme(axis.text = element_text(size = 26)) + # axis tick mark size
  theme(axis.title = element_text(size = 30)) #+ # axis label size
# ylim(2.5, 12) # manually setting y axis limits to include p-value brackets
m2_plot

# alternative plots
# model predictions for plotting
m2.df <- summarize_pollinator %>%
  dplyr::select(c("block", "patch", "fl.diversity"))
m2.df$div_pred <- predict(m2, re.form = NA)
m2.df$patch <- factor(m2.df$patch, levels = c("B", "W"))
# plotting
fl.div.pred <- m2.df %>%
    ggplot() +
    geom_point(aes(x = patch, y = fl.diversity, color = patch), size = 7, alpha = 0.7) + 
    geom_line(aes(x = patch, y = fl.diversity, group = block), size = 1.5, color = "black", alpha = 0.2) +
    geom_line(aes(x = patch, y = div_pred, group = block), size = 2) +
    scale_x_discrete(labels = c('Connected', 'Winged')) +
    scale_color_manual(values=c("#506D8F","#E2A03C")) +
    xlab("Patch Type") +
    ylab(expression(paste("Floral Foraging Diversity"))) +
  ylim(0.72, 0.95) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 26)) + # axis tick mark size
    theme(axis.title = element_text(size = 30)) #+ # axis label size
fl.div.pred

pdf(file = file.path("plots", "floral-diversity.pdf"), width = 10, height = 8)
fl.div.pred
dev.off()


