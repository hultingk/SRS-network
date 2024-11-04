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

#### Interaction abundance ####
# calculate abundance of interactions per patch and sampling round
abundance <- pollinator %>%
  count(block, patch, sampling_round)

# how does connectivity affect the abundance of plant-pollinator interactions? 
abund_mod <- glmmTMB(n ~ patch + (1|block) + (1|sampling_round), # block is random effect
              data = abundance,
              family = nbinom2) # negative binomial family for overdispersed count data
summary(abund_mod) # model output
plot(simulateResiduals(abund_mod)) # looks good
check_overdispersion(abund_mod) # no overdispersion for nbinom, overdispersed when family = poisson


#### Floral diversity ####
floral_wider <- pollinator %>% # converting floral interaction data into a community matrix
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  count(unique_ID, flower_species) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  select(!unique_ID) # removing unique ID in order to rarefy

fl.diversity <- diversity(floral_wider) # calculating diversity of flowers that are interacted with
summarize_pollinator <- pollinator %>%
  count(block, patch)
summarize_pollinator$fl.diversity <- fl.diversity

# how does connectivity affect the diversity of flowers that pollinators forage on?
fl_div_mod <- glmmTMB(fl.diversity ~ patch + (1|block), # block is random effect
              data = summarize_pollinator, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(fl_div_mod) # model summary
plot(simulateResiduals(fl_div_mod)) # looks good


# plotting richness model predictions
fl_div_mod.predict <- ggpredict(fl_div_mod, terms=c("patch [all]"), back.transform = T, allow.new.levels=TRUE)
fl_div_mod.stat.test <- tibble::tribble( # creating tibble of p-values for plotting
  ~group1, ~group2,   ~p.adj,
  "B",     "W", "0.01"
)
fl_div_mod_plot <- fl_div_mod.predict %>%
  ggplot(aes(x = x, y = predicted)) + # plotting predicted values
  geom_jitter(aes(x = patch, y = fl.diversity), data = summarize_pollinator, alpha = 0.2, width = 0.1, height = 0.1, size = 5)+ # adding raw data points
  geom_point(size = 5)+  # size of center prediction point
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, linewidth = 2.5) + # adding error bars of predicted confidence intervals
  theme_classic() + # choosing simplest ggplot theme
  stat_pvalue_manual( # adding p-values to plot
    fl_div_mod.stat.test, 
    size = 9,# size of letters
    bracket.size = 1.5,
    y.position = 3.2, step.increase = 0.12, # position of brakets
    label = "p.adj"
  ) +
  labs(title = NULL, # no title
       x = "Patch Type", # changing axis labels
       y = "Floral Foraging Richness") +
  scale_x_discrete(labels= c("Connected", "Winged")) + # replacing "B" and "W" with connected and winged
  theme(axis.text = element_text(size = 26)) + # axis tick mark size
  theme(axis.title = element_text(size = 30)) #+ # axis label size
# ylim(2.5, 12) # manually setting y axis limits to include p-value brackets
fl_div_mod_plot






