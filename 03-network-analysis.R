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

pollinator <- pollinator %>%
  filter(!pollinator_analysis %in% c("Apis mellifera")) %>%
  filter(!pollinator_species == "") %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "pollinator_analysis", "flower_species")) 



#### network analysis ####
pollinator_split <- pollinator %>%
  count(unique_ID, pollinator_analysis, flower_species) %>%
  group_by(unique_ID) %>%
  group_split() 
  

prepare_matrix <- function(df) {
  df_wide <- df %>% 
    pivot_wider(names_from = pollinator_analysis, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    column_to_rownames("flower_species") #convert years to rownames
}

webs <- pollinator_split %>%
  lapply(prepare_matrix)
  
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")
names(webs) <- webs.names


# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'nestedness') 
# Calculate network metric links per species for all plant-pollinator sites
net.metrics.h2 <- lapply(webs, networklevel, index = 'H2') 


# Make null models for all sites using the r2dtable null
#net.nulls.r2d <- lapply(webs, nullmodel, method = "r2dtable", N = 500) 
# Make null models for all sites using the vaznull null
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = 500) 
# Make null models for all sites using the swap.web null
#net.nulls.swap <- lapply(webs, nullmodel, method = "swap.web", N = 500)


# Null distribution function for nestedness - calculates the network nestedness for each null (using a particular null method) for each site 
net.null.nest = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'nestedness'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}

# Null distribution function for links per species - calculates the network links per species metric for each null (using a particular null method) for each site 
net.null.h2 = function(nulls){
  net.null.metric <- list()
  for (i in 1:length(nulls)) {
    net.null.metric[[i]] = do.call('rbind', 
                                   lapply(nulls[[i]], networklevel, index = 'H2'))
  }
  names(net.null.metric) <- webs.names
  return(net.null.metric)
}

vaz.nest <- net.null.nest(net.nulls.vaz)
vaz.h2 <- net.null.h2(net.nulls.vaz)

# Z-score function for comparing different networks
net.zscore = function(obsval, nullval) {
  (obsval - mean(nullval))/sd(nullval)  
} 

# Function that perform z-score calculation of nestedness using the observed and null networks
nest.zscore = function(nulltype){
  net.nest.zscore <- list() 
  for(i in 1:length(net.metrics.nest)){
    net.nest.zscore[[i]] = net.zscore(net.metrics.nest[[i]]['nestedness'], 
                                      nulltype[[i]][ ,'nestedness'])
  }
  names(net.nest.zscore) <- webs.names
  return(net.nest.zscore)
}

# Function that perform z-score calculation of links per species using the observed and null networks
h2.zscore = function(nulltype){
  net.h2.zscore <- list() 
  for(i in 1:length(net.metrics.h2)){
    net.h2.zscore[[i]] = net.zscore(net.metrics.h2[[i]]['H2'], 
                                       nulltype[[i]][ ,'H2'])
  }
  names(net.h2.zscore) <- webs.names
  return(net.h2.zscore)
}


vaz.nest.zscore <- nest.zscore(vaz.nest)
vaz.h2.zscore <- h2.zscore(vaz.h2)

nestedness <- as.data.frame(do.call('rbind', vaz.nest.zscore) )
nestedness <- nestedness %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into columns
  
h2 <- as.data.frame(do.call('rbind', vaz.h2.zscore) )
h2 <- h2 %>%
  rownames_to_column(var = "unique_ID") %>%
  separate(unique_ID, into = c("block", "patch"))  # seperating unique ID into column

network_metrics <- nestedness %>%
  left_join(h2, by = c("block", "patch"))

m.nested <- glmmTMB(nestedness ~ patch,
                    data = network_metrics, 
                    family = "gaussian")
summary(m.nested)
plot(simulateResiduals(m.nested))
check_model(m.nested)

m.h2 <- glmmTMB(H2 ~ patch + (1|block),
                    data = network_metrics, 
                    family = "gaussian")
summary(m.h2)
plot(simulateResiduals(m.h2))
check_model(m.h2)




#pollinator_wider <- pollinator %>%
#  filter(!pollinator_species == "") %>%
#  filter(patch == "B") %>%
#  count(pollinator_species, flower_species) %>%
#  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
 # column_to_rownames(var="flower_species")
#plotweb(pollinator_wider)
#### Interaction abundance ####
# calculate abundance of interactions per patch and sampling round
abundance <- pollinator %>%
  count(block, patch)

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
    geom_line(aes(x = patch, y = n, group = block), linewidth = 1.5, color = "black", alpha = 0.2) +
    geom_line(aes(x = patch, y = abund_pred, group = block), linewidth = 2) +
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
m2.predict <- ggpredict(m2, terms=c("patch [all]"), back_transform = T)
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
 # ylim(0.72, 0.95) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 26)) + # axis tick mark size
    theme(axis.title = element_text(size = 30)) #+ # axis label size
fl.div.pred

pdf(file = file.path("plots", "floral-diversity.pdf"), width = 10, height = 8)
fl.div.pred
dev.off()




##### butterfly diversity ####
butterfly_wider <- pollinator %>%
  filter(order == "Lepidoptera") %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  # filter(unique_ID != "57.W") %>%
  select(!unique_ID) # removing unique ID in order to rarefy

butterfly_data <- pollinator %>% # adding rarefied floral data to block/patches
  count(block, patch) %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) #%>% # creating a unique ID for each patch
butterfly.div <- diversity(butterfly_wider)
butterfly_data$butterfly.div <- butterfly.div

m3 <- glmmTMB(butterfly.div ~ patch + (1|block), # block is random effect
              data = butterfly_data, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(m3) # model summary
plot(simulateResiduals(m3)) # looks good




#### butterfly interaction diversity ####
interaction.div <- pollinator %>% # converting floral interaction data into a community matrix
  filter(order == "Lepidoptera") %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  mutate(interaction = paste(pollinator_species, flower_species, sep = "-")) %>%
  count(unique_ID, interaction) %>%
  pivot_wider(names_from = interaction, values_from = n, values_fill = 0) %>%
  #filter(unique_ID != "57.W") %>%
  select(!unique_ID) # removing unique ID in order to rarefy

interaction.div <- diversity(interaction.div)
butterfly_data$interaction.div <- interaction.div
m4 <- glmmTMB(interaction.div ~ patch + (1|block), # block is random effect
              data = butterfly_data, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(m4) # model summary
plot(simulateResiduals(m4)) # looks good







#filter(unique_ID != "57.W")
rare.data$rare <- sRare # adding rarefied floral data to block/patches
butterfly.div <- diversity(butterfly.div)
rare.data$butterfly.div <- butterfly.div
# diversity of butterfly interactions
interaction.div <- butterfly %>% # converting floral interaction data into a community matrix
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  mutate(interaction = paste(pollinator_species, flower_species, sep = "-")) %>%
  count(unique_ID, interaction) %>%
  pivot_wider(names_from = interaction, values_from = n, values_fill = 0) %>%
  #filter(unique_ID != "57.W") %>%
  select(!unique_ID) # removing unique ID in order to rarefy
spAbund <- rowSums(interaction.div) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare <- rarefy(interaction.div, 9) # now use function rarefy
rare.data <- butterfly %>% # adding rarefied floral data to block/patches
  count(block, patch) %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) #%>% # creating a unique ID for each patch
#filter(unique_ID != "57.W")
rare.data$rare <- sRare # adding rarefied floral data to block/patches
interaction.div <- diversity(interaction.div)
rare.data$interaction.div <- interaction.div
m2 <- glmmTMB(interaction.div ~ patch + (1|block), # block is random effect
              data = rare.data, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(m2) # model summary
plot(simulateResiduals(m2)) # looks good


# diversity of butterflies
butterfly.div <- butterfly %>% # converting floral interaction data into a community matrix
  mutate(unique_ID = paste(block, patch, sep = ".")) %>% # creating a unique ID for each patch
  count(unique_ID, pollinator_species) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  # filter(unique_ID != "57.W") %>%
  select(!unique_ID) # removing unique ID in order to rarefy
spAbund <- rowSums(butterfly.div) # calculating minimum # of observation 
min(spAbund) # 7 is the fewest interactions observed per patch
sRare <- rarefy(butterfly.div, 9) # now use function rarefy
rare.data <- butterfly %>% # adding rarefied floral data to block/patches
  count(block, patch) %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) #%>% # creating a unique ID for each patch
#filter(unique_ID != "57.W")
rare.data$rare <- sRare # adding rarefied floral data to block/patches
butterfly.div <- diversity(butterfly.div)
rare.data$butterfly.div <- butterfly.div
m3 <- glmmTMB(butterfly.div ~ patch + (1|block), # block is random effect
              data = rare.data, 
              family = "gaussian") # gaussian family for non-integer continuous data
summary(m3) # model summary
plot(simulateResiduals(m3)) # looks good





