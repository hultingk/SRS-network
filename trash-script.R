


###### TRASH SCRIPT #######
m.nested <- glmmTMB(`Shannon diversity` ~ patch + (1|block),
                    data = shannon, 
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







#### community composition? ####
community.comp <- pollinator %>%
  #mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  #mutate(interaction = paste(pollinator_analysis, flower_species, sep = "-")) %>%
  #count(unique_ID, interaction) %>%
  count(unique_ID, pollinator_analysis) %>%
  pivot_wider(names_from = pollinator_analysis, values_from = n, values_fill = 0)
community.comp.no.site <- community.comp[-1] # deleting site column

set.seed(300) # setting seed to reproduce later
# Choosing number of dimensions
# Specifying the raw data models for all dimensions k
mods.r <- list() # an empty list to hold model objects
stress.r <- numeric() # an empty numeric vector to hold final stress values 
conv.r <- character() # an empty character vector to hold convergence messages 
for(i in 1:7){
  mod <- metaMDS(community.comp.no.site, distance = 'bray', k = i, try = 20, trymax = 1000)
  mods.r[[as.character(i)]] <- mod 
  stress.r[i] <- mod$stress 
  conv.r[i] <- mod$converged
}
# scree plot
palette(c('red', 'black'))
plot(x = 1:7, y = stress.r, main = 'Raw Data', xlab = 'Dimensionality',
     ylab = 'Stress',
     pch = 19, col = factor(conv.r))
# red points are dimensions that there was no stable solution for


# Actual NMDS
comp.mds <- metaMDS(community.comp.no.site, distance = "bray", k = 3, autotransform = F, pc = T, trymax = 1000)

# Calculating how much variance in community data the NMDS accounts for
bcd <- vegdist(community.comp.no.site, method = 'bray')
r.ed <- dist(comp.mds$points[, 1:2]) # Calculating distance between sites according to the first two axes
cor(r.ed, bcd, method = 'spearman') # Calculating Spearman's rho

# Shepard plot, shows fit
stressplot(comp.mds) # looks good


# PERMANOVA patch type - testing if community composition differs between connected and winged patches
obs <- community.comp.no.site # community data
spp <- community.comp %>% 
  dplyr::select(unique_ID) %>%
  separate(unique_ID, c("block", "patch")) %>%
  dplyr::select(patch) # separating patch type

d.manova <- adonis2(obs ~ patch, method = "bray", data= spp) # testing difference between patch types
d.manova # summary

#### DIFFERENCES IN DISPERSION
# Compute distance to group centroid (dispersion analysis)
dispersion <- betadisper(vegdist(obs, method = "bray"), spp$patch)

# Perform a permutation test for homogeneity of dispersions
dispersion_test <- permutest(dispersion, permutations = 999)

# Print test results
print(dispersion_test)

# Plot the dispersion results
plot(dispersion, main = "Group Dispersion in NMDS Space")

# Enhanced ggplot visualization of NMDS with group dispersions
nmds_df <- as.data.frame(comp.mds$points)
nmds_df$Group <- spp$patch

ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal() +
  labs(title = "NMDS Plot with Group Dispersions",
       x = "NMDS1", y = "NMDS2")







# plotting NMDS
plot.community <- community.comp %>%
  separate(unique_ID, c("block", "patch")) %>%
  mutate(patch.color = if_else(patch == "B", 1, 2)) # adding #s that will correspond to plot colors

palette(c('thistle', 'lightblue')) # setting plot colors
#plotting
par(mar = c(5.1, 4.5, 4.1, 2.1))
p1 <- plot(x = comp.mds$points[, 1], # x axis = 1st dimension
           y = comp.mds$points[, 2], # y axis = 2nd dimension
           pch = 19, # point shapes
           cex = 2, # point size
           col = plot.community$patch.color, # point colors, correspond to patch type
           xlab = 'NMDS Axis 1', # x-axis label
           ylab = 'NMDS Axis 2', # y-axis label
           cex.axis = 1.5, # axis number size
           cex.lab = 1.9#, # axis label size
           #xlim = c(-2.4, 1.7), # plot x-axis size
           #ylim = c(-1.5, 1.4)
) # plot y-axis size
# Plotting ellipse polygons 
ordiellipse(comp.mds, 
            groups = plot.community$patch, # groups are patch type
            col = plot.community$patch.color, # ellipse colors, correspond to patch type
            draw = 'polygon', # type of ellipse
            kind = 'sd', # ellipse is drawing standard deviation of points
            alpha = 140, # transparency of ellipse
            border = plot.community$patch.color)# border color, correspond to patch type
legend("topleft", # adding legend to top left corner
       inset=c(0, -0.03), # exact legend locations
       legend = c('Connected Patch', 'Winged Patch'), # legend names
       bty = 'n', # no box drawn around legend
       pch = 19, # legend point shape
       col = plot.community$patch.color, # legend colors, correspond to patch type
       cex = 1.8, # legend size
       y.intersp = 0.4,
       x.intersp = 0.3) # vertical distance between legend lines
text(-2.1, -1.5, cex = 1.5, paste("stress =", round(comp.mds$stress, digits = 3))) # adding stress value to plot
nmds.plot <- recordPlot() # saving plot as an object











