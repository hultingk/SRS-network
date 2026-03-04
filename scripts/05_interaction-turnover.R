# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table)

# loading data 
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))

# creating unique ID
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 


#### network analysis: WITH ALL SPECIES ###
pollinator_split <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species, flower_species) %>%
  dplyr::group_by(unique_ID) %>%
  group_split() 

# getting each network into correct format
webs <- pollinator_split %>%
  lapply(prepare_matrix)

# adding patch names to each network
webs.names <- c("10.B","10.W","52.B","52.W", "53N.B", "53N.W",
                "53S.B", "53S.W", "54S.B", "54S.W", "57.B", "57.W", "8.B", "8.W")
names(webs) <- webs.names

#### interaction rewiring and turnover ####
network_dissimilarity <- betalinkr_multi(webs2array(list("10.B" = as.matrix(webs[[1]]), "10.W" = as.matrix(webs[[2]]),
                                                         "52.B" = as.matrix(webs[[3]]), "52.W" = as.matrix(webs[[4]]),
                                                         "53N.B" = as.matrix(webs[[5]]), "53N.W" = as.matrix(webs[[6]]),
                                                         "53S.B" = as.matrix(webs[[7]]), "53S.W" = as.matrix(webs[[8]]),
                                                         "54S.B" = as.matrix(webs[[9]]), "54S.W" = as.matrix(webs[[10]]),
                                                         "57.B" = as.matrix(webs[[11]]), "57.W" = as.matrix(webs[[12]]),
                                                         "8.B" = as.matrix(webs[[13]]), "8.W" = as.matrix(webs[[14]]))), 
                                         index = "sorensen", partitioning="commondenom")

# reformatting
network_dissimilarity <- network_dissimilarity %>%
  mutate(block_pair = paste(i, j, sep = "-")) %>%
  mutate(real_pairs = case_when(
    block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
                      "54S.B-54S.W", "57.B-57.W", "8.B-8.W") ~ "real", 
    .default = "random"
  )) %>%
  #filter(block_pair %in% c("10.B-10.W", "52.B-52.W", "53N.B-53N.W", "53S.B-53S.W",
  #                        "54S.B-54S.W", "57.B-57.W", "8.B-8.W")) %>%
  pivot_longer(cols = S:ST, names_to = "type", values_to = "dissimilarity")

# factor levels
network_dissimilarity$type <- factor(network_dissimilarity$type, levels=c("S", "WN", "ST", "OS"))

# filtering for just interaction rewiring and species turnover, filtering for pairs within a block
network_dissimilarity <- network_dissimilarity %>%
  filter(!type %in% c("S", "WN")) %>%
  filter(real_pairs == "real") 
# basic model for plotting
m.network_dissimilarity <- glmmTMB(dissimilarity ~ type,
                                   data = network_dissimilarity)
summary(m.network_dissimilarity)

# model predictions for plotting
network_dissimilarity.predict <- ggpredict(m.network_dissimilarity, terms = c("type"), back_transform = TRUE)

# plotting
dissimilarity_plot <- network_dissimilarity %>%
  ggplot() +
  geom_jitter(aes(x = type, y = dissimilarity, color = type), data = network_dissimilarity, size = 5, alpha = 0.55,
                width = 0.08, height = 0) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = x), color = "black",
                   data = network_dissimilarity.predict, width = 0, linewidth = 2.5) +
  geom_line(aes(x = x, y = predicted, group = group), data = network_dissimilarity.predict, linewidth = 2, linetype = 1) +
  geom_point(aes(x = x, y = predicted, fill = x), data = network_dissimilarity.predict, 
               size = 6, colour="black", pch=21, stroke = 2) +
  scale_x_discrete(labels = c('Species Turnover', 'Interaction Rewiring')) +
  scale_color_manual(values=c("#F5097C","#F7B3D4")) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4")) +
  #scale_x_discrete(labels = c(expression(beta[WN]), expression(beta[ST]), expression(beta[OS]))) +
  # scale_x_discrete(labels = c(expression("Species turnover"), expression("Interaction rewiring"))) +
  theme_classic(base_size = 20) +
  xlab("Dissimilarity component") +
  theme(legend.position = "none") +
  ylab(expression(paste("Dissimilarity between patch pairs"))) 
dissimilarity_plot

# exporting
# pdf(file = file.path("plots", "dissimilarity_plot.pdf"), width = 6, height = 5.5)
# dissimilarity_plot
# dev.off()






