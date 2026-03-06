# loading libraries
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table, betapart)

source(here::here(file.path("scripts", "00_functions.R")))

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
  #geom_line(aes(x = x, y = predicted, group = group), data = network_dissimilarity.predict, linewidth = 2, linetype = 1) +
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


#### turnover in plant species ####
floral_matrix <- pollinator %>%
  dplyr::count(unique_ID, flower_species) %>%
  dplyr::mutate(n = if_else(n == 0, 0, 1)) %>%
  pivot_wider(names_from = flower_species, values_from = n, values_fill = 0) %>%
  column_to_rownames("unique_ID")

floral_turnover <- beta.pair(floral_matrix, index.family = "sorensen")
floral_turnover_df <- as.data.frame(as.table(as.matrix(floral_turnover[["beta.sor"]])))
floral_turnover_df <- floral_turnover_df %>%
  filter(!(Var1 == Var2)) %>% # removing observations comparing patch to itself
  separate(Var1, into = c("block1", "patch1")) %>%
  separate(Var2, into = c("block2", "patch2")) %>%
  filter(block1 == block2) %>%
  dplyr::rename(floral_turnover = Freq, block = block1) %>%
  dplyr::select(block, floral_turnover)
floral_turnover_df <- floral_turnover_df[!duplicated(floral_turnover_df$floral_turnover), ]

#### turnover in pollinator species ####
pollinator_matrix <- pollinator %>%
  dplyr::count(unique_ID, pollinator_species) %>%
  dplyr::mutate(n = if_else(n == 0, 0, 1)) %>%
  pivot_wider(names_from = pollinator_species, values_from = n, values_fill = 0) %>%
  column_to_rownames("unique_ID")

pollinator_turnover <- beta.pair(pollinator_matrix, index.family = "sorensen")
pollinator_turnover_df <- as.data.frame(as.table(as.matrix(pollinator_turnover[["beta.sor"]])))
pollinator_turnover_df <- pollinator_turnover_df %>%
  filter(!(Var1 == Var2)) %>% # removing observations comparing patch to itself
  separate(Var1, into = c("block1", "patch1")) %>%
  separate(Var2, into = c("block2", "patch2")) %>%
  filter(block1 == block2) %>%
  dplyr::rename(pollinator_turnover = Freq, block = block1) %>%
  dplyr::select(block, pollinator_turnover)
pollinator_turnover_df <- pollinator_turnover_df[!duplicated(pollinator_turnover_df$pollinator_turnover), ]


## all together
species_turnover <- floral_turnover_df %>%
  left_join(pollinator_turnover_df, by = "block")


