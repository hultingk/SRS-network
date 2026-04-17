# install.packages("remotes")
# remotes::install_github("Pakillo/BayesianWebs")
librarian::shelf(tidyverse, plyr, vegan, bipartite, data.table, BayesianWebs)

# loading functions
source(here::here(file.path("scripts", "00_functions.R")))

# loading data
pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))
# getting into format
pollinator <- pollinator %>%
  mutate(unique_ID = paste(block, patch, sep = ".")) %>%
  dplyr::select(c("unique_ID", "order", "family", "pollinator_species", "flower_species")) 


#### network analysis: WITH ALL SPECIES ####
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


#### example WITHOUT APIS ####
web.53S.B <- as.matrix(webs$`53S.B`)

# plotting
plot_counts_obs(web.53S.B, sort = FALSE)

# preparing data for modeling
dt.53S.B <- prepare_data(mat = web.53S.B, sampl.eff = rep(20, nrow(web.53S.B)))

# model
set.seed(2)
options(mc.cores = 4)
fit.53S.B <- fit_model(dt.53S.B, refresh = 0, model = "Young2021", iter_sampling = 2000,
                       iter_warmup = 2000, beta = 0.1)

# checking model
BayesianWebs::check_model(fit.53S.B, data = dt.53S.B)

# posteriors
post.53S.B <- get_posterior(fit.53S.B, data = dt.53S.B)

# plotting interaction probability
plot_interaction_prob(post.53S.B)

# getting predicted counts
pred.df.53S.B <- predict_counts(fit.53S.B, data = dt.53S.B)

# plotting counts
plot_counts_pred(pred.df.53S.B, sort = FALSE)

# plotting predicted counts against observed counts
plot_counts_pred_obs(pred.df.53S.B, data = dt.53S.B)

# plotting residuals
plot_residuals(pred.df.53S.B, data = dt.53S.B, sort = FALSE)

# getting predicted counts, rounding
counts.53S.B <- pred.df.53S.B %>%
  mutate(interaction = paste(Plant, Animal, sep = "-")) %>%
  dplyr::group_by(interaction) %>%
  dplyr::summarize(pred.count = round(mean(count))) %>%
  filter(pred.count > 0)





#### example WITH APIS ####
web.54S.W <- as.matrix(webs$`54S.W`)

# plotting
plot_counts_obs(web.54S.W, sort = FALSE)

# preparing data for modeling
dt.54S.W <- prepare_data(mat = web.54S.W, sampl.eff = rep(20, nrow(web.54S.W)))

# model
set.seed(2)
options(mc.cores = 4)
fit.54S.W <- fit_model(dt.54S.W, refresh = 0, model = "Young2021", iter_sampling = 2000,
                 iter_warmup = 2000, beta = 0.005)

# checking model
BayesianWebs::check_model(fit.54S.W, data = dt.54S.W)
plot_prior(beta = 0.005, fit = fit.54S.W, data = dt.54S.W)
# posteriors
post.54S.W <- get_posterior(fit.54S.W, data = dt.54S.W)

# plotting interaction probability
plot_interaction_prob(post.54S.W)

# getting predicted counts
pred.df.54S.W <- predict_counts(fit.54S.W, data = dt.54S.W)

# plotting counts
plot_counts_pred(pred.df.54S.W, sort = FALSE)

# plotting predicted counts against observed counts
plot_counts_pred_obs(pred.df.54S.W, data = dt.54S.W)

# plotting residuals
plot_residuals(pred.df.54S.W, data = dt.54S.W, sort = FALSE)

# getting predicted counts, rounding
counts.54S.W <- pred.df.54S.W %>%
  mutate(interaction = paste(Plant, Animal, sep = "-")) %>%
  dplyr::group_by(interaction) %>%
  dplyr::summarize(pred.count = round(mean(count))) %>%
  filter(pred.count > 0)
  
