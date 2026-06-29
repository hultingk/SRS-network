install.packages("tidyverse")
remotes::install_github("valentinitnelav/bootstrapnet")
librarian::shelf(bootstrapnet, tidyverse)

# load webs
web.matrix <- readRDS(file = file.path("web.matrix.RData"))

# H2'
lst_h2 <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "H2",
                    level = "both",
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_h2, file = file.path("h2.RData"))

# NODF
lst_nodf <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "NODF",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_nodf, file = file.path("nodf.RData"))

# links per species 
lst_links <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "links per species",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_links, file = file.path("links.RData"))


# niche overlap
lst_niche <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "niche overlap",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_niche, file = file.path("niche.RData"))


# robustness
lst_robustness <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "robustness",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_robustness, file = file.path("robustness.RData"))


# shannon diversity of interactions
lst_shannon <- webs.matrix %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower", # column name for plants
                    col_higher = "higher", # column name for insects
                    index = "Shannon diversity",
                    level = "both", 
                    start = 90,
                    step = 1,
                    n_boot = 1000,
                    n_cpu = 4,
                    weighted = F)
saveRDS(lst_shannon, file = file.path("shannon.RData"))

