librarian::shelf(tidyverse, sf)

pollinator <- read.csv(file = file.path("data", "L1_wrangled", "cleaned-SRS-plant-pollinator.csv"))
pollinator <- pollinator %>%
  filter(location != "20B")
# get # of interactions per grid cell per patch type
interactions <- pollinator %>%
  dplyr::count(patch, location) %>%
  pivot_wider(names_from = patch, values_from = n)

# interaction_c <- pollinator %>%
#   filter(patch == "B") %>%
#   dplyr::count(location)
# 
# interaction_w <- pollinator %>%
#   filter(patch == "W") %>%
#   dplyr::count(location)

# making 100x100 patch
site <- st_polygon(list(rbind(
  c(0, 0),
  c(100, 0),
  c(100, 100),
  c(0, 100),
  c(0, 0)
))) |> st_sfc()

# adding grids
grid <- st_make_grid(site, cellsize = 12.5, square = TRUE)
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

# plotting to figure out how grids are labeled
grid_sf %>%
  ggplot() +
  geom_sf(fill = "grey80") +
  geom_sf(
    data = dplyr::filter(grid_sf, grid_id == 9),
    fill = NA,
    color = "red",
    linewidth = 1.2
  ) +
  theme_minimal()

# vector of grid cell names
grid_id <- c("13C", "13D", "14C", "14D", "15C", "15D", "16C", "16D", 
             "13A", "13B", "14A", "14B", "15A", "15B", "16A", "16B",
             "9C", "9D", "10C", "10D", "11C", "11D", "12C", "12D",
             "9A", "9B", "10A", "10B", "11A", "11B", "12A", "12B",
             "5C", "5D", "6C", "6D", "7C", "7D", "8C", "8D",
             "5A", "5B", "6A", "6B", "7A", "7B", "8A", "8B",
             "1C", "1D", "2C", "2D", "3C", "3D", "4C", "4D",
             "1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B")

# adding grid names
grid_sf$grid_names <- grid_id

# joining to interaction data
grid_interactions <- grid_sf %>%
  left_join(interactions, by = c("grid_names" = "location"))%>%
  pivot_longer(cols = c("B", "W"), names_to = "patch", values_to = "interactions")
# grid_interactions_c <- grid_sf %>%
#   left_join(interaction_c, by = c("grid_names" = "location")) 
# 
# grid_interactions_w <- grid_sf %>%
#   left_join(interaction_w, by = c("grid_names" = "location")) 
  
# plotting
grid_interactions %>%
  ggplot() +
  geom_sf(aes(fill = interactions)) +
  scale_fill_viridis_c(name = "Interactions") +
  theme_minimal() +
  facet_grid(~patch)

# 
# grid_interactions_c %>%
#   ggplot() +
#   geom_sf(aes(fill = n)) +
#   scale_fill_viridis_c(name = "Interactions") +
#   theme_minimal() +
#   labs(
#     title = "Connected Patches"
#   )
# 
# 
# grid_interactions_w %>%
#   ggplot() +
#   geom_sf(aes(fill = n)) +
#   scale_fill_viridis_c(name = "Interactions") +
#   theme_minimal() +
#   labs(
#     title = "Winged Patches"
#   )








