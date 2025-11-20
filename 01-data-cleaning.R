# -------------------------------------- #
#### Data cleaning script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse)


# loading data 
pollinator <- read.csv(file = file.path("data", "2024-SRS-plant-pollinator.csv"))

# --------------------------- #
#### pollinator species cleaning ####
# --------------------------- #
# removing wasps, lost insects, unknowns, no observations
pollinator <- pollinator %>%
  filter(!pollinator_common %in% c("Wasp", "unknown butterfly", "unknown skipper", "No associated insect", "0", "WASP", "LOST", "beetle", "NOT A POLLINATOR", "Butterfly", "", " "))

# grouping Erynnis species 
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis horatius", "Erynnis sp.")
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis zarucco", "Erynnis sp.")

pollinator <- pollinator %>%
  mutate(pollinator_species = dplyr::case_when(
    pollinator_species %in% c("Zodion sp.1", "Zodion sp.2") ~ "Zodion sp.",
    .default = pollinator_species
  )) %>%
  filter(pollinator_species != "Lasioglossum sp.")


# adding family and order 
pollinator <- pollinator %>%
  separate(pollinator_species, into = c("genus", "species"), sep = " ") %>%
  mutate(family = dplyr::case_when(
    genus %in% c("Abaeis", "Eurema", "Phoebis", "Pyrisitia") ~ "Pieridae",
    genus %in% c("Agapostemon", "Lasioglossum", "Augochlora", "Augochlorella", "Augochloropsis", "Halictus") ~ "Halictidae",
    genus %in% c("Agraulis", "Danaus", "Euptoieta", "Junonia", "Phyciodes") ~ "Nymphalidae",
    genus %in% c("Allograpta", "Copestylum", "Milesia", "Ocyptamus", "Palpada", "Platycheirus", 
                 "Syritta", "Toxomerus", "Xylota") ~ "Syrphidae",
    genus %in% c("Anthrax", "Anastoechus", "Chrysanthrax", "Exoprosopa", "Geron", "Lepidophora", 
                 "Poecilognathus", "Systoechus", "Villa") ~ "Bombyliidae",
    genus %in% c("Ancyloxypha", "Burnsius", "Cecropterus", "Epargyreus", "Erynnis", "Hylephila",
                 "Lerema", "Nastra", "Panoquina", "Polites", "Problema", "Thorybes", "Urbanus", "Wallengrenia") ~ "Hesperiidae",
    genus %in% c("Anthidiellum", "Anthidium", "Coelioxys", "Megachile", "Trachusa", "Hoplitis") ~ "Megachilidae",
    genus %in% c("Apis", "Bombus", "Epimelissodes", "Melissodes", "Xylocopa", "Ceratina") ~ "Apidae",
    genus %in% c("Archytas", "Tachinidae") ~ "Tachinidae",
    genus %in% c("Atlides", "Calycopis", "Cupido", "Satyrium", "Strymon") ~ "Lycaenidae",
    genus %in% c("Battus", "Papilio") ~ "Papilionidae",
    genus %in% c("Cochliomyia", "Lucilia") ~ "Calliphoridae",
    genus %in% c("Colletes", "Hylaeus") ~ "Colletidae",
    genus %in% c("Physoconops", "Zodion") ~ "Conopidae",
    genus %in% c("Tabanidae") ~ "Tabanidae",
    genus %in% c("Calyptratae") ~ "Calyptratae", 
    genus %in% c("Perdita") ~ "Andrenidae",
    genus %in% c("Rivellia") ~ "Platystomatidae",
    genus %in% c("Dolichopodidae") ~ "Dolichopodidae",
    genus %in% c("Miltogramminae") ~ "Sarcophagidae"
  )) %>% 
  mutate(order = dplyr::case_when(
    family %in% c("Apidae", "Colletidae", "Halictidae", "Megachilidae", "Andrenidae") ~ "Hymenoptera",
    family %in% c("Bombyliidae", "Calliphoridae", "Conopidae", "Syrphidae", "Tabanidae", "Tachinidae", "Calyptratae", "Platystomatidae", "Sarcophagidae", "Dolichopodidae") ~ "Diptera",
    family %in% c("Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae") ~ "Lepidoptera"
  ))


pollinator <- pollinator %>%
  mutate(pollinator_species = paste(genus, species, sep = " "))


# --------------------------- #
#### flower species cleaning ####
# --------------------------- #
# combining Eupatorium glaucescens and Eupatorium linearifolium as one species
pollinator$flower_species <- str_replace(pollinator$flower_species, "Eupatorium glaucescens", "Eupatorium linearifolium")
pollinator <- pollinator %>%
  mutate(flower_species = dplyr::case_when(
    flower_species %in% c("Desmodium 1", "Desmodium 2", "Desmodium ciliare",
                          "Desmodium fernaldii", "Desmodium marilandicum",
                          "Desmodium obtusum", "Desmodium viridiflorum") ~ "Desmodium sp.",
    flower_species %in% c("Sericocarpus asteroides", "Sericocarpus tortifolius") ~ "Sericocarpus sp.",
    flower_species %in% c("Lespedeza angustifolia", "Lespedeza hirta",
                          "Lespedeza repens", "Lespedeza stuevei", "Lespedeza virginica", "Lespedeza cuneata") ~ "Lespedeza sp.",
    flower_species %in% c("Tephrosia florida", "Tephrosia spicata") ~ "Tephrosia sp.",
    flower_species %in% c("Solidago nemoralis", "Solidago odora") ~ "Solidago sp.",
    .default = flower_species
  ))





# --------------------------- #
#### writing cleaned file ####
# --------------------------- #
write.csv(pollinator, file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


