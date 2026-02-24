# -------------------------------------- #
#### Data cleaning script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse)


# loading data 
pollinator <- read.csv(file = file.path("data", "L0_original", "2024-SRS-plant-pollinator.csv"))

# --------------------------- #
#### pollinator species cleaning ####
# --------------------------- #
# removing wasps, lost insects, unknowns, no observations
pollinator <- pollinator %>%
  filter(!pollinator_common %in% c("Wasp", "unknown butterfly", "unknown skipper", "No associated insect", 
                                   "0", "WASP", "LOST", "beetle", "NOT A POLLINATOR", "Butterfly", "", " ",
                                   "Mosquito", "damaged - missing from pin", "not a pollinator"))

# grouping some sets of species 
pollinator <- pollinator %>%
  mutate(pollinator_species = dplyr::case_when(
    pollinator_species %in% c("Erynnis horatius", "Erynnis zarucco") ~ "Erynnis sp.", # grouping Erynnis species 
    pollinator_species %in% c("Melissodes communis", "Melissodes sp.") ~ "Melissodes sp.01",
    .default = pollinator_species
  )) %>%
  filter(!pollinator_species %in% c("Lasioglossum sp.", "Skipper sp.", "Lasioglossum batya")) # not confident on my ID of Lasioglossum batya

# adding family and order 
pollinator <- pollinator %>%
  separate(pollinator_species, into = c("genus", "species"), sep = " ", remove = F) %>%
  mutate(family = dplyr::case_when(
    genus %in% c("Eurema", "Phoebis") ~ "Pieridae",
    genus %in% c("Agapostemon", "Lasioglossum", "Augochlora", "Augochlorella", "Augochloropsis", "Halictus") ~ "Halictidae",
    genus %in% c("Agraulis", "Danaus", "Euptoieta", "Junonia", "Phyciodes") ~ "Nymphalidae",
    genus %in% c("Allograpta", "Copestylum", "Milesia", "Ocyptamus", "Palpada", "Platycheirus", 
                 "Syritta", "Toxomerus", "Xylota", "Orthonevra", "Syrphidae") ~ "Syrphidae",
    genus %in% c("Anthrax", "Anastoechus", "Chrysanthrax", "Exoprosopa", "Geron", "Lepidophora", 
                 "Poecilognathus", "Systoechus", "Villa", "Toxophora") ~ "Bombyliidae",
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
    genus %in% c("Chloropidae") ~ "Chloropidae",
    genus %in% c("Tabanidae") ~ "Tabanidae",
    genus %in% c("Calyptratae") ~ "Calyptratae", 
    genus %in% c("Perdita") ~ "Andrenidae",
    genus %in% c("Rivellia") ~ "Platystomatidae",
    genus %in% c("Dolichopodidae") ~ "Dolichopodidae",
    genus %in% c("Miltogramminae") ~ "Sarcophagidae",
    genus %in% c("Acalyptrate") ~ "Acalyptrate",
    genus %in% c("Milichiidae") ~ "Milichiidae"
  )) %>% 
  mutate(order = dplyr::case_when(
    family %in% c("Apidae", "Colletidae", "Halictidae", "Megachilidae", "Andrenidae") ~ "Hymenoptera",
    family %in% c("Bombyliidae", "Calliphoridae", "Conopidae", "Chloropidae", "Syrphidae", "Tabanidae", "Tachinidae", "Calyptratae", "Platystomatidae", "Sarcophagidae", "Dolichopodidae", "Acalyptrate", "Milichiidae") ~ "Diptera",
    family %in% c("Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae") ~ "Lepidoptera"
  ))


# removing fly not IDed to family level
pollinator <- pollinator %>%
  filter(!family %in% c("Acalyptrate"))


# --------------------------- #
#### flower species cleaning ####
# --------------------------- #
# combining Eupatorium glaucescens and Eupatorium linearifolium as one species
pollinator <- pollinator %>%
  mutate(flower_species = dplyr::case_when(
    flower_species %in% c("Desmodium 1", "Desmodium 2", "Desmodium ciliare",
                          "Desmodium fernaldii", "Desmodium marilandicum",
                          "Desmodium obtusum", "Desmodium viridiflorum", "Desmodium marilandicum complex") ~ "Desmodium sp.",
    flower_species %in% c("Eupatorium glaucescens", "Eupatorium linearifolium") ~ "Eupatorium linearifolium",
    flower_species %in% c("Sericocarpus asteroides", "Sericocarpus tortifolius") ~ "Sericocarpus sp.",
    flower_species %in% c("Lespedeza angustifolia", "Lespedeza hirta",
                          "Lespedeza repens", "Lespedeza stuevei", "Lespedeza virginica", "Lespedeza cuneata") ~ "Lespedeza sp.",
    flower_species %in% c("Tephrosia florida", "Tephrosia spicata") ~ "Tephrosia sp.",
    flower_species %in% c("Solidago nemoralis", "Solidago odora") ~ "Solidago sp.",
    flower_species %in% c("Liatris elegans", "Liatris graminifolia", "Liatris secunda", "Liatris tenuifolia",
                          "Liatris virgata") ~ "Liatris sp.",
    flower_species %in% c("Helianthus hirsutus") ~ "Helianthus sp.",
    .default = flower_species
  ))


# --------------------------- #
#### writing cleaned file ####
# --------------------------- #
write.csv(pollinator, file = file.path("data", "L1_wrangled","cleaned-SRS-plant-pollinator.csv"))


