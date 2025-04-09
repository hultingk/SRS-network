# -------------------------------------- #
#### Data cleaning script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libaries
library(tidyverse)

# loading data 
pollinator <- read.csv(file = file.path("data", "2024-SRS-plant-pollinator.csv"))

# removing wasps, lost insects, unknowns, no observations
pollinator <- pollinator %>%
  filter(!pollinator_common %in% c("Wasp", "unknown butterfly", "unknown skipper", "No associated insect", "0", "WASP", "LOST", "beetle", "NOT A POLLINATOR", "Butterfly"))


pollinator <- pollinator %>%
  mutate(family = dplyr::case_when(
    pollinator_species %in% c("Abaeis nicippe", "Eurema nicippe", "Phoebis sennae", 
                              "Pyrisitia lisa") ~ "Pieridae",
    pollinator_species %in% c("Agapostemon splendens", "Lasioglossum sp.", "Augochlora pura",
                              "Augochlorella aurata", "Augochloropsis metallica", "Halictus sp.") ~ "Halictidae",
    pollinator_species %in% c("Agraulis vanillae", "Danaus plexippus", "Euptoieta claudia",
                              "Junonia coenia", "Phyciodes tharos") ~ "Nymphalidae",
    pollinator_species %in% c("Allograpta exotica", "Allograpta obliqua", "Copestylum vittatum",
                              "Milesia virginiensis", "Ocyptamus fuscipennis", "Palpada furcata",
                              "Platycheirus sp.", "Syritta flaviventris", "Toxomerus geminatus",
                              "Toxomerus marginatus", "Xylota sp.") ~ "Syrphidae",
    pollinator_species %in% c("Agraulis vanillae", "Anthrax argyropygus", "Anastoechus barbatus",
                              "Anthrax irroratus", "Chrysanthrax cypris", "Exoprosopa fasciata",
                              "Exoprosopa fascipennis", "Geron sp.", "Lepidophora lepidocera",
                              "Poecilognathus sulphureus", "Systoechus solitus", "Villa lateralis") ~ "Bombyliidae",
    pollinator_species %in% c("Ancyloxypha numitor", "Burnsius albezens", "Cecropterus lyciades",
                              "Epargyreus clarus", "Erynnis horatius", "Erynnis sp.", "Erynnis zarucco",
                              "Hylephila phyleus", "Lerema accius", "Nastra lherminier",
                              "Panoquina ocola", "Polites vibex", "Problema byssus", "Thorybes bathyllus",
                              "Urbanus proteus", "Wallengrenia otho") ~ "Hesperiidae",
    pollinator_species %in% c("Anthidiellum notatum", "Anthidiellum perplexum",
                              "Anthidium maculifrons", "Coelioxys sp.", "Megachile albitarsis",
                              "Megachile deflexa", "Megachile georgica", "Megachile mendica",
                              "Megachile petulans", "Megachile sp.", "Megachile texana", 
                              "Megachile xylocopoides", "Trachusa ridingsii", "Hoplitis sp.") ~ "Megachilidae",
    pollinator_species %in% c("Apis mellifera", "Bombus fraternus", "Bombus griseocollis",
                              "Bombus impatiens", "Bombus pensylvanicus", "Epimelissodes atripes",
                              "Epimelissodes obliqua", "Melissodes bimaculatus", "Melissodes boltoniae",
                              "Melissodes communis", "Melissodes sp.", "Xylocopa micans", 
                              "Xylocopa virginica", "Ceratina sp.") ~ "Apidae",
    pollinator_species %in% c("Archytas metallicus", "Tachinidae sp.1") ~ "Tachinidae",
    pollinator_species %in% c("Atlides halesus", "Calycopis cecrops", "Cupido comyntas", 
                              "Satyrium calanus", "Satyrium titus", "Strymon melinus") ~ "Lycaenidae",
    pollinator_species %in% c("Battus philenor", "Papilio glaucus", "Papilio palamedes",
                              "Papilio troilus") ~ "Papilionidae",
    pollinator_species %in% c("Cochliomyia macellaria", "Lucilia sp.") ~ "Calliphoridae",
    pollinator_species %in% c("Colletes sp.", "Hylaeus sp.") ~ "Colletidae",
    pollinator_species %in% c("Physoconops bulbirostris", "Zodion sp.1", "Zodion sp.2") ~ "Conopidae",
    pollinator_species %in% c("Tabanidae sp.") ~ "Tabanidae",
    pollinator_species %in% c("Calyptratae family") ~ "Calyptratae", 
    pollinator_species %in% c("Perdita sp.") ~ "Andrenidae",
    pollinator_species %in% c("Rivellia sp.") ~ "Platystomatidae"
  )) %>% 
  mutate(order = dplyr::case_when(
    family %in% c("Apidae", "Colletidae", "Halictidae", "Megachilidae", "Andrenidae") ~ "Hymenoptera",
    family %in% c("Bombyliidae", "Calliphoridae", "Conopidae", "Syrphidae", "Tabanidae", "Tachinidae", "Calyptratae", "Platystomatidae") ~ "Diptera",
    family %in% c("Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae") ~ "Lepidoptera"
  ))


# --------------------------- #
#### flower species cleaning ####
# --------------------------- #
# combining Eupatorium glaucescens and Eupatorium linearifolium as one species
pollinator$flower_species <- str_replace(pollinator$flower_species, "Eupatorium glaucescens", "Eupatorium linearifolium")

# checking flower species
pollinator %>%
  dplyr::count(flower_species)

# --------------------------- #
#### pollinator species cleaning ####
# --------------------------- #
# grouping Erynnis species 
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis horatius", "Erynnis sp.")
pollinator$pollinator_species <- str_replace(pollinator$pollinator_species, "Erynnis zarucco", "Erynnis sp.")

# checking pollinator species
pollinator %>%
  dplyr::count(pollinator_species)

pollinator <- pollinator %>%
  mutate(pollinator_analysis = case_when(
    pollinator_species %in% c("Megachile albitarsis", "Megachile deflexa", "Megachile georgica",
                              "Megachile mendica", "Megachile petulans", "Megachile texana", "Megachile xylocopoides") ~ "Megachile sp.",
    pollinator_species %in% c("Melissodes bimaculatus", "Melissodes boltoniae", "Melissodes communis") ~ "Melissodes sp.",
    .default = pollinator_species
  ))

# --------------------------- #
#### writing cleaned file ####
# --------------------------- #
write.csv(pollinator, file = file.path("data", "cleaned-SRS-plant-pollinator.csv"))


