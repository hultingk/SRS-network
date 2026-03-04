# -------------------------------------- #
#### Data cleaning script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries
librarian::shelf(tidyverse)


# loading data 
floral <- read.csv(file = file.path("data", "L0_original", "2024-SRS-floral-plot.csv"))
floral_count <- read.csv(file = file.path("data", "L0_original", "2024-SRS-floral-species-count.csv"))


# unifying species names and fixing typos
floral <- floral %>%
  mutate(flower_species = dplyr::case_when(
    flower_species %in% c("Desmodium paniculatum", "Desmodium sp.", "Desmodium strictum") ~ "Desmodium sp.",
    flower_species %in% c("Dyschoriste linearis", "Dyschoriste oblongifolia") ~ "Dyschoriste sp.",
    flower_species %in% c("Eupatorium glaucescens", "Eupatorium linearifolium") ~ "Eupatorium linearifolium",
    flower_species %in% c("Sericocarpus asteroides", "Sericocarpus tortifolius") ~ "Sericocarpus sp.",
    flower_species %in% c("Lespedeza angustifolia", "Lespedeza hirta",
                          "Lespedeza repens", "Lespedeza stuevei", "Lespedeza stuvei", "Lespedeza virginica", 
                          "Lespedeza cuneata", "Lespedeza violacea") ~ "Lespedeza sp.",
    flower_species %in% c("Tephrosia florida", "Tephrosia spicata") ~ "Tephrosia sp.",
    flower_species %in% c("Solidago nemoralis", "Solidago odora") ~ "Solidago sp.",
    flower_species %in% c("Liatris elegans", "Liatris graminifolia", "Liatris secunda", "Liatris tenuifolia",
                          "Liatris virgata", "Liatris squarrulosa") ~ "Liatris sp.",
    flower_species %in% c("Helianthus hirsutus") ~ "Helianthus sp.",
    flower_species %in% c("Galactia sp.", "Galactia regularis", "Galatica sp.") ~ "Galactia sp.",
    flower_species %in% c("Galium Sp.", "Galium sp.") ~ "Galium sp.",
    .default = flower_species
  ))

# unifying species names and fixing typos
floral_count <- floral_count %>%
  mutate(flower_species = dplyr::case_when(
    flower_species %in% c("Desmodium paniculatum", "Desmodium sp.", "Desmodium strictum") ~ "Desmodium sp.",
    flower_species %in% c("Dyschoriste linearis", "Dyschoriste oblongifolia") ~ "Dyschoriste sp.",
    flower_species %in% c("Eupatorium glaucescens", "Eupatorium linearifolium") ~ "Eupatorium linearifolium",
    flower_species %in% c("Sericocarpus asteroides", "Sericocarpus tortifolius") ~ "Sericocarpus sp.",
    flower_species %in% c("Lespedeza angustifolia", "Lespedeza hirta",
                          "Lespedeza repens", "Lespedeza stuevei", "Lespedeza stuvei", "Lespedeza virginica", 
                          "Lespedeza cuneata", "Lespedeza violacea") ~ "Lespedeza sp.",
    flower_species %in% c("Tephrosia florida", "Tephrosia spicata") ~ "Tephrosia sp.",
    flower_species %in% c("Solidago nemoralis", "Solidago odora") ~ "Solidago sp.",
    flower_species %in% c("Liatris elegans", "Liatris graminifolia", "Liatris secunda", "Liatris tenuifolia",
                          "Liatris virgata", "Liatris squarrulosa") ~ "Liatris sp.",
    flower_species %in% c("Helianthus hirsutus") ~ "Helianthus sp.",
    flower_species %in% c("Galactia sp.", "Galactia regularis", "Galatica sp.") ~ "Galactia sp.",
    flower_species %in% c("Galium Sp.", "Galium sp.") ~ "Galium sp.",
    .default = flower_species
  ))

# floral abundance 
floral_avg <- floral_count %>%
  dplyr::mutate(fl_avg = rowMeans(across(avg_1:avg_5), na.rm = TRUE)) %>%
  dplyr::mutate(floral_month = paste(flower_species, sampling_round, sep = "-")) %>%
  dplyr::select(floral_month, fl_avg) %>%
  dplyr::group_by(floral_month) %>%
  dplyr::summarize(fl_avg = mean(fl_avg)) 


floral_patch <- floral %>%
  dplyr::mutate(unique_ID = paste(block, patch, sampling_round, sep = "-")) %>%
  dplyr::group_by(unique_ID, sampling_round, flower_species) %>%
  dplyr::summarize(individuals = sum(no_individuals)) %>%
  dplyr::filter(flower_species != 0) %>%
  dplyr::mutate(floral_month = paste(flower_species, sampling_round, sep = "-"))


floral_patch <- floral_patch %>%
  left_join(floral_avg, by = "floral_month") %>%
  dplyr::mutate(fl_avg = if_else(is.na(fl_avg), individuals, fl_avg))

floral_patch <- floral_patch %>%
  dplyr::group_by(unique_ID) %>%
  dplyr::summarize(fl_abund = sum(fl_avg)) %>%
  separate(unique_ID, into = c("block", "patch", "sampling_round"))

# log transforming floral abundance
floral_patch$fl_abund_log <- log(floral_patch$fl_abund)

# floral abundance model
m.floral_abund <- glmmTMB(fl_abund_log ~ patch * sampling_round + (1|block),
              data = floral_patch)
summary(m.floral_abund)

# plotting
floral_patch$sampling_round <- factor(floral_patch$sampling_round, levels = c("June", "July", "August", "September"))
floral_abund.plot <- floral_patch %>%
  ggplot() + 
  geom_boxplot(aes(sampling_round, fl_abund_log, fill = patch)) +
  theme_classic(base_size = 16) +
  scale_fill_manual(values=c("#F5097C","#F7B3D4"), labels = c("Connected", "Unconnected"), name = "Patch Type") +
  xlab("Month") +
  ylab("Log floral abundance")
floral_abund.plot


# exporting
# pdf(file = file.path("plots", "floral_abundance.pdf"), width = 6, height = 4)
# floral_abund.plot
# dev.off()
