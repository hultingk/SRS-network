# -------------------------------------- #
#### Analysis script ####
# Author: Katherine Hulting, hultingk@msu.edu
# -------------------------------------- #

# loading libraries 
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggeffects)
library(ggpubr)

# load data
network_metrics <- read.csv(file = file.path("data", "network_metrics.csv"))


#### nestedness ####
m.nestedness <- glmmTMB(nestedness ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.nestedness)
plot(simulateResiduals(m.nestedness))
check_model(m.nestedness)

#### H2' ####
m.h2 <- glmmTMB(H2 ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.h2)
plot(simulateResiduals(m.h2))
check_model(m.h2)

#### interaction diversity ####
m.interaction.div <- glmmTMB(Shannon.diversity ~ patch + (1|block),
                data = network_metrics,
                family = "gaussian")
summary(m.interaction.div)
plot(simulateResiduals(m.interaction.div))
check_model(m.interaction.div)

#### interaction diversity ####
m.links <- glmmTMB(links.per.species ~ patch + (1|block),
                             data = network_metrics,
                             family = "gaussian")
summary(m.links)
plot(simulateResiduals(m.links))
check_model(m.links)


#### flower diversity ####
# with Apis
m.floral.div <- glmmTMB(floral_diversity ~ patch + (1|block),
                   data = network_metrics,
                   family = "gaussian")
summary(m.floral.div)
plot(simulateResiduals(m.floral.div))
check_model(m.floral.div)

# no Apis
m.floral.div <- glmmTMB(floral_div_noApis ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.floral.div)
plot(simulateResiduals(m.floral.div))
check_model(m.floral.div)


#### pollinator diversity ####
m.pollinator.div <- glmmTMB(pollinator_diversity ~ patch + (1|block),
                        data = network_metrics,
                        family = "gaussian")
summary(m.pollinator.div)
plot(simulateResiduals(m.pollinator.div))
check_model(m.pollinator.div)

# no Apis
m.pollinator.div <- glmmTMB(pollinator_div_noApis ~ patch + (1|block),
                            data = network_metrics,
                            family = "gaussian")
summary(m.pollinator.div)
plot(simulateResiduals(m.pollinator.div))
check_model(m.pollinator.div)