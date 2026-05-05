librarian::shelf(glmmTMB, DHARMa, tidyverse, performance, ggeffects)

network_metrics_boot <- read.csv(file = file.path("data", "L2_boot_metrics", "network_metrics_boot.csv"))

# nestedness
m_nodf <- glmmTMB(nodf ~ patch + (1|block),
                  data = network_metrics_boot)
summary(m_nodf)

# links per species 
m_links <- glmmTMB(links_per_sp ~ patch + (1|block),
                   data = network_metrics_boot)
summary(m_links)

# pollinator niche overlap
m_HL_niche <- glmmTMB(HL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_HL_niche)

# plant niche overlap
m_LL_niche <- glmmTMB(LL_niche ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_LL_niche)

# web asymmetry
m_asymmetry <- glmmTMB(asymmetry ~ patch + (1|block),
                      data = network_metrics_boot)
summary(m_asymmetry)

# connectance
m_connectance <- glmmTMB(connectance ~ patch + (1|block),
                       data = network_metrics_boot)
summary(m_connectance)

# linkage density
m_density <- glmmTMB(linkage_density ~ patch + (1|block),
                         data = network_metrics_boot)
summary(m_density)



