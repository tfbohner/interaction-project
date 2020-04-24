# Interaction project
## 2. Analyze
#### Purpose: Model effects of distance, site characteristics on pairwise synchrony
#### Author: Teresa Bohner
#### Date Modified: 21 April 2020

## load packages
library(tidyverse)
library(brms)
library(tidybayes)

## load necessary data
dist_corr <- read.csv("Processed Data/synchrony_distance.csv")

rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)

tree_data <- rings_field %>% 
  group_by(tree.uniqueID) %>% 
  summarise(DBH=first(DBH))

spp_diff <- read.csv("Processed Data/pred_spp_diff.csv")
site_grow <- read.csv("Processed Data/pred_site_grow.csv") %>% 
  dplyr::select(-c(DBH, totBA))

## join all data ----
alldata <- dist_corr %>% 
  left_join(tree_data, by=c("periphID"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x)) %>% 
  left_join(spp_diff, by="pair") %>% 
  left_join(site_grow, by="Site")

## brms models----
pair_mod <- brm(.value ~ dist + sizeratio + abs(spp_diff) + BAI_pred, data=alldata, cores=4)
pair_mod2 <- brm(.value ~ dist + sizeratio + pair + BAI_pred, data=alldata, cores=4)
pair_mod3 <- brm(.value ~ dist*pair + sizeratio + pair + BAI_pred, data=alldata, cores=4)

l2 <- loo(pair_mod2)
l3 <- loo(pair_mod3)

loo_model_weights(l2, l3)
