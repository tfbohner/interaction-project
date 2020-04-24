# Interaction project
## 3. Post-processing
#### Purpose: Match correlations with distance relationships
#### Author: Teresa Bohner
#### Date Modified: 21 April 2020

## Inputs: tree correlations, field data (distance)
## Outputs: matched tree correlation dataframe

## load packages
library(tidyverse)
library(modelr)
library(brms)
library(tidybayes)

## load and clean necessary data
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)

dist <- read.csv("Processed Data/dist_focal.csv") %>% 
  mutate(tree1=str_remove(focalID, "_"),
         tree2=str_remove(periphID, "_"))

tree_cor <- read.csv("Processed Data/synchrony_multivar.csv")

## wrangle data----
tree_site <- dplyr::select(rings_field, c(tree.uniqueID,Region, hilo)) %>% 
  group_by(tree.uniqueID) %>% 
  summarise_all(first) %>% 
  ungroup() %>% 
  dplyr::select(-Species)

dist_sub <- dist %>% 
  group_by(tree1, tree2)%>% 
  summarise_all(first)


matched_cor <- inner_join(dist_sub, tree_cor, by=c("tree1", "tree2"))


write.csv(matched_cor, "Processed Data/synchrony_distance.csv", row.names = F)




