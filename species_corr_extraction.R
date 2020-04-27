# Interaction project
## 3. Post-processing
#### Purpose: Intra and interspecific correlation comparisons
#### Author: Teresa Bohner
#### Date Modified: 24 April 2020

## Inputs: tree correlation exports both time-invariant and decadal
## Outputs: Focal species expanded dataframe, intra/inter relationship dataframe

## load packages
library(tidyverse)
library(modelr)
library(brms)
library(tidybayes)

## load necessary data
bigvars2 <- read.csv("Processed Data/synchrony_multivar.csv")

## load tree data 
tree_data <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  mutate(tree.uniqueID=as.factor(str_remove(tree.uniqueID, "_")),
         Species=tolower(Species)) %>% 
  dplyr::select(c(Site, Species, Neighborhood, Region, hilo, tree.uniqueID, DBH)) %>% 
  group_by(tree.uniqueID) %>% 
  summarise_all(first) %>% 
  filter(Site!="ic")

## 1a. Functions for time invariant models----
## subset each species
spp_subsetter <- function(data, sp) {
  a <- filter(data, pair==paste(sort(c(sp, "ac")), collapse="-")) %>%
    mutate(focal_sp=sp,
           comp_sp="ac")
  b <- filter(data, pair==paste(sort(c(sp, "pj")), collapse="-")) %>%
    mutate(focal_sp=sp,
           comp_sp="pj")
  c <- filter(data, pair==paste(sort(c(sp, "pl")), collapse="-")) %>%
    mutate(focal_sp=sp,
           comp_sp="pl")
  d <- filter(data, pair==paste(sort(c(sp, "pp")), collapse="-")) %>%
    mutate(focal_sp=sp,
           comp_sp="pp")
  
  bind_rows(a,b,c,d) %>%
    mutate(comp_sp=ifelse(comp_sp==focal_sp, "self", comp_sp))
}

## subset species----
acdat <- spp_subsetter(bigvars2, "ac")
pjdat <- spp_subsetter(bigvars2, "pj")
pldat <- spp_subsetter(bigvars2, "pl")
ppdat <- spp_subsetter(bigvars2, "pp")

spp_expanded <- bind_rows(acdat, pjdat, pldat, ppdat)

## export spp expanded dataframe----
write.csv(spp_expanded, "Processed Data/synchrony_spp_expanded.csv", row.names = F)


length(unique(acdat$tree1[which(acdat$comp=="intra")])) + 
length(unique(pjdat$tree1[which(pjdat$comp=="intra")]))
length(unique(pldat$tree1[which(pldat$comp=="intra")]))
length(unique(ppdat$tree1[which(ppdat$comp=="intra")]))


## gridding to compare individual pairwise intra versus interspecific correlations ----

## grid all interactions at neighborhood level
tree_grid <- tree_data %>% 
  group_by(tree1=tree.uniqueID, Site, Neighborhood, spp1=Species, Region, hilo) %>%
  data_grid(tree2=tree.uniqueID) %>% 
  left_join(dplyr::select(tree_data, c(tree2=tree.uniqueID, Site2=Site, nb2=Neighborhood, spp2=Species))) %>% 
  filter(Site==Site2, Neighborhood==nb2) %>% 
  mutate(comp=ifelse(spp1==spp2,"intra", "inter")) %>% 
  dplyr::select(-c(Site2, nb2)) %>% 
  filter(tree1!=tree2)


intra_grid <- filter(tree_grid, comp=="intra") %>% 
  rename_at(vars(-c(Site, Neighborhood, Region, hilo)), function(x) str_c(x, ".x"))
inter_grid <- filter(tree_grid, comp=="inter") %>% 
  rename_all(function(x) str_c(x, ".y"))

## grid intra versus interspecific correlations at neighborhood level
pairwise_grid <- expand_grid(intra_grid, inter_grid) %>% 
  filter(Site==Site.y, Neighborhood==Neighborhood.y, tree1.x==tree1.y) %>% 
  dplyr::select(-c(Site.y, Neighborhood.y, Region.y, hilo.y))

## now match each pairwise interaction with the proper correlatios
cor_data <- bigvars2 %>% 
  dplyr::select(c(.variable, .value, .lower, .upper, tree1, tree2))

## match intraspecific correlations
match1 <- inner_join(pairwise_grid, cor_data, by=c("tree1.x"="tree1", "tree2.x"="tree2"))
match2 <- inner_join(pairwise_grid, cor_data, by=c("tree1.x"="tree2", "tree2.x"="tree1"))
match3 <- bind_rows(match1, match2)


## match interspecific correlations
match4 <- inner_join(match3, cor_data, by=c("tree1.y"="tree1", "tree2.y"="tree2"))
match5 <- inner_join(match3, cor_data, by=c("tree1.y"="tree2", "tree2.y"="tree1"))
match6 <- bind_rows(match4, match5)


all_cor <- match6 %>% 
  mutate(med_diff=.value.x-.value.y,
         lower_diff=.lower.x-.lower.y,
         upper_diff=.upper.x-.upper.y,
         intra_greater=ifelse(med_diff>0, 1, 0),
         interval_intra_greater=ifelse(lower_diff>0&upper_diff>0, 1,0)) %>% 
  dplyr::select(-c(spp2.x, spp1.y)) %>% 
  rename(intra_spp=spp1.x, inter_spp=spp2.y)

num_or_char <- function(x)  ifelse(is.numeric(x), mean(x), first(x))

id_cor_summ <- all_cor %>% 
  dplyr::select(-c(tree2.x, tree1.y, tree2.y, .variable.x, .variable.y)) %>% 
  group_by(Region, Site, hilo, Neighborhood, intra_spp, inter_spp, tree1.x) %>% 
  mutate(intra_greater=sum(intra_greater),
         interval_intra_greater=sum(interval_intra_greater),
         trials=length(intra_greater),
         prop=intra_greater/trials) %>% 
  summarise_all(num_or_char)

## export Raa > Rab dataframe----
write.csv(all_cor, "Processed Data/inta_inter_pairwise_all.csv", row.names = F)
write.csv(id_cor_summ, "Processed Data/inta_inter_pairwise_summary.csv", row.names = F)

