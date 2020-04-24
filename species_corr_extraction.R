# Interaction project
## 3. Post-processing
#### Purpose: Intra and interspecific correlation comparisons
#### Author: Teresa Bohner
#### Date Modified: 20 April 2020

## Inputs: tree correlation exports both time-invariant and decadal
## Outputs: Focal species expanded dataframe, intra/inter relationship dataframe

## load packages
library(tidyverse)
library(brms)
library(tidybayes)

## load necessary data
bigvars2 <- read.csv("Processed Data/synchrony_multivar.csv")

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


unique(pjdat$tree1[which(pjdat$comp=="intra")]) %in% unique(pjdat$tree2[which(pjdat$comp=="intra")])
## expand to intra-interspecific correlation comparisons----
acgrid <- expand_grid(acdat$.variable[which(acdat$comp=="intra")], acdat$.variable[which(acdat$comp=="inter")]) %>% 
  rename(var1=1, var2=2) %>% 
  inner_join(dplyr::select(acdat, c("var1"=".variable", ".value", "Site", "Neighborhood", "intra_sp"="spp1", "tree1", "pair"))) %>% 
  inner_join(dplyr::select(acdat, c(".variable", ".value", "Site", "Neighborhood", "inter_sp"="spp2", "tree1", "pair")), by=c("var2"=".variable")) %>% 
  filter(Site.x==Site.y, Neighborhood.x==Neighborhood.y, tree1.x==tree1.y) %>% 
  mutate(diff=.value.x-.value.y,
         intra_greater=ifelse(diff>0, 1,0))

length(acgrid$intra_greater[which(acgrid$intra_greater==1)])/nrow(acgrid)

ac_summ <- acgrid %>% 
  group_by(var1, Site.x, Neighborhood.x, pair.y, intra_sp, inter_sp) %>% 
  summarise(intra_greater=sum(intra_greater),
            size=length(intra_greater)) %>% 
  mutate(nb=as.factor(Neighborhood.x))


  
# test_mod <- brm(intra_greater|trials(size)~ 1, data=ac_summ, family=binomial, cores=4)
# 
# test_mod2 <- brm(intra_greater|trials(size)~ Site.x +pair.y, data=ac_summ, family=binomial, cores=4)


pjgrid <- expand_grid(pjdat$.variable[which(pjdat$comp=="intra")], pjdat$.variable[which(pjdat$comp=="inter")]) %>% 
  rename(var1=1, var2=2) %>% 
  inner_join(dplyr::select(pjdat, c("var1"=".variable", ".value", "Site", "Neighborhood", "intra_sp"="spp1","tree1", "pair"))) %>% 
  inner_join(dplyr::select(pjdat, c(".variable", ".value", "Site", "Neighborhood", "spp1", "spp2", "tree1", "tree2", "pair")), by=c("var2"=".variable")) %>% 
  mutate(inter_sp=ifelse(spp1=="pj", as.character(spp2), as.character(spp1))) %>% 
  filter(Site.x==Site.y, Neighborhood.x==Neighborhood.y, 
         as.character(tree1.x)==as.character(tree1.y)|as.character(tree1.x)==as.character(tree2)) %>% 
  mutate(diff=.value.x-.value.y,
         intra_greater=ifelse(diff>0, 1,0)) 

plgrid <- expand_grid(pldat$.variable[which(pldat$comp=="intra")], pldat$.variable[which(pldat$comp=="inter")]) %>% 
  rename(var1=1, var2=2) %>% 
  inner_join(dplyr::select(pldat, c("var1"=".variable", ".value", "Site", "Neighborhood", "intra_sp"="spp1","tree1", "pair"))) %>% 
  inner_join(dplyr::select(pldat, c(".variable", ".value", "Site", "Neighborhood", "spp1", "spp2", "tree1", "tree2", "pair")), by=c("var2"=".variable")) %>% 
  mutate(inter_sp=ifelse(spp1=="pj", as.character(spp2), as.character(spp1))) %>% 
  filter(Site.x==Site.y, Neighborhood.x==Neighborhood.y, 
         as.character(tree1.x)==as.character(tree1.y)|as.character(tree1.x)==as.character(tree2)) %>% 
  mutate(diff=.value.x-.value.y,
         intra_greater=ifelse(diff>0, 1,0))

ppgrid <- expand_grid(ppdat$.variable[which(ppdat$comp=="intra")], ppdat$.variable[which(ppdat$comp=="inter")]) %>% 
  rename(var1=1, var2=2) %>% 
  inner_join(dplyr::select(ppdat, c("var1"=".variable", ".value", "Site", "Neighborhood", "intra_sp"="spp1","tree1", "pair"))) %>% 
  inner_join(dplyr::select(ppdat, c(".variable", ".value", "Site", "Neighborhood", "spp1", "spp2", "tree1", "tree2", "pair")), by=c("var2"=".variable")) %>% 
  mutate(inter_sp=ifelse(spp1=="pj", as.character(spp2), as.character(spp1))) %>% 
  filter(Site.x==Site.y, Neighborhood.x==Neighborhood.y, 
         as.character(tree1.x)==as.character(tree1.y)|as.character(tree1.x)==as.character(tree2)) %>% 
  mutate(diff=.value.x-.value.y,
         intra_greater=ifelse(diff>0, 1,0))


allgrid <- bind_rows(list(acgrid, pjgrid, plgrid, ppgrid)) %>% 
  dplyr::select(var1, var2, intra_cor=.value.x, Site=Site.x, Neighborhood=Neighborhood.x, 
                intra_sp, tree1=tree1.x, inter_cor=.value.y, inter_sp, intra_greater)


intra_inter_summ <- allgrid %>% 
  mutate(Neighborhood=as.factor(Neighborhood)) %>% 
  group_by(tree1, Site, Neighborhood, intra_sp, inter_sp) %>% 
  summarise(med_intra_corr=median(intra_cor),
            med_inter_corr=median(inter_cor),
            intra_greater=sum(intra_greater),
            trials=length(intra_sp),
            prop=intra_greater/trials)

test_mod <- brm(intra_greater|trials(trials)~ 1, data=intra_inter_summ, family=binomial, cores=4)
test_mod2 <- brm(intra_greater|trials(trials)~ intra_sp + Site, data=intra_inter_summ, family=binomial, cores=4)
marginal_effects(test_mod2)

test_mod23 <- brm(intra_greater|trials(trials)~ (1|intra_sp) + (1|Site), data=intra_inter_summ, family=binomial, cores=4)



