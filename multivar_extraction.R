# Interaction project
## 3. Post-processing
#### Purpose: Functions and execution of correlation extraction from multivariate brms models
#### Author: Teresa Bohner
#### Date Modified: 20 April 2020

## Inputs: brms multivariate models
## Outputs: tree correlation exports both time-invariant and decadal

## load packages
library(tidyverse)
library(brms)
library(tidybayes)

## load necessary data
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)

## 1a. Functions for time invariant models----
## fun: extract summaries of residual correlations
cor_extract <- function(mod) {
  vars <- mod %>% 
    gather_draws(`rescor.*`, regex=T) %>% 
    median_qi() %>% 
    filter(str_detect(.variable, "pptmm")==FALSE) %>% 
    mutate(
      tree1=str_sub(str_split_fixed(.variable, "__", 3)[,2], start=2),
      spp1=ifelse(str_detect(tree1, "ac")==TRUE, "ac", 
                  ifelse(str_detect(tree1, "pj")==TRUE, "pj", 
                         ifelse(str_detect(tree1, "pl")==TRUE, "pl", 
                                ifelse(str_detect(tree1, "pp")==TRUE, "pp", "am")))),
      tree2=str_sub(str_split_fixed(.variable, "__", 3)[,3], start=2),
      spp2=ifelse(str_detect(tree2, "ac")==TRUE, "ac", 
                  ifelse(str_detect(tree2, "pj")==TRUE, "pj", 
                         ifelse(str_detect(tree2, "pl")==TRUE, "pl", 
                                ifelse(str_detect(tree2, "pp")==TRUE, "pp", "am")))),
      pair=str_c(spp1, spp2, sep="-"),
      comp=ifelse(spp1==spp2, "intra", "inter"),
      Site=i,
      Neighborhood=j) 
}

## fun: extract all draws of residual correlations
keep_draws <- function(mod) {
  vars <- mod %>% 
    gather_draws(`rescor.*`, regex=T) %>% 
    # filter(str_detect(.variable, "pptmm")==FALSE) %>% 
    mutate(
      tree1=str_sub(str_split_fixed(.variable, "__", 3)[,2], start=2),
      spp1=ifelse(str_detect(tree1, "ac")==TRUE, "ac",
                  ifelse(str_detect(tree1, "pj")==TRUE, "pj",
                         ifelse(str_detect(tree1, "pl")==TRUE, "pl",
                                ifelse(str_detect(tree1, "pp")==TRUE, "pp", "am")))),
      tree2=str_sub(str_split_fixed(.variable, "__", 3)[,3], start=2),
      spp2=ifelse(str_detect(tree2, "ac")==TRUE, "ac",
                  ifelse(str_detect(tree2, "pj")==TRUE, "pj",
                         ifelse(str_detect(tree2, "pl")==TRUE, "pl",
                                ifelse(str_detect(tree2, "pp")==TRUE, "pp", "am")))),
      pair=str_c(spp1, spp2, sep="-"),
      comp=ifelse(spp1==spp2, "intra", "inter"),
      Site=i,
      Neighborhood=j)
}


## 1b. Functions for time slice models----
## fun: extract summaries of residual correlations
cor_extract2 <- function(mod) {
  vars <- mod %>% 
    gather_draws(`rescor.*`, regex=T) %>% 
    median_qi() %>% 
    filter(str_detect(.variable, "pptmm")==FALSE) %>% 
    mutate(
      tree1=str_sub(str_split_fixed(.variable, "__", 3)[,2], start=2),
      spp1=ifelse(str_detect(tree1, "ac")==TRUE, "ac", 
                  ifelse(str_detect(tree1, "pj")==TRUE, "pj", 
                         ifelse(str_detect(tree1, "pl")==TRUE, "pl", 
                                ifelse(str_detect(tree1, "pp")==TRUE, "pp", "am")))),
      tree2=str_sub(str_split_fixed(.variable, "__", 3)[,3], start=2),
      spp2=ifelse(str_detect(tree2, "ac")==TRUE, "ac", 
                  ifelse(str_detect(tree2, "pj")==TRUE, "pj", 
                         ifelse(str_detect(tree2, "pl")==TRUE, "pl", 
                                ifelse(str_detect(tree2, "pp")==TRUE, "pp", "am")))),
      pair=str_c(spp1, spp2, sep="-"),
      comp=ifelse(spp1==spp2, "intra", "inter"),
      Site=i,
      Neighborhood=j,
      time=t) 
}


## 2a. Extract correlations from time-invariant models----
## load all of the brms models
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    assign(paste0("fit", i, j), readRDS(paste0("saved models/with covar/fit", i, j, ".rds")), envir = .GlobalEnv)
  }
}

## extract the residual correlations for each model - summary
for(i in sort(unique(rings_field$Site))) {
  for(j in sort(unique(rings_field$Neighborhood))){
    # print(paste0(i,j))
    mod <- get(paste0("fit", i, j))
    vars <- cor_extract(mod)
    
    if(i=="bm"&j==1){
      bigvars <- vars
    } else bigvars <- bind_rows(bigvars, vars)
  }
}

bigvars2 <- filter(bigvars, Site!="ic")

## extract the residual correlations for each model - draws
for(i in sort(unique(rings_field$Site))) {
  for(j in sort(unique(rings_field$Neighborhood))){
    # print(paste0(i,j))
    mod <- get(paste0("fit", i, j))
    draws <- keep_draws(mod)
    
    if(i=="bm"&j==1){
      bigdraws <- draws
    } else bigdraws <- bind_rows(bigdraws, draws)
  }
}

bigdraws2 <- filter(bigdraws, Site!="ic")

## clean up workspace
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    rm(list=paste0("fit", i, j))
  }
}

## 2b. Extract correlations from time slice models----
## extract the residual correlations for each model
for(file in list.files("saved models/timeslice")) {
  print(file)
  mod <- readRDS(paste0("saved models/timeslice/", file))
  
  i <- str_sub(file, 4,5)
  j <- as.numeric(str_sub(file, 6,6))
  t <- as.numeric(str_sub(file, 7,10))
  
  vars <- cor_extract2(mod)
  
  if(i=="bm"&j==1&t==1910){
    bigvars_t <- vars
  } else bigvars_t <- bind_rows(bigvars_t, vars)
  
}

bigvars_t2 <- filter(bigvars_t, Site!="ic")


## 3a. Save output for time-invariant models----
write.csv(bigvars2, "Processed Data/synchrony_multivar.csv", row.names = F)
write.csv(bigdraws2, "Processed Data/synchrony_multivar_draws.csv", row.names = F)
## 3b. Save output for time slice models----
write.csv(bigvars_t2, "Processed Data/synchrony_decadal.csv", row.names = F)



## 2b. DON'T EVALUATE: SLOW - keep draws for decadal models

# ## function to keep all draws of residual correlations
# keep_draws2 <- function(mod) {
#   vars <- mod %>% 
#     gather_draws(`rescor.*`, regex=T) %>% 
#     # filter(str_detect(.variable, "pptmm")==FALSE) %>% 
#     mutate(
#       tree1=str_sub(str_split_fixed(.variable, "__", 3)[,2], start=2),
#       spp1=ifelse(str_detect(tree1, "ac")==TRUE, "ac",
#                   ifelse(str_detect(tree1, "pj")==TRUE, "pj",
#                          ifelse(str_detect(tree1, "pl")==TRUE, "pl",
#                                 ifelse(str_detect(tree1, "pp")==TRUE, "pp", "am")))),
#       tree2=str_sub(str_split_fixed(.variable, "__", 3)[,3], start=2),
#       spp2=ifelse(str_detect(tree2, "ac")==TRUE, "ac",
#                   ifelse(str_detect(tree2, "pj")==TRUE, "pj",
#                          ifelse(str_detect(tree2, "pl")==TRUE, "pl",
#                                 ifelse(str_detect(tree2, "pp")==TRUE, "pp", "am")))),
#       pair=str_c(spp1, spp2, sep="-"),
#       comp=ifelse(spp1==spp2, "intra", "inter"),
#       Site=i,
#       Neighborhood=j,
#       time=t)
# }
# 
# for(file in list.files("saved models/timeslice")) {
#   print(file)
#   mod <- readRDS(paste0("saved models/timeslice/", file))
#   
#   i <- str_sub(file, 4,5)
#   j <- as.numeric(str_sub(file, 6,6))
#   t <- as.numeric(str_sub(file, 7,10))
#   
#   draws <- keep_draws2(mod)
#   
#   if(i=="bm"&j==1&t==1910){
#     big_draws <- draws
#   } else big_draws <- bind_rows(big_draws, draw
#                                 
# }




