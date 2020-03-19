## This script will run the multivariate brms model
library(tidyverse)
library(brms)
library(tidybayes)

## Data ---
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) ## counter id if we want it rather than unique ID

## Example neighborhood----
neigh1 <- rings_field %>% 
  filter(siteno==11) %>% 
  filter(Neighborhood==1)

grow_mat <- neigh1 %>% 
  pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
  na.omit() %>% ## can include missing data if we want
  mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)

names(grow_mat) <- c("year", "spei12", "mean_temp_C","total_ppt_mm","a1", "a2", "a3", "a4", "b5", "b6", "b7", "b8", "c9", "c10", "c11", "c12")


## no predictors, corr mat is close to pearson
fit <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0, sigma~0), data=grow_mat, cores=4)

# mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) is the section that I can't figure out how to turn into a variable
# that would allow for different names for each neighborhood and speed up the looping process or allow me to create a custom function

## with SPEI and annual temp
fit2 <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0 + spei12 + mean_temp_C), data=grow_mat, cores=4)

## with missing data 
bf1 <- bf(a1|mi()~0)
bf2 <- bf(a2|mi()~0)
bf3 <- bf(a3|mi()~0) #... and so on 

# fit3 <- brm(mvbf(bf1, bf2, bf3...), data=grow_mat, cores=4


## Extract parameters from no predictor model
vars <- fit %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

## group by interacting species
spp_groups <- vars %>% 
  group_by(pair) %>% 
  summarize_all(median)

ggplot(spp_groups,aes(x=.value, y=pair)) +
  geom_point()

## Extract parameters from clim predictor model
vars2 <- fit2 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )


spp_groups2 <- vars2 %>% 
  group_by(pair) %>% 
  summarize_all(median)

ggplot(spp_groups2,aes(x=.value, y=pair)) +
  geom_point()

#### Looping or function ----

# I'm liking the function route lately...I did it a lot in the drought analysis because you can specify 
# parameters...may be the solution for coding the response vars...this would need some changing for it to actually work

form_maker <- function(col, predictors) {
  form <- bf(as.formula(paste0(names(grow_mat[,col]), "~", predictors)))
}

# and/or loop
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    data_sub <- filter(rings_field, Site==i&Neighborhood==j)
    
    grow_mat <- data_sub %>% 
      pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
      na.omit() %>% 
      mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
    
    for(n in 5:ncol(grow_mat)) {
      form <- form_maker(n, "spei12 + mean_temp_C")
      assign(paste0("bf", n), form)
    }

    if(ncol(grow_mat)==12){
      mod <- brm(mvbf(bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==13){
      mod <- brm(mvbf(bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==14){
      mod <- brm(mvbf(bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==15){
      mod <- brm(mvbf(bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==16){
      mod <- brm(mvbf(bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15, bf16), data=grow_mat, cores=4)
    }}}}}
    
  

    
    
    assign(paste0("fit", i, j), mod)
    saveRDS(mod, file=paste0("saved models/with covar/", "fit", i, j, ".rds"))
  }
}


## Check lengths of columns for if statement, 12-16 (so 8-12 response vars)
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    data_sub <- filter(rings_field, Site==i&Neighborhood==j)
    
    grow_mat <- data_sub %>% 
      pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
      na.omit() %>% 
      mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
    
    print(ncol(grow_mat))
  }
}

