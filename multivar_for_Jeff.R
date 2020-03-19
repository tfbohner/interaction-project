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

model_fun <- function(data, resp, params) {
  fit2 <- brm(bf(mvbind(resp) ~ 0 + params ), data=data, cores=4)
}

# and/or loop
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    data_sub <- filter(rings_field, Site==i&Neighborhood==j)
    
    grow_mat <- data_sub %>% 
      pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
      na.omit() %>% 
      mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
    
    mod <- model_fun(grow_mat, resp, params)
    
    assign(paste("fit", i, j, sep="_"), mod)
  }
}
