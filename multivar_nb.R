## This script will run the multivariate brms model and calculate residual individual tree correlations
## Both the models for the entire time series and the moving window approach are in this script.
library(tidyverse)
library(brms)
library(tidybayes)

## Data ---
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)


## Example neighborhood----
neigh1 <- rings_field %>% 
  filter(siteno==11) %>% 
  filter(Neighborhood==1) 

ggplot(neigh1, aes(Year, nexp_growth, color=tree.uniqueID)) + geom_line() + facet_wrap(~tree.uniqueID) + theme(legend.position = 'none')
ggplot(neigh1, aes(Year, raw_growth, color=tree.uniqueID)) + geom_line() + facet_wrap(~tree.uniqueID) + theme(legend.position = 'none')
ggplot(neigh1, aes(Year, BAI, color=tree.uniqueID)) + geom_line() + facet_wrap(~tree.uniqueID) + theme(legend.position = 'none')
ggplot(neigh1, aes(Year, spline_growth, color=tree.uniqueID)) + geom_line() + facet_wrap(~tree.uniqueID) + theme(legend.position = 'none')


count_na <- function(x) sum(is.na(x))

grow_mat <- neigh1 %>% 
  pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
  # na.omit() %>% ## can include missing data if we want
  mutate(nacount=apply(., 1, count_na)) %>% 
  filter(nacount<4) %>% 
  select(-nacount) #%>% 
  #mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)

names(grow_mat) <- c("year", "spei12", "mean_temp_C","total_ppt_mm","a1", "a2", "a3", "a4", "b5", "b6", "b7", "b8", "c9", "c10", "c11", "c12")
apply(grow_mat, 2, count_na)

## no predictors, corr mat is close to pearson
fit.a <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0, sigma~0), data=grow_mat, cores=4)
fit.b <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 1), data=grow_mat, cores=4)
fit.c <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 1), data=grow_mat, cores=4)
fit.d <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 1), data=grow_mat, cores=4)

fit.e <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 1 +(1|year)), data=grow_mat, cores=4)

fit4 <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0), 
            family=student, 
            prior = prior(gamma(2, .1), class = nu),
            data=grow_mat, cores=4)

# mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) is the section that I can't figure out how to turn into a variable
# that would allow for different names for each neighborhood and speed up the looping process or allow me to create a custom function

## with SPEI and annual temp
fit2 <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0 + spei12 + mean_temp_C), data=grow_mat, cores=4)

## with missing data 
bf1 <- bf(a1|mi()~0)
bf2 <- bf(a2~0)
bf3 <- bf(a3~0)
bf4 <- bf(a4~0)
bf5 <- bf(b5~0)
bf6 <- bf(b6|mi()~0)
bf7 <- bf(b7|mi()~0)
bf8 <- bf(b8~0)
bf9 <- bf(c9~0)
bf10 <- bf(c10~0)
bf11 <- bf(c11~0)
bf12 <- bf(c12~0)


fit3 <- brm(bf1 + bf2 + bf3 + bf4 + bf5 + bf6 + bf7 + bf8 + bf9 + bf10 + bf11 + bf12 + set_rescor(TRUE),
            family=student, 
            prior = prior(gamma(2, .1), class = nu),
            data=grow_mat, cores=4)



## Extract parameters from no predictor model
vars <- fit.e %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

vars4 <- fit4 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

plot(vars$.value, vars4$.value)
abline(0,1)

vars3 <- fit3 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

plot(vars$.value, vars4$.value)
abline(0,1)

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

#### Looping through neighborhoods ----

# I'm liking the function route lately...I did it a lot in the drought analysis because you can specify 
# parameters...may be the solution for coding the response vars...this would need some changing for it to actually work

form_maker <- function(col, predictors) {
  form <- bf(as.formula(paste0(names(grow_mat[,col]), "~", predictors)))
}

# loop through neighborhoods
for(i in unique(rings_field$Site)) {
  for(j in unique(rings_field$Neighborhood)){
    data_sub <- filter(rings_field, Site==i&Neighborhood==j)
    
    grow_mat <- data_sub %>% 
      pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
      na.omit() %>% 
      mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
    
    for(n in 5:ncol(grow_mat)) {
      form <- form_maker(n, "0") #spei12 + mean_temp_C
      assign(paste0("bf", n), form)
    }

    if(ncol(grow_mat)==12){
      mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==13){
      mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==14){
      mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==15){
      mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15), data=grow_mat, cores=4)
    } else{if(ncol(grow_mat)==16){
      mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15, bf16), data=grow_mat, cores=4)
    }}}}}
    
    assign(paste0("fit", i, j), mod)
    saveRDS(mod, file=paste0("saved models/", "fit", i, j, ".rds")) #with covar/
  }
}


#### Looping through timme chunks ----
timeseq <- seq(1910, 2010, by=10)

for(i in unique(rings_field$Site)[4:8]) {
  for(j in unique(rings_field$Neighborhood)){
    for(t in timeseq) {
      data_sub <- filter(rings_field, Site==i&Neighborhood==j) %>% 
        filter(Year>t-11, Year<t+11)
      
      if(t==2010){
        grow_mat <- data_sub %>% 
          pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
          na.omit() %>%
          # select_if(~ !any(is.na(.))) %>%
          mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
      } else{
        grow_mat <- data_sub %>% 
          pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
          # na.omit() %>%
          select_if(~ !any(is.na(.))) %>%
          mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)
      }
      
      if(nrow(grow_mat)>1){
        
        for(n in 5:ncol(grow_mat)) {
          form <- form_maker(n, "0") #spei12 + mean_temp_C
          assign(paste0("bf", n), form)
        }
        
        if(ncol(grow_mat)==9){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==10){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==11){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==12){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==13){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==14){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==15){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        } else{if(ncol(grow_mat)==16){
          mod <- brm(mvbf(bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12, bf13, bf14, bf15, bf16), family=student, data=grow_mat, cores=4, control=list(adapt_delta=0.99))
        }}}}}}}}
        
        
        # assign(paste0("fit", i, j, t), mod)
        saveRDS(mod, file=paste0("saved models/timeslice/", "fit", i, j, t, ".rds")) #with covar/
      }
    }
  }
}


vars1 <- fitbm11950 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

vars2 <- fitbm12010 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

plot(vars1$.value, vars2$.value)
abline(0,1)
median(vars1$.value)
median(vars2$.value)

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

