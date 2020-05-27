library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(modelr)
library(DescTools)

## Data----
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(comp_BA=BA_con+BA_het) %>%
  arrange(Site, Neighborhood, Species, ID)

region <- rings_field %>% 
  ungroup() %>% 
  dplyr::select(c(Site, Region, hilo)) %>% 
  group_by(Site) %>% 
  summarize_all(first) %>% 
  mutate(Region=factor(Region, levels=c("San Jac", "SENF", "SEKI", "Mammoth")))

ids <- ungroup(rings_field) %>% 
  dplyr::select(c(Site, Neighborhood, tree.uniqueID)) %>% 
  group_by(Site, Neighborhood, tree.uniqueID) %>%  summarize_all(first)

tree_data <- rings_field %>% 
  group_by(tree.uniqueID) %>% 
  summarise(DBH=first(DBH),
            comp_BA=first(comp_BA))

precip <- read.csv("Processed Data/precip_temp_spei.csv")

decade <- seq(1910, 2020, by=10)

precip_summary_raw <- expand_grid(precip, decade) %>% 
  filter(Year-decade>-29, Year<=decade) %>% 
  group_by(Site, siteno, Neighborhood, hilo, decade) %>% 
  summarise_at(vars(total_ppt_mm, mean_temp_C, spei12), list(mean=function(x) mean(x, na.rm=T), sd=function(x) sd(x, na.rm=T), cv=function(x) sd(x, na.rm=T)/mean(x, na.rm=T))) %>% 
  left_join(region) %>% 
  ungroup()

precip_summary <- precip_summary_raw %>% 
  mutate_at(vars(total_ppt_mm_mean, mean_temp_C_mean, spei12_mean, 
                 total_ppt_mm_sd, mean_temp_C_sd, spei12_sd,
                 total_ppt_mm_cv, mean_temp_C_cv, spei12_cv), .funs = function(x) scale(x))

spp_pairs <- expand_grid(Species.x=c("ac", "pj", "pl", "pp"), Species.y=c("ac", "pj", "pl", "pp")) %>% 
  mutate(pair_first=str_c(Species.x, Species.y, sep="-"),
         pair=c("ac-ac", "ac-pj", "ac-pl", "ac-pj", "ac-pj", "pj-pj", "pj-pl", "pj-pp", 
                "ac-pl", "pj-pl", "pl-pl", "pl-pp", "ac-pp", "pj-pp", "pl-pp", "pp-pp")) %>% 
  dplyr::select(-c(pair_first))

pearson <- read.csv("Processed Data/synchrony_dat.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate_at(vars(Species.x, Species.y), tolower) %>% 
  left_join(spp_pairs) %>% 
  mutate(spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5),
         comp=ifelse(Species.x==Species.y, "intra", "inter"))


pearson_t <- read.csv("Processed Data/synchrony_decadal_pearson.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate_at(vars(Species.x, Species.y), tolower) %>% 
  left_join(spp_pairs) %>% 
  mutate(spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5),
         comp=ifelse(Species.x==Species.y, "intra", "inter"),
         pearson_r=round(pearson_r, 4),
         pearson_r=ifelse(pearson_r>=1, 0.9999, ifelse(pearson_r<=-1,-0.9999, pearson_r)))

summary_t <- pearson_t %>% 
  dplyr::select(c(pearson_r, decade, Site.x, Neighborhood.x, pair, comp)) %>% 
  group_by(Site.x, Neighborhood.x, pair, comp, decade) %>% #
  summarise(samp.depth=length(pearson_r),
            pearson_sd=sd(pearson_r, na.rm=T),
            pearson_r=mean(pearson_r, na.rm=T)) %>% 
  left_join(precip_summary, by=c("Site.x"="Site", "Neighborhood.x"="Neighborhood", "decade")) #

## join all data ----
alldata <- pearson %>% 
  left_join(tree_data, by=c("tree2"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x),
         Region=factor(Region, levels=c("San Jac", "SENF", "SEKI", "Mammoth"))) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale)

alldata_t <- pearson_t %>% 
  left_join(tree_data, by=c("tree2"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x),
         Region=factor(Region, levels=c("San Jac", "SENF", "SEKI", "Mammoth"))) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale) %>% 
  left_join(dplyr::select(precip_summary, -c(Region, hilo)), by=c('Site.x'='Site', "Neighborhood.x"="Neighborhood", 'decade'))


### Part 1: time invariant models----
mod0 <- brm(pearson_r~ 1 + (1|Region:hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))
mod0b <- brm(pearson_r~ 1 + (1|Region/hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))

loo_model_weights(list(loo(mod0), loo(mod0b)), method="pseudobma")
l0 <- loo(mod0)

mod1 <- update(mod0, formula= ~. + dist, newdata=pearson, cores=4, control = list(adapt_delta=0.99))
l1 <- loo(mod1)

loo_model_weights(list(l0, l1))

mod2<- brm(pearson_r~dist + (1|pair) + (1|Region:hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))
l2 <- loo(mod2)

mod3<- brm(pearson_r~dist + (dist|pair) + (1|Region:hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))
mod4<- brm(pearson_r~dist + (1|pair) + (dist|Region:hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))
mod5<- brm(pearson_r~dist + (dist|pair) + (dist|Region:hilo), data=pearson, cores=4, control = list(adapt_delta=0.99))
mod6 <- update(mod4, formula= ~. + sizeratio, newdata=alldata, cores=4, control = list(adapt_delta=0.99))

mod7<- brm(pearson_r~dist + sizeratio + (1|pair) + (dist+sizeratio|Region:hilo), data=alldata, cores=4, control = list(adapt_delta=0.99))
mod8 <- update(mod6, formula=~. + comp_BA, newdata=alldata, cores=4, control = list(adapt_delta=0.99))

mod9<- brm(pearson_r~dist + sizeratio + (1|pair) + (dist|Region:hilo:Neighborhood.x), data=alldata, cores=4, control = list(adapt_delta=0.99))

mod10<- brm(pearson_r~dist + sizeratio + (1|pair/Region:hilo:Neighborhood.x) + (dist|Region:hilo:Neighborhood.x), data=alldata, cores=4, control = list(adapt_delta=0.99))

mod11<- brm(pearson_r~dist + sizeratio + (1|pair/Region:hilo) + (dist|Region:hilo) + (1|Region:hilo:Neighborhood.x), data=alldata, cores=4, control = list(adapt_delta=0.99))

mod12<- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~dist + sizeratio + (1|pair/Region:hilo) + (1|Region:hilo/Neighborhood.x), data=alldata, cores=4, control = list(adapt_delta=0.99))

for(i in 12) {
  mod <- paste0("mod",i)
  l <- loo(get(mod))
  assign(paste0("l", i), l)
}

loo_model_weights(list(l11, l12))

saveRDS(mod11, "saved models/pearson models/mod11.rds")
saveRDS(mod12, "saved models/pearson models/mod12.rds")


### Part 2: through time----
## spline models----
mod_s_t <- brm(bf(pearson_r ~ s(decade)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
conditional_effects(mod_s_t)
plot(conditional_effects(mod_s_t), points=TRUE)

saveRDS(mod_s_t, "saved models/spline.rds")


mod_s_t3 <- brm(bf(pearson_r ~ pair + s(decade, by=pair)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
mod_s_t4 <- brm(bf(pearson_r ~ Site.x + s(decade, by=Site.x)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
conditional_effects(mod_s_t4)

saveRDS(mod_s_t3, "saved models/spline_spp.rds")
saveRDS(mod_s_t4, "saved models/spline_site.rds")


## climate models----
# clim0 <- brm(pearson_r~total_ppt_mm_mean + total_ppt_mm_sd + mean_temp_C_mean + mean_temp_C_sd, data=summary_t, cores=4, control = list(adapt_delta=0.99))
# ## high predictor correlations
# vcov(clim0, correlation = T) %>% 
#   round(digits = 3)

clim1 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1b <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (1|Site.x/pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1c <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (1|Site.x + pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1d <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1e <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x) + (1|pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))

loo_model_weights(list(loo(clim1), loo(clim1e)))

clim2 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_cv + (total_ppt_mm_cv|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim3 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~mean_temp_C_mean + (mean_temp_C_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim4 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~mean_temp_C_cv + (mean_temp_C_cv|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim5 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~spei12_mean + (spei12_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim6 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~spei12_cv + (spei12_cv|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))

for(i in 1:6) {
  mod <- paste0("clim",i)
  l <- loo(get(mod))
  assign(paste0("l", i), l)
}

loo_model_weights(list(l1, l2, l3, l4, l5, l6))

pp_check(clim2)

saveRDS(clim2, "saved models/pearson models/clim2.rds")

