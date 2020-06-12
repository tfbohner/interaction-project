library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(modelr)

## Data----
treeclim <- read.csv("Processed Data/treeclim_spp.csv")

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

precip_order <- data.frame(Site=c("bm", "cm", "lc", "pp", "pr", "sl", "sp"), 
                           precip_order=as.factor(c(3, 6, 4, 5, 7, 1, 2)))
compos <- rings_field %>% 
  group_by(tree.uniqueID, Region, hilo, Site) %>% 
  filter(Site!="ic") %>% 
  summarise(Species=tolower(first(Species))) %>% 
  left_join(precip_order)


site_clim <- read.csv("Processed Data/precip_temp_spei.csv")%>% 
  filter(Site!="ic", Year>1895) %>% 
  dplyr::select(c(Year, Site, Region, hilo, Neighborhood, total_ppt_mm, mean_temp_C, spei12)) %>% #, spei12
  group_by(Site, Region, hilo, Neighborhood, Year) %>% #x
  summarize_all(mean, na.rm=T) %>% 
  group_by(Site, Region, hilo) %>%
  dplyr::select(-c(Neighborhood, Year)) %>%
  summarize_all(list(mean, sd)) %>%
  ungroup() %>% 
  mutate(ppt_cv=total_ppt_mm_fn2/total_ppt_mm_fn1,
         temp_cv=mean_temp_C_fn2/mean_temp_C_fn1)

climdat <- as.matrix(dplyr::select(site_clim, -c(Site, Region, hilo))) 
climdat <- as.matrix(dplyr::select(site_clim, c(total_ppt_mm_fn1, mean_temp_C_fn1))) 
rownames(climdat) <- site_clim$Site
pca <- prcomp(climdat, scale=T)

site_clim <- site_clim %>%  
  mutate(pca1=pca$x[,1],
         pca2=pca$x[,2],
         Region=str_replace(Region, pattern = " ", replacement = "."),
         site=paste(Region, hilo, sep="_"))

precip <- read.csv("Processed Data/precip_temp_spei.csv")

decade <- seq(1910, 2020, by=10)

precip_summary_raw <- expand_grid(precip, decade) %>% 
  filter(Year-decade>-29, Year<=decade) %>% 
  filter(Site!="ic") %>% 
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

site_corr <- pearson %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  group_by(Site.x) %>% 
  summarise(site_r=mean(pearson_r)) %>% 
  left_join(site_clim, by=c("Site.x"="Site")) 


pearson_t <- read.csv("Processed Data/synchrony_decadal_pearson.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate_at(vars(Species.x, Species.y), tolower) %>% 
  left_join(spp_pairs) %>% 
  mutate(spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5),
         comp=ifelse(Species.x==Species.y, "intra", "inter"))
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
         realsizeratio=sizeratio,
         realdist=dist,
         Region=factor(Region, levels=c("San Jac", "SENF", "SEKI", "Mammoth"))) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale)

alldata_t <- pearson_t %>% 
  left_join(tree_data, by=c("tree2"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x),
         Region=factor(Region, levels=c("San Jac", "SENF", "SEKI", "Mammoth"))) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale) %>% 
  left_join(dplyr::select(precip_summary, -c(Region, hilo)), by=c('Site.x'='Site', 'decade'='decade'))

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

# clim1 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1b <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (1|Site.x/pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1c <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (1|Site.x + pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1d <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim1e <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x) + (1|pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))

loo_model_weights(list(loo(clim1), loo(clim1e)))

clim1 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_mean + (total_ppt_mm_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim2 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~total_ppt_mm_cv + (total_ppt_mm_cv|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim3 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~mean_temp_C_mean + (mean_temp_C_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim4 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~mean_temp_C_cv + (mean_temp_C_cv|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim5 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~spei12_mean + (spei12_mean|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim6 <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~spei12_sd + (spei12_sd|Site.x) + (1|Site.x:pair), data=summary_t, cores=4, control = list(adapt_delta=0.99))
clim6b <- brm(pearson_r|trunc(lb=-1.001, ub=1.001)~spei12_sd + (spei12_sd|Site.x:pair) , data=summary_t, cores=4, control = list(adapt_delta=0.99))

for(i in 1:6) {
  mod <- get(paste0("clim",i))
  l <- loo(mod)
  assign(paste0("l", i), l)
  
  mod <- add_criterion(mod, "waic")
  assign(paste0("clim", i), mod)
  
  # saveRDS(mod, paste0("saved models/pearson models/clim", i, ".rds"))
}

loo_model_weights(list(l1, l2, l3, l4, l5, l6))
loo_model_weights(list(l6, l6b))

w <- loo_compare(clim1, clim2, clim3, clim4, clim5, clim6, criterion='waic')
cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] * 2) %>% 
  round(digits=2)

model_weights(clim1, clim2, clim3, clim4, clim5, clim6,
              weights = "waic") %>% 
  round(digits=2)

model_weights(clim6, clim6b,
              weights = "waic") %>% 
  round(digits=2)

pp_check(clim2)

for(i in 1:6) {
  mod <- readRDS(paste0("saved models/pearson models/clim", i, ".rds"))
  assign(paste0("clim",i), mod)
  
}
