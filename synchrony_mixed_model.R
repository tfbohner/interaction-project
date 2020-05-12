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
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)

tree_data <- rings_field %>% 
  group_by(tree.uniqueID) %>% 
  summarise(DBH=first(DBH))

spp_diff <- read.csv("Processed Data/pred_spp_diff.csv") %>% 
  mutate(spp_abs_diff = abs(spp_diff))
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
# pair_mod3 <- brm(.value ~ dist*pair + sizeratio + BAI_pred, data=alldata, cores=4)
pair_mod3 <- update(pair_mod2, formula. =  ~ .  + dist:pair, newdata=alldata, cores=4)
# pair_mod4 <- brm(.value ~ dist*pair + sizeratio + Region:hilo, data=alldata, cores=4)
pair_mod4 <- update(pair_mod3, formula. =  ~ . - BAI_pred + Site, newdata=alldata, cores=4)

pair_mod4a <- update(pair_mod3, formula. =  ~ . - pair -dist:pair + comp +dist:comp, newdata=alldata, cores=4)

pair_mod5 <- brm(.value ~ dist + sizeratio + (dist|pair) + (1|Site), 
                 data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod6 <- brm(.value ~ dist + sizeratio + (dist|pair) + (1|Region/hilo), 
                 data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod7 <- brm(.value ~ dist + sizeratio + (1|pair) + (1|Region/hilo), 
                 data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

pair_mod8 <- brm(.value ~ dist*pair + sizeratio + (1|Region/hilo), 
                 data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

pair_mod9 <-  brm(.value ~ dist + sizeratio + (1|Region/hilo/pair), 
              data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod10 <-  brm(.value ~ dist + sizeratio + (1|Site/pair), 
                  data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

pair_mod11 <-  brm(.value ~ dist + sizeratio + (1|Site/spp_abs_diff), 
                   data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

pair_mod12 <-  brm(.value ~ dist + sizeratio + (1|Region/pair), 
                   data=alldata, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

l2 <- loo(pair_mod2)
l3 <- loo(pair_mod3)
l4 <- loo(pair_mod4)
l5 <- loo(pair_mod5)
l6 <- loo(pair_mod6)
l7 <- loo(pair_mod7)
l8 <- loo(pair_mod8)
l9 <- loo(pair_mod9)
l10 <- loo(pair_mod10)
l11 <- loo(pair_mod11)
l12 <- loo(pair_mod12)

loo_model_weights(list(l9, l12), method = "stack")

## model10 or 9

dat <- alldata %>% 
  filter(!is.na(sizeratio)) %>% 
  mutate(cor=.value)

write.csv(dat, "Processed Data/data_for_sync_mod.csv")

pair_mod10a <-  update(pair_mod10, formula= ~ . - sizeratio, newdata=dat, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod10b <-  update(pair_mod10, formula= ~ . - dist, newdata=dat, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod10c <-  update(pair_mod10, formula= ~ . - sizeratio - dist, newdata=dat, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))
pair_mod10d <- update(pair_mod10, formula= ~ . + sizeratio:dist, newdata=dat, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

l10a <- loo(pair_mod10a)
l10b <- loo(pair_mod10b)
l10c <- loo(pair_mod10c)
l10d <- loo(pair_mod10d)

loo_model_weights(list(l10, l10a, l10b, l10c), method = "stack")

dat <- alldata %>% 
  filter(!is.na(sizeratio)) %>% 
  mutate(cor=.value,
         dist=scale(dist),
         sizeratio=scale(sizeratio))

pair_mod10 <-  brm(cor ~ dist + sizeratio + (1|Site/pair), 
                   data=dat, cores=4, control=list(adapt_delta=0.99, max_treedepth=13))

saveRDS(pair_mod9, "saved models/sync_mods/pair_mod9.rds")
saveRDS(pair_mod10, "saved models/sync_mods/pair_mod10.rds")

dat_post <- pair_mod10 %>% 
  spread_draws(b_Intercept, `r_Site:pair`[comp,]) %>% 
  mutate(comp_med=inv_logit_scaled(b_Intercept + `r_Site:pair`)) %>% 
  dplyr::select(-c(b_Intercept, `r_Site:pair`)) %>% 
  median_qi()


ggplot(dat_post, aes(x=comp_med, y=comp)) +
  geom_pointintervalh()

pred_dat <- dat %>% 
  group_by(Site, pair) %>% 
  data_grid(cor, dist=round(quantile(dat$dist, probs=c(0.5)), digits=2), 
            sizeratio=round(quantile(dat$sizeratio, probs=c(0.5)), digits=2)) %>% 
  add_fitted_draws(pair_mod10, n=100, re_formula = NULL, allow_new_levels=TRUE) %>% 
  mutate(spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5))

## check out predicted data 
ggplot(pred_dat, aes(x = Site, y = .value, fill=Site)) +
  geom_violin(aes(y = .value)) +
  facet_grid(spp1~spp2) +
  scale_fill_viridis(discrete=TRUE) +
  geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted correlation R") 

 pred_dat <- dat %>% 
  group_by(Site, pair) %>% 
  data_grid(cor, dist=round(quantile(dat$dist, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), digits=2), 
            sizeratio=round(quantile(dat$sizeratio, probs=c(0.5)), digits=2)) %>% 
  add_fitted_draws(pair_mod10, n=100, re_formula = NULL, allow_new_levels=TRUE) 
  
ggplot(pred_dat, aes(x = dist, y = .value)) +
  stat_lineribbon(aes(y = .value), .width = c(.95), alpha=0.2) +
  stat_lineribbon(aes(y = .value), .width = c(.01)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_viridis(discrete=TRUE) +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted correlation R") 
