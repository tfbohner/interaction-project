# Interaction project
## 2. Analyze
#### Purpose: Model intra-interpsecific correlation comparison across species and sites
#### Author: Teresa Bohner
#### Date Modified: 24 April 2020
library(tidyverse)
library(brms)

## Load data----
all_cor <- read.csv("Processed Data/inta_inter_pairwise_all.csv") %>% 
  mutate(pair2 =str_c(intra_spp,inter_spp, sep="-"))
id_cor_summ <- read.csv("Processed Data/inta_inter_pairwise_summary.csv") %>% 
  mutate(pair2 =str_c(intra_spp,inter_spp, sep="-"))

tree_data <- rings_field %>% 
  group_by(tree.uniqueID) %>% 
  summarise(DBH=first(DBH))

spp_diff <- read.csv("Processed Data/pred_spp_diff.csv")
site_grow <- read.csv("Processed Data/pred_site_grow.csv") %>% 
  dplyr::select(-c(DBH, totBA))

big_dat <- left_join(all_cor, site_grow)

## Explore----
sub <- filter(id_cor_summ, intra_spp!="pp", inter_spp!="pp", Region!="Mammoth")
ggplot(sub, aes(Region, prop, fill=hilo)) +
  geom_boxplot()

summary <- id_cor_summ %>% 
  group_by(Region, inter_spp, intra_spp) %>% 
  summarise(n=length(intra_greater))

## Bernoulli model----

bern_mod0 <- brm(intra_greater~ 1 + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod1 <- brm(intra_greater~ 1 + (1|Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod2 <- brm(intra_greater~ 1 + (1|Region/hilo) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))

bern_mod3 <- brm(intra_greater~ 1 + intra_spp + (1|Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod4 <- brm(intra_greater~ 1 + intra_spp + intra_spp:inter_spp + (1|Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod5 <- brm(intra_greater~ 1 + pair2 + (1|Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod6 <- brm(intra_greater~ 1 + (1|intra_spp/pair2) + (1|Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod7 <- brm(intra_greater~ 1 + (1|Region/intra_spp/pair2) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod8 <- brm(intra_greater~ 1 + (1|intra_spp/pair2/Region) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))
bern_mod9 <- brm(intra_greater~ 1 + (1|Site/intra_spp/pair2) + (1|tree1.x), data=all_cor, family=bernoulli, cores=4, control=list(adapt_delta=0.99))

for(i in 6:7) {
  mod <- paste0("bern_mod",i)
  l <- loo(get(mod))
  assign(paste0("l", i), l)
}

loo_model_weights(list(l6, l7, l9), method='pseudobma')

saveRDS(bern_mod6, "saved models/intra_greater_mods/bern_mod6.rds")

saveRDS(bern_mod7, "saved models/intra_greater_mods/bern_mod7.rds")
saveRDS(bern_mod9, "saved models/intra_greater_mods/bern_mod9.rds")

## Visualize model results----library(viridis)
library(tidybayes)
all_cor <- read.csv("Processed Data/inta_inter_pairwise_all.csv") %>% 
  mutate(pair2 =str_c(intra_spp,inter_spp, sep="-"))

bern_mod6 <- readRDS("saved models/intra_greater_mods/bern_mod6.rds")
bern_mod7 <- readRDS("saved models/intra_greater_mods/bern_mod7.rds")

dat_post <- bern_mod6 %>% 
  spread_draws(b_Intercept, r_intra_spp[species,]) %>% 
  mutate(spp_med=inv_logit_scaled(b_Intercept + r_intra_spp)) %>% 
  dplyr::select(-c(b_Intercept, r_intra_spp)) %>% 
  median_qi()


ggplot(dat_post, aes(x=spp_med, y=species)) +
  geom_pointintervalh()


dat_post <- bern_mod6 %>% 
  spread_draws(b_Intercept, `r_intra_spp:pair2`[comp,]) %>% 
  mutate(spp_med=inv_logit_scaled(b_Intercept + `r_intra_spp:pair2`)) %>% 
  dplyr::select(-c(b_Intercept, `r_intra_spp:pair2`)) %>% 
  median_qi()


ggplot(dat_post, aes(x=spp_med, y=comp)) +
  geom_pointintervalh()

dat_post <- bern_mod6 %>% 
  spread_draws(b_Intercept, r_Region[region,]) %>% 
  mutate(reg_med=inv_logit_scaled(b_Intercept + r_Region)) %>% 
  dplyr::select(-c(b_Intercept, r_Region)) %>% 
  median_qi()


ggplot(dat_post, aes(x=reg_med, y=region)) +
  geom_pointintervalh()


dat_post <- bern_mod7 %>% 
  spread_draws(b_Intercept, `r_Region:intra_spp:pair2`[comp,]) %>% 
  mutate(spp_med=inv_logit_scaled(b_Intercept + `r_Region:intra_spp:pair2`)) %>% 
  dplyr::select(-c(b_Intercept, `r_Region:intra_spp:pair2`)) %>% 
  median_qi()


ggplot(dat_post, aes(x=spp_med, y=comp)) +
  geom_pointintervalh()


## predicted posterior----
pred_dat <- all_cor %>% 
  group_by(Region, intra_spp, pair2) %>% 
  data_grid(intra_greater) %>% 
  add_fitted_draws(bern_mod7, re_formula = ~(1|Region/intra_spp/pair2), allow_new_levels=TRUE) ## re_formula=NA includes no group-level effects (Year)

## check out predicted data 
ggplot(pred_dat, aes(x = intra_spp, y = .value, fill=pair2, alpha=Region)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_alpha_ordinal(range=c(0.3,1)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted proportion Raa > Rab") 


pred_dat <- all_cor %>% 
  group_by(Region, intra_spp, pair2) %>% 
  data_grid(intra_greater) %>% 
  add_fitted_draws(bern_mod6, re_formula = ~(1|Region) + (1|intra_spp/pair2), allow_new_levels=TRUE) ## re_formula=NA includes no group-level effects (Year)

## check out predicted data 
ggplot(pred_dat, aes(x = intra_spp, y = .value, fill=pair2, alpha=Region)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_alpha_ordinal(range=c(0.3,1)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted proportion Raa > Rab") 
