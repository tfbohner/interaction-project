## Script to estimate average site productivity (growth from rings)

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

sub <- filter(rings_field, Species!="AM") %>% 
  ungroup() %>% 
  mutate(totBA=BA_con+BA_het,
         # Species = ifelse(Species=="PP", "PJ", as.character(Species)),
         Species=droplevels(Species),
         Region=droplevels(Region))

summary <- sub %>% 
  group_by(Site, Species, tree.uniqueID) %>% 
  summarize_at(vars(totBA, DBH), mean)

hist(sub$BAI)

## Model site-species BAI----
site_diff <- brm(BAI~DBH + totBA + Region*hilo + (1|Year), family=lognormal, data=sub, 
                 cores=4, control=list(adapt_delta=0.90))

site_diff2 <- brm(BAI~DBH + totBA + Species + Region*hilo + (1|Year), family=lognormal, data=sub, 
                 cores=4, control=list(adapt_delta=0.90))


site_diff3 <- brm(BAI~DBH + totBA*Species + Region*hilo + (1|Year), family=lognormal, data=sub,
                  cores=4, control=list(adapt_delta=0.90))


site_diff4 <- brm(BAI~DBH + totBA + Species + Region + totBA:Species + Species:Region +Region:hilo + (1|Year), family=lognormal, data=sub,
                  cores=4, control=list(adapt_delta=0.99)) ## incomplete factor crosses so problematic estimation


site_diff5 <- brm(BAI~DBH + totBA + (totBA|Species) + (1|Species/Region/hilo) + (1|Year), family=lognormal, data=sub, 
                  cores=4, control=list(adapt_delta=0.99, max_treedepth=15))



l1 <- loo(site_diff)
l2 <- loo(site_diff2)
l3 <- loo(site_diff3)

l4 <- loo(site_diff4)
l5 <- loo(site_diff5)


loo_model_weights(list(l1, l2, l3))

# moving forward with site_diff3, can't get group terms to converge for 5 and 
#### 4 estimates pairings that have no data...not sure if this is actually a problem though

## Predict site-species BAI----
pred_spp <- data.frame(Species=c("AC", "PJ", "PL", "PP"), Region="SEKI", hilo="high", DBH=50, totBA=1) %>% 
  add_fitted_draws(site_diff3, re_formula = NA, allow_new_levels=TRUE) %>% 
  ungroup() %>% 
  dplyr::select(c(Species, .draw, .value)) %>% 
  mutate(.draw=rep(seq(1,4000, by=1), times=4))

pred_spatial <- data.frame(Species="AC", Region=c("Mammoth", "San Jac", "SEKI", "SENF"), DBH=50, totBA=1) %>% 
  expand_grid(hilo=c("high", "low")) %>% 
  mutate(Site=c("sl", "ic", "bm", "sp", "pr", "cm", "pp", "lc")) %>% 
  add_fitted_draws(site_diff3, re_formula = NA, allow_new_levels=TRUE)

## check out predicted data 
ggplot(pred_spp, aes(x = Species, y = .value, fill=Species)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE) +
  # geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("Predicted BAI") 
  
ggplot(pred_spatial, aes(x = Region, y = .value, fill=hilo)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.7) +
  # geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("Predicted BAI") 
  
## Summarize draws----
## Species differences

spp_summary <- expand_grid(spp1=c("AC", "PJ", "PL", "PP"), spp2=c("AC", "PJ", "PL", "PP")) %>% 
  left_join(pred_spp, by=c("spp1"="Species")) %>% 
  left_join(pred_spp, by=c("spp2"="Species", ".draw")) %>% 
  mutate(pair=tolower(paste(spp1, spp2, sep="-")),
         spp_diff=.value.x-.value.y) %>% 
  group_by(pair) %>% 
  summarize(spp_diff=median(spp_diff),
            spp_diff_upper=quantile(spp_diff, probs = 0.975),
            spp_diff_lower=quantile(spp_diff, probs = 0.025))
    
site_summary <- pred_spatial %>% 
  group_by(Region, hilo, Site, DBH, totBA) %>% 
  summarize(BAI_pred=median(.value),
            BAI_upper=quantile(.value, probs = 0.975),
            BAI_lower=quantile(.value, probs = 0.025))


## Export low density site-species growth averages----
write.csv(site_summary, "Processed Data/pred_site_grow.csv", row.names = F)
write.csv(spp_summary, "Processed Data/pred_spp_diff.csv", row.names = F)

pred_dat <- sub %>% 
  group_by(Species, Region, hilo, Site) %>% 
  data_grid(BAI, DBH=50,
            totBA=1) %>% 
  add_fitted_draws(site_diff3, n=100, re_formula = NA, allow_new_levels=TRUE) ## re_formula=NA includes no group-level effects (Year)

## check out predicted data 
ggplot(pred_dat, aes(x = Species, y = .value, fill=Region, alpha=hilo)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_alpha_ordinal(range=c(0.3,1)) +
  # geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("Predicted BAI") 


pred_summary <- pred_dat %>% 
  group_by(Species, Region, hilo, Site, DBH, totBA) %>% 
  summarize(BAI_pred=median(.value),
            BAI_upper=quantile(.value, probs = 0.975),
            BAI_lower=quantile(.value, probs = 0.025))
