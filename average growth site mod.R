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

sub <- filter(rings_field, Site!="ic", Site!="sl", Year>1990) %>% 
  ungroup() %>% 
  mutate(totBA=BA_con+BA_het,
         Species = ifelse(Species=="PP", "PJ", as.character(Species)),
         Region=droplevels(Region))

summary <- sub %>% 
  group_by(Site, Species, tree.uniqueID) %>% 
  summarize_at(vars(totBA, DBH), mean)

hist(sub$BAI)


site_diff <- brm(BAI~DBH + totBA + Region*hilo + (1|Year), family=lognormal, data=sub, 
                 cores=4, control=list(adapt_delta=0.99))

site_diff2 <- brm(BAI~DBH + totBA + Region*hilo, family=lognormal, data=sub, 
                 cores=4, control=list(adapt_delta=0.99))


site_diff3 <- brm(BAI~DBH + totBA + Species*Region +Region:hilo + (1|Year), family=lognormal, data=sub, 
                  cores=4, control=list(adapt_delta=0.99))


site_diff4 <- brm(BAI~DBH + totBA + Species + Region + totBA:Species + Species*Region +Region:hilo + (1|Year), family=lognormal, data=sub, 
                  cores=4, control=list(adapt_delta=0.99))


site_diff5 <- brm(BAI~DBH + totBA  + (1|Species/Region) + (1|Year), family=lognormal, data=sub, 
                  cores=4, control=list(adapt_delta=0.99))


site_diff6 <- brm(BAI| trunc(lb = 2, ub = 15000) ~ DBH + totBA + Species + Region + totBA:Species + Species*Region +Region:hilo + (1|Year), 
                  family=lognormal, data=sub, 
                  cores=4, control=list(adapt_delta=0.99))

pp_check(site_diff6)

l3 <- loo(site_diff3)
l4 <- loo(site_diff4)
l5 <- loo(site_diff5)
l6 <- loo(site_diff6)

loo_model_weights(list(l3, l4, l5, l6))

## site_diff4 is best model...here's the code to get the simulated conditional posteriors

data %>%
  group_by(Region) %>% 
  data_grid(Y, spei12 = seq_range(spei12, n = 10),
            DBH = c(min(data$DBH), 0, max(data$DBH)),
            comp = c(min(data$comp), 0, max(data$comp)),
            tree.uniqueID=1) %>%
  add_fitted_draws(model, n = 50, re_formula = NULL, allow_new_levels=TRUE) %>% 
  ungroup() %>% 
  mutate(DBH=as.factor(round(DBH*data$sd_dbh+data$mean_dbh, digits=2)),
         comp=as.factor(round(comp*data$sd_comp+data$mean_comp, digits=2)))

pred_dat <- sub %>% 
  group_by(Species, Region, hilo) %>% 
  data_grid(BAI, DBH=50,
            totBA=1) %>% 
  add_fitted_draws(site_diff4, n=100, re_formula = NA, allow_new_levels=TRUE)

ggplot(pred_dat, aes(x = Species, y = .value, fill=Region, alpha=hilo)) +
  geom_violin(aes(y = .value)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_alpha_ordinal(range=c(0.3,1)) +
  # geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("Predicted BAI") 




