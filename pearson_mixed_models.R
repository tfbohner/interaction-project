library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(modelr)

## Data----
pearson <- read.csv("Processed Data/synchrony_dat.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate(pair=paste(Species.x, Species.y, sep="-"),
         pair=ifelse(pair=="PJ-PP", "PP-PJ", ifelse(pair=="PL-PP", "PP-PL", pair)),
         Species.x=ifelse(pair=="PP-PJ", "PP", as.character(Species.x)),
         Species.y=ifelse(pair=="PP-PJ", "PJ", as.character(Species.y)),
         comp=ifelse(Species.x==Species.y, "intra", "inter"))

rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(comp_BA=BA_con+BA_het) %>%
  arrange(Site, Neighborhood, Species, ID)

tree_data <- rings_field %>% 
  group_by(tree.uniqueID) %>% 
  summarise(DBH=first(DBH),
            comp_BA=first(comp_BA))

spp_diff <- read.csv("Processed Data/pred_spp_diff.csv") %>% 
  mutate(spp_abs_diff = abs(spp_diff))
site_grow <- read.csv("Processed Data/pred_site_grow.csv") %>% 
  dplyr::select(-c(DBH, totBA))

## join all data ----
alldata <- pearson %>% 
  left_join(tree_data, by=c("tree2"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x)) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale)


## Visualize----

ggplot(pearson, aes(dist, pearson_r)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x, 2), se=FALSE) +
  facet_grid(Site.x~Neighborhood.x)

ggplot(pearson, aes(pearson_r, fill=comp)) +
  geom_density(aes(y = ..scaled..), alpha=0.5) +
  facet_grid(Species.x~Species.y) +
  theme_test() +
  geom_vline(xintercept=0, linetype="dashed")

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

mod12<- brm(pearson_r~dist + sizeratio + (1|pair/Region:hilo) + (1|Region:hilo/Neighborhood.x), data=alldata, cores=4, control = list(adapt_delta=0.99))

conditional_effects(mod12)

for(i in 12) {
  mod <- paste0("mod",i)
  l <- loo(get(mod))
  assign(paste0("l", i), l)
}

loo_model_weights(list(l11, l12))

saveRDS(mod11, "saved models/pearson models/mod11.rds")
saveRDS(mod12, "saved models/pearson models/mod12.rds")

## Manip----
get_variables(mod11)

## distance versus sizeratio?
mod11 %>% 
  gather_draws(b_dist, b_sizeratio) %>% 
  median_qi() %>% 
  ggplot(aes(y=.variable, x=.value)) +
  geom_pointintervalh() +
  geom_vline(xintercept=0, linetype="dashed")


alldata %>%
  data_grid(pearson_r, dist = rnorm(5, mean(alldata$dist), sd(alldata$dist)), 
            sizeratio = 1) %>%
  add_fitted_draws(mod11, allow_new_levels=TRUE) %>%
  ggplot(aes(x = dist, y = .value)) +
  stat_lineribbon(aes(y = .value), .width = c(.95), alpha=0.2) +
  stat_lineribbon(aes(y = .value), .width = c(.01)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_viridis(discrete=TRUE) 


mod11 %>% 
  spread_draws(b_Intercept, r_pair[pair,]) %>% 
  mutate(spp_mean=b_Intercept + r_pair) %>% 
  median_qi() %>% 
  ggplot(aes(y=pair, x=spp_mean, xmin=spp_mean.lower, xmax=spp_mean.upper)) +
  geom_pointintervalh() +
  geom_vline(xintercept=0, linetype="dashed")


mod11 %>% 
  spread_draws(r_pair[pair,]) %>% 
  compare_levels(r_pair, by = pair) %>%
  median_qi() %>%
  mutate(pair1=str_sub(pair, 1,5),
         pair1_spp1=str_sub(pair, 1,2),
         pair1_spp2=str_sub(pair, 4,5),
         pair1_comp=ifelse(pair1_spp1==pair1_spp2, "intra", "inter"),
         pair2=str_sub(pair, 9,13),
         pair2_spp1=str_sub(pair, 9,10),
         pair2_spp2=str_sub(pair, 12,13),
         pair2_comp=ifelse(pair2_spp1==pair2_spp2, "intra", "inter"),
         valid=ifelse(pair1_comp=="intra", str_detect(pair2, pair1_spp1), str_detect(pair1, pair2_spp1))) %>% 
  filter(pair1_comp!=pair2_comp, valid==TRUE) %>% 
  mutate(new_pair=ifelse(pair1_comp=="intra", str_c(pair1, pair2, sep=" - "), str_c(pair2, pair1, sep=" - ")),
         r_pair=ifelse(pair1_comp=="intra", r_pair, r_pair*-1),
         .upper=ifelse(pair1_comp=="intra", .upper, .upper*-1),
         .lower=ifelse(pair1_comp=="intra", .lower, .lower*-1)) %>% 
  ggplot(aes(y = new_pair, x = r_pair)) +
  geom_pointintervalh() +
  geom_vline(xintercept=0, linetype="dashed")

mod11 %>% 
  spread_draws(b_Intercept, r_pair[pair,], `r_pair:Region:hilo`[site,]) %>% 
  mutate(spp_mean=b_Intercept + r_pair + `r_pair:Region:hilo`,
         valid=str_detect(site, pair)) %>% 
  filter(valid==1) %>% 
  median_qi() %>% 
  mutate(newsite=str_sub(site, start=7)) %>% 
  ggplot(aes(y=pair, x=spp_mean, xmin=spp_mean.lower, xmax=spp_mean.upper, color=newsite)) +
  geom_pointintervalh(position = position_dodgev(height=0.3)) +
  geom_vline(xintercept=0, linetype="dashed") +
  theme(legend.position = 'none')

mod11 %>% 
  spread_draws(r_pair[pair,], `r_pair:Region:hilo`[site,]) %>% 
  mutate(spp_mean=r_pair + `r_pair:Region:hilo`,
         valid=str_detect(site, pair)) %>% 
  filter(valid==1) %>% 
  mutate(newsite=str_sub(site, start=7)) %>% 
  group_by(newsite) %>% 
  compare_levels(spp_mean, by = pair) %>%
  median_qi() %>%
  mutate(pair1=str_sub(pair, 1,5),
         pair1_spp1=str_sub(pair, 1,2),
         pair1_spp2=str_sub(pair, 4,5),
         pair1_comp=ifelse(pair1_spp1==pair1_spp2, "intra", "inter"),
         pair2=str_sub(pair, 9,13),
         pair2_spp1=str_sub(pair, 9,10),
         pair2_spp2=str_sub(pair, 12,13),
         pair2_comp=ifelse(pair2_spp1==pair2_spp2, "intra", "inter"),
         valid=ifelse(pair1_comp=="intra", str_detect(pair2, pair1_spp1), str_detect(pair1, pair2_spp1))) %>% 
  filter(pair1_comp!=pair2_comp, valid==TRUE) %>% 
  mutate(new_pair=ifelse(pair1_comp=="intra", str_c(pair1, pair2, sep=" - "), str_c(pair2, pair1, sep=" - ")),
         spp_mean=ifelse(pair1_comp=="intra", spp_mean, spp_mean*-1),
         .upper=ifelse(pair1_comp=="intra", .upper, .upper*-1),
         .lower=ifelse(pair1_comp=="intra", .lower, .lower*-1)) %>% 
  ggplot(aes(y = new_pair, x = spp_mean, color=newsite)) +
  geom_pointintervalh(position = position_dodgev(height=0.5)) +
  geom_vline(xintercept=0, linetype="dashed")

library(modelr)
alldata %>%
  data_grid(pair, dist=-1, sizeratio=0) %>%
  add_predicted_draws(mod11, allow_new_levels=TRUE, re_formula = ~(1|pair)) %>%
  median_qi() %>% 
  ggplot(aes(x = .prediction, y = pair)) +
  geom_pointintervalh()






### Part 2: through time----
## Data----

pearson_t <- read.csv("Processed Data/synchrony_decadal_pearson.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate(pair=paste(Species.x, Species.y, sep="-"),
         pair=ifelse(pair=="PJ-PP", "PP-PJ", ifelse(pair=="PL-PP", "PP-PL", pair)))

precip <- read.csv("Processed Data/precip_temp_spei.csv")

timeseq <- seq(1910, 2020, by=10)

precip_summary <- expand_grid(precip, timeseq) %>% 
  filter(Year-timeseq>-29, Year<=timeseq) %>% 
  group_by(Site, siteno, Neighborhood, hilo, timeseq) %>% 
  summarise_at(vars(total_ppt_mm, mean_temp_C, spei12), list(mean=mean, sd=sd))


## join all data ----
alldata_t <- pearson_t %>% 
  left_join(tree_data, by=c("tree2"="tree.uniqueID")) %>% 
  mutate(sizeratio=ifelse(DBH.x>DBH.y, DBH.x/DBH.y, DBH.y/DBH.x)) %>% 
  mutate_at(vars(dist, sizeratio, comp_BA), scale) %>% 
  left_join(precip_summary, by=c('Site.x'='Site', 'decade'='timeseq'))

library(ggridges)

ggplot(precip_summary, aes(timeseq, mean_temp_C_mean, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(precip_summary, aes(timeseq, mean_temp_C_sd, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(precip_summary, aes(timeseq, total_ppt_mm_mean, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(precip_summary, aes(timeseq, total_ppt_mm_sd, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(precip_summary, aes(timeseq, spei12_mean, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(precip_summary, aes(timeseq, spei12_sd, group=interaction(Site, Neighborhood))) +
  geom_line(aes(color=Site))

ggplot(alldata_t, aes(x=pearson_r, y=as.factor(decade))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~pair) +
  theme_test() +
  ggtitle("Pearson 30 year blocks, 10 year lag") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("30 year block (year - 29)")

summary <- alldata_t %>% 
  group_by(pair, decade) %>% 
  summarise(pearson_r=median(pearson_r))

ggplot(summary, aes(decade, pearson_r, color=pair)) +
  # geom_point() +
  geom_line()

summary <- alldata_t %>% 
  group_by(Site.x, decade) %>% 
  summarise(pearson_r=median(pearson_r))

ggplot(summary, aes(decade, pearson_r, color=Site.x)) +
  # geom_point() +
  geom_line()

summary <- alldata_t %>% 
  group_by(pair, decade, Site.x) %>% 
  summarise(pearson_r=median(pearson_r),
            spei12_sd=median(spei12_sd, na.rm=T))

ggplot(summary, aes(decade, pearson_r, group=interaction(pair,Site.x), color=pair)) +
  # geom_point() +
  geom_line() +
  facet_wrap(~Site.x)


ggplot(summary, aes(spei12_sd, pearson_r, color=Site.x)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~pair)


mod_s_t <- brm(bf(pearson_r ~ s(decade)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
conditional_effects(mod_s_t)
plot(conditional_effects(mod_s_t), points=TRUE)

mod_s_t2 <- brm(bf(pearson_r ~ s(decade) + (decade|pair)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
conditional_effects(mod_s_t)

mod_s_t3 <- brm(bf(pearson_r ~ pair + s(decade, by=pair)), data=alldata_t, cores=4, control = list(adapt_delta=0.99))

mod11_t<- brm(pearson_r~decade + (decade|pair/Region:hilo), 
              data=alldata_t, cores=4, control = list(adapt_delta=0.99))
mod11_t<- brm(pearson_r~dist + sizeratio + decade + (1|pair/Region:hilo) + (dist|Region:hilo) + (1|Region:hilo:Neighborhood.x), 
              data=alldata_t, cores=4, control = list(adapt_delta=0.99))

saveRDS(mod_s_t3, "saved models/spline.rds")


clim1 <- brm(pearson_r~spei12_mean + (spei12_mean|Site.x:pair), data=alldata_t, cores=4, control = list(adapt_delta=0.99))
