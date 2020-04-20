## Script to calculate group level ICC's (sensu Shestakova)
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

sub <- filter(rings_field, Site!="ic")

mixed_mod <- brm(raw_growth ~ 1 + (1|Year:Species:Region:hilo), data=sub, cores=4)


sig_vars <- mixed_mod %>% 
  spread_draws(sigma, `r_Year:Species:Region:hilo`[group,]) %>% 
  median_qi() 

sig_vars2 <- sig_vars %>% 
  mutate(Year=str_sub(group,1,4),
         Species=str_sub(group, 6,7),
         Region=str_split_fixed(group, "_", 4)[,3],
         elev=str_split_fixed(group, "_", 4)[,4]) 
  # group_by(Species) %>% 
  # mutate(species_sd=sd(`r_Year:Species:Site`)) %>% 
  # ungroup() %>% 
  # group_by(Site) %>% 
  # mutate(site_sd=sd(`r_Year:Species:Site`))

sig_vars_sp <- sig_vars2 %>% 
  group_by(Species) %>% 
  summarize(sp_sigma = sd(`r_Year:Species:Region:hilo`),
            n=length(sigma),
            sigma= mean(sigma)) 

sig_vars_sp2 <- sig_vars_sp %>% 
  rename(Species1=Species, sp1_sigma=sp_sigma, sp1_n=n) %>% 
  expand_grid( sig_vars_sp[,1:3]) %>% 
  mutate(group_sig=sqrt((sp1_sigma^2*(sp1_n-1)+sp_sigma^2*(n-1))/(sp1_n+n-2)),
         a_hat=ifelse(Species1==Species, group_sig^2/(group_sig^2+sigma^2), 
                      group_sig^2/sqrt((sp1_sigma^2+sigma^2)*(sp_sigma^2+sigma^2))))


sig_vars_loc <- sig_vars2 %>% 
  group_by(Region) %>% 
  summarize(loc_sigma = sd(`r_Year:Species:Region:hilo`),
            n=length(sigma),
            sigma= mean(sigma)) 

sig_vars_loc2 <- sig_vars_loc %>% 
  rename(Region1=Region, loc1_sigma=loc_sigma, loc1_n=n) %>% 
  expand_grid(sig_vars_loc[,1:3]) %>% 
  mutate(group_sig=sqrt((loc1_sigma^2*(loc1_n-1)+loc_sigma^2*(n-1))/(loc1_n+n-2)),
         a_hat=ifelse(Region1==Region, group_sig^2/(group_sig^2+sigma^2), 
                      group_sig^2/sqrt((loc1_sigma^2+sigma^2)*(loc_sigma^2+sigma^2))))

sig_vars_elev <- sig_vars2 %>% 
  group_by(Region, elev) %>% 
  summarize(elev_sigma = sd(`r_Year:Species:Region:hilo`),
            n=length(sigma),
            sigma= mean(sigma)) 

sig_vars_elev2 <- sig_vars_elev %>% 
  rename(elev1=elev, Region1=Region, elev1_sigma=elev_sigma, elev1_n=n) %>% 
  expand_grid(sig_vars_elev[,1:4]) %>% 
  mutate(group_sig=sqrt((elev1_sigma^2*(elev1_n-1)+elev_sigma^2*(n-1))/(elev1_n+n-2)),
         a_hat=ifelse(elev1==elev&Region1==Region, group_sig^2/(group_sig^2+sigma^2), 
                      group_sig^2/sqrt((elev1_sigma^2+sigma^2)*(elev_sigma^2+sigma^2))))

sig_vars_all <- sig_vars2 %>% 
  group_by(Species, Region, elev) %>% 
  summarize(SRE_sigma = sd(`r_Year:Species:Region:hilo`),
            n=length(sigma),
            sigma= mean(sigma)) 

sig_vars_all2 <- sig_vars_all %>% 
  rename(Species1=Species, elev1=elev, Region1=Region, SRE1_sigma=SRE_sigma, SRE1_n=n) %>% 
  expand_grid(sig_vars_all[,1:5]) %>% 
  mutate(group_sig=sqrt((SRE1_sigma^2*(SRE1_n-1)+SRE_sigma^2*(n-1))/(SRE1_n+n-2)),
         a_hat=ifelse(Species1==Species&elev1==elev&Region1==Region, group_sig^2/(group_sig^2+sigma^2), 
                      group_sig^2/sqrt((SRE1_sigma^2+sigma^2)*(SRE_sigma^2+sigma^2))))






dist_long <- read.csv("Processed Data/dist_long.csv") %>% 
  mutate(focalID=str_replace(focalID, "_", ""),
         periphID=str_replace(periphID, "_", "")) %>% 
  rename(tree1=focalID, tree2=periphID)

synchrony_dat <- read.csv("Processed Data/synchrony_multivar.csv") %>% 
  left_join(dist_long)
sync_dist <- filter(synchrony_dat, !is.na(dist))

distmod <- brm(.value ~ pair*dist, data=sync_dist, cores=4)


install.packages("performance")
library("performance")

h <- paste("sd_Region:hilo:Neighborhood__Intercept^2 / (sd_Region:hilo:Neighborhood__Intercept^2 + sigma^2) = 0")

(hyp2 <- hypothesis(mixed_mod2, h, class = NULL))

PPD_0 <- predict(mixed_mod2, re.form = NA)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)
