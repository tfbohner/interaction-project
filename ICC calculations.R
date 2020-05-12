## Script to calculate group level ICC's Adlay/Shestakova package and brms equivalent
library(tidyverse)
library(brms)
library(tidybayes)
library(ggstance)
library(performance)

## Shestakova package
library(DendroSync)

## Data ---
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)


## Dendrosync----
## Example neighborhood
neigh1 <- rings_field %>% 
  filter(siteno==11) %>%
  # filter(Neighborhood==1) %>%
  # filter(Species%in%c("AC", "PJ")) %>% 
  group_by(tree.uniqueID) %>% 
  mutate(min=min(Year),
         max=max(Year)) %>% 
  ungroup() %>% 
  filter(min<1930, max>2010) %>%
  # filter(Year>=max(min), Year<min(max))
  # mutate(treeID=as.factor(str_c("t", tree.uniqueID))) %>% 
  as.data.frame() %>% 
  mutate_if(is.factor, droplevels)

## all time
## homogeneous error
ModHm <- dendro.varcov(spline_growth ~ tree.uniqueID, varTime = "Year", varGroup = "Species",
                       data = neigh1, homoscedastic = TRUE)

summary(ModHm)
gen.aSE(ModHm$mBE)#Broad evaluation
gen.aSE(ModHm$mUN)#Unstructured

sync(ModHm, modname = "mBE")
sync(ModHm, modname = "mCS")
sync(ModHm, modname = "mUN")


sync.plot(sync(ModHm, modname = "mUN"))

## heterogenous error
ModHet <- dendro.varcov(spline_growth ~ tree.uniqueID, varTime = "Year", varGroup = "Species",
                       data = neigh1, homoscedastic = FALSE)

summary(ModHet)
gen.het.aSE(ModHet$mHeNE)#
gen.het.aSE(ModHet$mHeUN)#Unstructured


sync(ModHet, modname = "mHeNE")
sync(ModHet, modname = "mHeCS")
sync(ModHet, modname = "mHeUN")

sync.plot(sync(ModHet, modname = "mHeUN"))



## through time
mBE.trend <- sync.trend(spline_growth ~ tree.uniqueID, varTime = "Year", varGroup = "Species",
                        data = neigh1,  window = 30, lag = 5,
                        null.mod = TRUE)

sync.trend.plot(mBE.trend)# Broad evaluation synchrony linechart

## homogeneous error
spp.trend <- sync.trend(spline_growth ~ tree.uniqueID, varTime = "Year", varGroup = "Species",
                        data = neigh1,  window = 30, lag = 5,
                        null.mod = FALSE, homoscedastic = TRUE)

sync.trend.plot(spp.trend)


spp2.trend <- sync.trend(spline_growth ~ tree.uniqueID, varTime = "Year", varGroup = "Species",
                        data = neigh1,  window = 30, lag = 5,
                        null.mod = FALSE, homoscedastic = TRUE, between.group=TRUE)

sync.trend.plot(spp2.trend)




## Brms version----
## Species fixed effect---- 
## I think this is the equivalent to the unstructured model???
icc_test <- brm(spline_growth ~ Species + (Species||Year), data=neigh1, cores=4)

performance::icc(icc_test)

brms::VarCorr(icc_test)

## the overall ungrouped icc is similar but not the same (.19 versus .20) which makes sense because in this model
## year effects are partially pooled whereas in the dendrosync package they are not.

## Extract species and year effects
sig_vars <- icc_test %>% 
  spread_draws(sigma, `r_Year`[Year, Species]) %>%
  # median_qi()
  ungroup() %>% 
  mutate(Species=ifelse(Species=="Intercept", "SpeciesAC", Species)) 

sig_vars_sp <- sig_vars %>% 
  ungroup() %>% 
  group_by(Species, .chain, .iteration, .draw) %>%
  summarise(sp_sigma = sd(`r_Year`),
            n=length(sigma),
            sigma= mean(sigma),
            icc=sp_sigma^2/(sp_sigma^2 + sigma^2))

## Intra-Class Correlations based on posterior, within species
intra_summ <- sig_vars_sp %>% 
  ungroup() %>% 
  dplyr::select(c(Species, icc)) %>% 
  group_by(Species) %>% 
  median_qi()

ggplot(intra_summ, aes(y=Species, x=icc)) +
  geom_pointintervalh()

## Inter-class correlations based on posterior, between species

inter_grid <- sig_vars_sp %>% 
  rename_all(.funs = function(x)str_c(x, ".x")) %>% 
  expand_grid(sig_vars_sp[,1:8]) %>% 
  filter(.draw.x==.draw)
  
sig_vars_sp2 <- inter_grid %>% 
  mutate(group_sig = sqrt((sp_sigma.x^2*(n.x-1)+sp_sigma^2*(n-1))/(n.x+n-2)),  ## I've tried calculating the group level variance a few different ways...
         a_hat = ifelse(Species.x==Species, group_sig^2/(group_sig^2+sigma^2), 
                        group_sig^2/sqrt((sp_sigma.x^2+sigma^2)*(sp_sigma^2+sigma^2))))

inter_summ <- sig_vars_sp2 %>% 
  ungroup() %>% 
  dplyr::select(c(Species.x, Species, a_hat)) %>% 
  group_by(Species.x, Species) %>% 
  median_qi()

ggplot(inter_summ, aes(y=Species, x=a_hat, color=Species.x)) +
  geom_pointintervalh(position = position_dodgev(height=0.3))


icc_test <- brm(spline_growth ~ 1 + (1|Species/Year) +(1|tree.uniqueID), data=neigh1, cores=4)


##### Way different results from either of the dendrosync calculations
## can include correlations:
icc_test2 <- brm(spline_growth ~ Species + (Species|Year), data=neigh1, cores=4)
brms::VarCorr(icc_test2)

## or we can model Species hierarchically too:
icc_test3 <- brm(spline_growth ~ (1|Species:Year), data=neigh1, cores=4, control = list(adapt_delta=0.99))
brms::VarCorr(icc_test3)

## ICC based on posterior
sig_vars <- icc_test3 %>% 
  spread_draws(sigma, `r_Species:Year`[group,]) %>%
  mutate(Year=str_sub(group,4,7),
         Species=str_sub(group, 1,2)) 

sig_vars_sp <- sig_vars %>% 
  ungroup() %>% 
  group_by(Species, .chain, .iteration, .draw) %>%
  summarise(sp_sigma = sd(`r_Species:Year`),
            n=length(sigma),
            sigma= mean(sigma),
            icc=sp_sigma^2/(sp_sigma^2 + sigma^2)) 

## Intra-Class Correlations based on posterior, within species
intra_summ <- sig_vars_sp %>% 
  ungroup() %>% 
  dplyr::select(c(Species, icc)) %>% 
  group_by(Species) %>% 
  median_qi()

ggplot(intra_summ, aes(y=Species, x=icc)) +
  geom_pointintervalh()

inter_grid <- intra_summ[,1] %>%
  rename(Species.x=Species) %>% 
  expand_grid(intra_summ[,1] ) %>% 
  rename(Species.y=Species) %>% 
  mutate(pair=str_c(Species.x, Species.y, sep = "-")) %>% 
  pivot_longer(cols=c(Species.x, Species.y), names_to = "xy", values_to="Species") %>% 
  left_join(sig_vars[,3:8]) %>% 
  group_by(pair, .draw) %>% 
  mutate(sigma = mean(sigma),
            group_sigma = sd(`r_Species:Year`)) %>% 
  group_by(pair, xy, Species, .draw) %>% 
  mutate(sp_sigma = sd(`r_Species:Year`)) %>% 
  summarise_all(mean)

  
sig_vars_sp2 <- inter_grid %>% 
  pivot_wider(id_cols = c(pair, .draw, sigma, group_sigma), names_from=xy, values_from=c(Species, sp_sigma)) %>% 
  ungroup() %>% 
  group_by(pair, .draw) %>% 
  summarise(a_hat=ifelse(Species_Species.x==Species_Species.y,
                       sp_sigma_Species.x^2/(sp_sigma_Species.x^2+sigma^2), 
                       group_sigma^2/sqrt((sp_sigma_Species.x^2+sigma^2)*(sp_sigma_Species.y^2+sigma^2))))

inter_summ <- sig_vars_sp2 %>% 
  ungroup() %>% 
  dplyr::select(c(pair, a_hat)) %>% 
  group_by(pair) %>% 
  median_qi()

ggplot(inter_summ, aes(y=pair, x=a_hat)) +
  geom_pointintervalh()





## ICC based on predicted posterior ----
## equivalent to performance::icc() but cannot figure out if I'm doing the group covariances correctly
dat <- neigh1

PPD <- brms::posterior_predict(icc_test2, re.form = ~(Species|Year), summary = FALSE, newdata=dat)
var_total <- apply(PPD, MARGIN = 1, FUN = stats::var)

PPD_0 <- brms::posterior_predict(icc_test2, re.form = NA, summary = FALSE, newdata=dat)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

var_icc_intra <- var_rand_intercept / var_total
var_residual <- var_total - var_rand_intercept

mean(var_residual)
1 - mean(var_icc_intra) 

## is it just as simple as subsetting speceis or groups of species? This returns more similar values to modHm unstructured...
## there are no other predictors in the model other than species and year so actual data input shouldn't matter except for equal year representation?

dat <- filter(neigh1, Species=="AC") ## or
dat <- mutate(neigh1,
              Species="AC")
PPD <- brms::posterior_predict(icc_test2, re.form = ~(Species|Year), summary = FALSE, newdata=dat)
var_total <- apply(PPD, MARGIN = 1, FUN = stats::var)

PPD_0 <- brms::posterior_predict(icc_test2, re.form = NA, summary = FALSE, newdata=dat)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

var_icc_intra <- var_rand_intercept / var_total
var_residual <- var_total - var_rand_intercept

mean(var_residual)
1 - mean(var_icc_intra) 

dat <- mutate(neigh1,
              Species="PJ")
PPD <- brms::posterior_predict(icc_test2, re.form = ~(Species|Year), summary = FALSE, newdata=dat)
var_total <- apply(PPD, MARGIN = 1, FUN = stats::var)

PPD_0 <- brms::posterior_predict(icc_test2, re.form = NA, summary = FALSE, newdata=dat)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

var_icc_intra <- var_rand_intercept / var_total
var_residual <- var_total - var_rand_intercept

mean(var_residual)
1 - mean(var_icc_intra) 

dat <- mutate(neigh1,
              Species="PL")
PPD <- brms::posterior_predict(icc_test2, re.form = ~(Species|Year), summary = FALSE, newdata=dat)
var_total <- apply(PPD, MARGIN = 1, FUN = stats::var)

PPD_0 <- brms::posterior_predict(icc_test2, re.form = NA, summary = FALSE, newdata=dat)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

var_icc_intra <- var_rand_intercept / var_total
var_residual <- var_total - var_rand_intercept

mean(var_residual)
1 - mean(var_icc_intra) 


sub <- filter(neigh1, Species%in%c("AC", "PJ"))
PPD <- brms::posterior_predict(icc_test2, re.form = ~(Species|Year), summary = FALSE, newdata=sub)
var_total <- apply(PPD, MARGIN = 1, FUN = stats::var)

PPD_0 <- brms::posterior_predict(icc_test2, re.form = NA, summary = FALSE, newdata=sub)
var_rand_intercept <- apply(PPD_0, MARGIN = 1, FUN = stats::var)

var_icc_intra <- var_rand_intercept / var_total
var_residual <- var_total - var_rand_intercept

mean(var_residual)
1 - mean(var_icc_intra) 


