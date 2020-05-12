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


dist <- read.csv("Processed Data/dist_focal.csv") %>% 
  filter(focalID%in%unique(neigh1$tree.uniqueID)) %>% 
  pivot_wider(id_cols = focalID, names_from = periphID, values_from = dist) %>% 
  select(sort(tidyselect::peek_vars()))

dist_mat <- as.matrix(dplyr::select(dist, -focalID))
rownames(dist_mat) <- dist$focalID


mod <- brm(spline_growth~ 1 + (1|tree.uniqueID) + (1|Species), data=neigh1, cores=4,)


raw_loc <- read.csv("Raw data/all_loc.csv")  

neigh1_loc <- raw_loc %>% 
  filter(tree.uniqueID%in%unique(neigh1$tree.uniqueID))

neigh1 <- left_join(neigh1, neigh1_loc) %>% 
  filter(Year>1950) 

neigh1$nb <- as.factor(neigh1$Neighborhood)

str(neigh1)



gp_mod <- brm(data=neigh1,
              spline_growth ~ 1 + gp(X, Y, gr=TRUE, scale=T),
              prior = c(prior(normal(0, 10), class = Intercept)),
              iter = 1e4, warmup = 2000, chains = 4, cores = 4,
              seed = 13,
              control = list(adapt_delta = 0.999,
                             max_treedepth = 12))

gp_mod0 <- brm(data=neigh1,
               spline_growth ~ 1 ,
               prior = c(prior(normal(0, 10), class = Intercept)),
               iter = 1e4, warmup = 2000, chains = 4, cores = 4,
               seed = 13,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 12))

l_gp <- loo(gp_mod)
l_0 <- loo(gp_mod0)


fit2 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD,
            autocor = cor_errorsar(COL.nb),
            chains = 2, cores = 2)

tree_weight <- 1/dist_mat

diag(tree_weight) <- 0


car_mod <- brm(data=neigh1,
               spline_growth ~ 1 , autocor = cor_icar(tree_weight),
               prior = c(prior(normal(0, 10), class = Intercept)),
               iter = 1e4, warmup = 2000, chains = 4, cores = 4,
               seed = 13,
               control = list(adapt_delta = 0.999,
                              max_treedepth = 12))

loo_model_weights(list(l_gp, l_0))

# for `sample_n()`
set.seed(13)

# wrangle

post_gp <- posterior_samples(gp_mod)
post_gp %>% 
  transmute(iter  = 1:n(),
            etasq = sdgp_gpXY^2,
            rhosq = lscale_gpXY^2 * .5) %>% 
  sample_n(100) %>% 
  expand(nesting(iter, etasq, rhosq),
         x = seq(from = 0, to = 10, by = 0.1)) %>% 
  mutate(covariance = etasq * exp(-rhosq * x^2)) %>% 
  
  # plot
  ggplot(aes(x = x, y = covariance)) +
  geom_line(aes(group = iter),
            size = 1/4, alpha = 1/4, color = "darkblue") +
  stat_function(fun = function(x) median(post_gp$sdgp_gpXY)^2 *
                  exp(-median(post_gp$lscale_gpXY)^2 *.5 * x^2),
                color = "darkblue", size = 1.1) +
  scale_x_continuous("distance", expand = c(0, 0),
                     breaks = seq(from = 0, to = 10, by = 1)) +
  # coord_cartesian(xlim = 0:5,
  # ylim = 0:1) +
  theme_classic()



dist_mat[is.na(dist_mat)]=0

k <- matrix(0, nrow = 12, ncol = 12)
for (i in 1:12)
  for (j in 1:12)
    k[i, j] <- median(post_gp$sdgp_gpXY^2) * 
  exp(-median(post_gp$lscale_gpXY)^2 *.5 *dist_mat[i, j]^2)

diag(k) <- median(post_gp$sdgp_gpXY^2) + 0.01

k %>% round(4)

rho <- round(cov2cor(k), 2)
