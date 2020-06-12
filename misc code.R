## miscellaneous code from exploring different methods
library(tidyverse)
library(tidybayes)
library(Hmisc)
library(cocor)
library(ggfortify)
library(corrr)
library(codyn)

## Data ---
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)


## Example neighborhood----
neigh1 <- rings_field %>% 
  filter(siteno==11) 

grow1 <- neigh1 %>% 
  group_by(Site, Species, Neighborhood, Year) %>% 
  summarise(spline_growth=mean(spline_growth))

neigh1_sync <- synchrony(df=grow1,
                          time.var="Year",
                          species.var="Species",
                          abundance.var="spline_growth",
                          replicate.var = "Neighborhood")

neigh1_variance_ratio <- variance_ratio(df = grow1, 
                                     species.var = "Species", 
                                     time.var = "Year",
                                     abundance.var = "spline_growth", 
                                     bootnumber = 10, 
                                     replicate.var = "Neighborhood")


count_na <- function(x) sum(is.na(x))
  

grow_mat <- neigh1 %>% 
  pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
  na.omit() %>% ## can include missing data if we want
  mutate(nacount=apply(., 1, count_na)) %>% 
  filter(nacount<4) %>% 
  select(-nacount)
  # mutate_at(.vars = vars(-c(Year, spei12, mean_temp_C, total_ppt_mm)), .funs = scale)


names(grow_mat) <- c("year", "spei12", "mean_temp_C","total_ppt_mm","a1", "a2", "a3", "a4", "b5", "b6", "b7", "b8", "c9", "c10", "c11", "c12")
apply(grow_mat, 2, count_na)

test <- cocor(~a1 + a2 | a1 + b5, as.data.frame(grow_mat), "williams1959", alternative="greater")



mat <- round(rcorr(as.matrix(grow_mat[,5:16]), type='pearson')$r,2)
pr <- as.dist(round(cor(grow_mat[,5:13], method = 'pearson'),2))


mat <- correlate(grow_mat[,5:16])
network_plot(mat)

rearrange(mat)

pc <- prcomp(grow_mat[,5:16], scale=T)

pca.plot <- autoplot(pc, data = grow_mat, colour="spei12")

fviz_pca_var(pc,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

pc <- prcomp(t(grow_mat[,5:16]), scale=T)

fviz_pca_ind(pc, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

R <- cor(grow_mat[,5:16])
r.eigen <- eigen(R)

for (r in r.eigen$values) {
  print(r / sum(r.eigen$values))
}



## Data ---
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)


## Example neighborhood----
neigh1 <- rings_field %>% 
  filter(siteno==8) %>% 
  filter(Neighborhood==1)


dist <- read.csv("Processed Data/dist_focal.csv") %>% 
  filter(focalID%in%unique(neigh1$tree.uniqueID)) %>% 
  pivot_wider(id_cols = focalID, names_from = periphID, values_from = dist) %>% 
  dplyr::select(sort(everything()))

dist_mat <- dplyr::select(dist, -focalID) %>% 
  dplyr::select(sort(everything())) %>% 
  as.matrix()

dist_mat <- dist_mat[,order(colnames(dist_mat))]
rownames(dist_mat) <- dist$focalID

dist_mat[is.na(dist_mat)]=0




mod <- brm(spline_growth~ 1 + (1|tree.uniqueID) + (1|Species), data=neigh1, cores=4)

## Example data----
site1 <- rings_field %>% 
  filter(siteno==8) 


raw_loc <- read.csv("Raw data/all_loc.csv") %>% 
  filter(tree.uniqueID%in%unique(site1$tree.uniqueID))

IDs <- as.character(raw_loc$tree.uniqueID)

nb_weight <- dnearneigh(as.matrix(raw_loc[,2:3]), d1 = 0, d2 = 0.04, longlat=TRUE, row.names = IDs)


car_mod <- brm(spline_growth ~ 1, data=site1, cores=4,
               autocor=cor_car(nb_weight, formula=~1|tree.uniqueID))

oneyear <- filter(site1, Year==2000)
sar_mod <- brm(spline_growth ~ 1, data=oneyear, cores=4,
               autocor=cor_sar(nb_weight, type = "lag"))

car_mod <- brm(spline_growth ~ 1, data=oneyear, cores=4,
               autocor=cor_car(nb_weight))


dist <- read.csv("Processed Data/dist_focal.csv") %>% 
  filter(focalID%in%unique(site1$tree.uniqueID)) %>% 
  pivot_wider(id_cols = focalID, names_from = periphID, values_from = dist) %>% 
  dplyr::select(order(colnames(.))) %>% 
  arrange(focalID)

dist_mat <- dplyr::select(dist, -focalID) %>% 
  as.matrix()

dist_mat[is.na(dist_mat)] <-0
rownames(dist_mat) <- dist$focalID

car_mod <- brm(spline_growth ~ 1, data=site1, cores=4,
               autocor=cor_car(dist_mat, formula=~1|tree.uniqueID))

car_mod <- brm(spline_growth ~ DBH + Species, data=site1, cores=4,
               autocor=cor_car(dist_mat, formula=~1|tree.uniqueID))

nocar <- brm(spline_growth ~ DBH + Species, data=site1, cores=4)

loo_model_weights(list(loo(car_mod), loo(nocar)))


sar_mod <- brm(spline_growth ~ (1|Year), data=site1, cores=4,
               autocor = cor_sar(dist_mat, type = 'lag' ))



dist_mat <- dist_mat[,order(colnames(dist_mat))]
rownames(dist_mat) <- dist$focalID

dist_mat[is.na(dist_mat)]=0

neigh1_loc <- raw_loc %>% 
  filter(tree.uniqueID%in%unique(neigh1$tree.uniqueID))

neigh1 <- left_join(neigh1, neigh1_loc) %>% 
  filter(Year>1950)


m.gp2 <- map2stan(
  alist(
    spline_growth ~ dnorm(mu, sigma),
    mu <- a + g[tree] ,
    g[tree] ~ GPL2( dist_mat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1),
    sigma ~ dunif(0,10)
  ),
  data=list(
    spline_growth=neigh1$spline_growth,
    tree=neigh1$tree.uniqueID,
    dist_mat=dist_mat),
  warmup=2000 , iter=1e4 , chains=4, cores=4)


data(Kline2)
data(islandsDistMatrix)
d <- 
  Kline2 %>%
  mutate(society = 1:10)
d_mat <- islandsDistMatrix

postm.gp <- rethinking::extract.samples(m.gp)[2:5] %>% as.data.frame()

# wrangle
postm.gp %>% 
  transmute(iter  = 1:n(),
            etasq = etasq,
            rhosq = rhosq) %>% 
  sample_n(100) %>% 
  group_by(iter, etasq, rhosq) %>% 
  expand_grid(x = seq(from = 0, to = 10, by = .1)) %>% 
  mutate(covariance = etasq * exp(-rhosq * x^2)) %>% 
  
  # plot
  ggplot(aes(x = x, y = covariance)) +
  geom_line(aes(group = iter),
            size = 1/4, alpha = 1/4, color = "blue") +
  stat_function(fun = function(x) median(postm.gp$etasq) *
                  exp(-median(postm.gp$rhosq) * x^2),
                color = "blue", size = 1.1) +
  scale_x_continuous("distance (m)", expand = c(0, 0),
                     breaks = seq(from = 0, to = 10, by = 1)) +
  theme_test()

k <- matrix(0, nrow = 9, ncol = 9)
for (i in 1:9)
  for (j in 1:9)
    k[i, j] <- median(postm.gp$etasq) * 
  exp(-median(postm.gp$rhosq) * 
        dist_mat[i, j]^2)

diag(k) <- median(postm.gp$etasq) + 0.01

k %>% round(4)

rho <- round(cov2cor(k), 2)

m13.7 <- map2stan(
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + g[society] + bp*logpop,
    g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=d$total_tools,
    logpop=d$logpop,
    society=d$society,
    Dmat=islandsDistMatrix),
  warmup=2000 , iter=1e4 , chains=4)


data(oldcol, package = "spdep")
fit1 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "lag"), 
            data = COL.OLD, data2 = list(COL.nb = COL.nb),
            chains = 2, cores = 2)
summary(fit1)
plot(fit1)

fit2 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "error"), 
            data = COL.OLD, data2 = list(COL.nb = COL.nb),
            chains = 2, cores = 2)

fit2 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD,
            autocor = cor_errorsar(COL.nb),
            chains = 2, cores = 2)
summary(fit2)
plot(fit2)

summary(fit2)
plot(fit2)





library(rethinking)
data(islandsDistMatrix)

# display short column names, so fits on screen
d_mat <- islandsDistMatrix
colnames(d_mat) <- c("Ml", "Ti", "SC", "Ya", "Fi", 
                     "Tr", "Ch", "Mn", "To", "Ha")
round(d_mat, 1)

data(Kline2) # load the ordinary data, now with coordinates

d <- 
  Kline2 %>%
  mutate(society = 1:10)

rm(Kline2)

d %>% glimpse()
b13.7 <- 
  brm(data = d, family = poisson,
      total_tools ~ 1 + gp(lat, lon2) + logpop,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(inv_gamma(2.874624, 0.393695), class = lscale),
                prior(cauchy(0, 1), class = sdgp)),
      iter = 1e4, warmup = 2000, chains = 4, cores = 4,
      seed = 13,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 12))
library(rethinking)
data(islandsDistMatrix)

# display short column names, so fits on screen
d_mat <- islandsDistMatrix
colnames(d_mat) <- c("Ml", "Ti", "SC", "Ya", "Fi", 
                     "Tr", "Ch", "Mn", "To", "Ha")
round(d_mat, 1)

data(Kline2) # load the ordinary data, now with coordinates

d <- 
  Kline2 %>%
  mutate(society = 1:10)

rm(Kline2)

d %>% glimpse()
b13.7 <- 
  brm(data = d, family = poisson,
      total_tools ~ 1 + gp(lat, lon2) + logpop,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(inv_gamma(2.874624, 0.393695), class = lscale),
                prior(cauchy(0, 1), class = sdgp)),
      iter = 1e4, warmup = 2000, chains = 4, cores = 4,
      seed = 13,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 12))

post_b13.7 <- posterior_samples(b13.7)

post_b13.7 %>% 
  transmute(iter  = 1:n(),
            etasq = sdgp_gplatlon2^2,
            rhosq = lscale_gplatlon2^2 * .5) %>% 
  sample_n(100) %>% 
  group_by(iter, etasq, rhosq) %>% 
  expand_grid(x = seq(from = 0, to = 55, by = 1)) %>% 
  mutate(covariance = etasq * exp(-rhosq * x^2)) %>% 
  
  # plot
  ggplot(aes(x = x, y = covariance)) +
  geom_line(aes(group = iter),
            size = 1/4, alpha = 1/4, color = "#EEDA9D") +
  stat_function(fun = function(x) median(post_b13.7$sdgp_gplatlon2)^2 *
                  exp(-median(post_b13.7$lscale_gplatlon2)^2 *.5 * x^2),
                color = "#EEDA9D", size = 1.1) +
  scale_x_continuous("distance (thousand km)", expand = c(0, 0),
                     breaks = seq(from = 0, to = 50, by = 10)) +
  coord_cartesian(xlim = 0:50,
                  ylim = 0:1)




install.packages("synchrony")
library(synchrony)

data(pisco.data)
d=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "sst"))
semiv=vario(data=d)
moran=vario(data=d, type="moran", nrands=100)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv$mean.bin.dist, semiv$vario, xlab="Lag distance (km)", ylab="Semivariance")
plot(moran$mean.bin.dist, moran$vario, xlab="Lag distance (km)", ylab="Moran's I", t="l")
points(moran$mean.bin.dist[moran$pvals >= 0.05], moran$vario[moran$pvals >= 0.05],
       bg="white", pch=21)
points(moran$mean.bin.dist[moran$pvals < 0.05], moran$vario[moran$pvals < 0.05],
       bg="black", pch=21)
abline(h=0, lty=2)
# Compute spatial synchrony
d.upw=subset(pisco.data, select=c("latitude", "longitude", "year", "upwelling"))
d.cov=subset(pisco.data, select=c("latitude", "longitude", "year", "mussel_abund"))
# Reshape the data
d.upw.wide=reshape(data=d.upw, idvar=c("latitude", "longitude"), timevar=c("year"),
                   direction="wide")
d.cov.wide=reshape(data=d.cov, idvar=c("latitude", "longitude"), timevar=c("year"),
                   direction="wide")
# Generate variograms
v.upw=vario(n.bins=12, data=d.upw.wide, type="pearson", extent=1, nrands=999)
v.cov=vario(n.bins=12, data=d.cov.wide, type="pearson", extent=1, nrands=999)
## Fit variograms
v.cov.per=vario.fit(v.cov$vario, v.cov$mean.bin.dist, type="period",
                    start.vals=list(a=1, b=3, c=0))
v.upw.lin=vario.fit(v.upw$vario, v.upw$mean.bin.dist, type="linear")
par(mfrow=c(2,1))
plot(v.cov, xlab="Lag distance (km)", bg.sig="red", col.nonsig="red",
     main="Mussel cover",
     rug=TRUE, ylim=c(-0.3, 0.3))
lines(v.cov$mean.bin.dist, v.cov.per$fit, col="red")
plot(v.upw, xlab="Lag distance (km)", bg.sig="blue", col.nonsig="blue",
     main="Upwelling", rug=TRUE)
lines(v.upw$mean.bin.dist, v.upw.lin$fit, col="blue")


d.tree <- dplyr::select(ungroup(rings_field), c(Year, spline_growth, tree.uniqueID)) %>% 
  group_by(tree.uniqueID) %>% 
  mutate(min=min(Year),
         max=max(Year)) %>% 
  left_join(raw_loc, by=c("tree.uniqueID")) %>% 
  filter(min<1980, max>2015) %>% 
  filter(Year>1980, Year<2015) %>% 
  pivot_wider(id_cols=c(Y, X), names_from = Year, values_from = spline_growth)

v.tree=vario(n.bins=4, data=d.tree, type='pearson', extent=1, nrands=999)

v.tree.lin=vario.fit(v.tree$vario, v.tree$mean.bin.dist, type = 'linear')

plot(v.tree, xlab="Lag distance (km)", bg.sig="blue", col.nonsig="blue",
     main="Tree", rug=TRUE)
lines(v.tree$mean.bin.dist, v.tree.lin$fit, col="blue")


install.packages("spdep")
library(spdep)
boreal_dists <- boreal %>%
  dplyr::select(x,y) %>% 
  dist() %>%
  as.matrix()

#inverse of distance is weight
boreal_dists_weights <- 1/boreal_dists
diag(boreal_dists_weights) <- 0

#turn into a weights list
boreal_w <- mat2listw(boreal_dists_weights)

tree_weight <- 1/dist_mat
diag(tree_weight) <- 0
tree_weight[is.infinite(tree_weight)] <-0

tree_w <- mat2listw(tree_weight)


mod <- lm(data=neigh1, spline_growth~1)
moran.test(site1$spline_growth[neigh1$Year==2000], tree_w)
moran.plot(neigh1$spline_growth[neigh1$Year==2000], tree_w)

library(ncf)

grow_mat <- neigh1 %>% 
  pivot_wider(id_cols=c("tree.uniqueID"), names_from="Year", values_from="spline_growth", names_prefix="y") %>%
  dplyr::select(-tree.uniqueID) %>% 
  as.matrix


locs2 <- loc_data %>% 
  filter(tree.uniqueID%in%unique(neigh1$tree.uniqueID))

test <- Sncf(x=locs2$X, y=locs2$Y, z = grow_mat, na.rm=T)
summary(test)
plot(test)
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv")  %>% 
  left_join(precip) %>% 
  group_by(Site, Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% ## counter id if we want it rather than unique ID
  arrange(Site, Neighborhood, Species, ID)


## Example neighborhood----
neigh1 <- rings_field %>% 
  filter(siteno==8) %>% 
  filter(Neighborhood==1)


dist <- read.csv("Processed Data/dist_focal.csv") %>% 
  filter(focalID%in%unique(neigh1$tree.uniqueID)) %>% 
  pivot_wider(id_cols = focalID, names_from = periphID, values_from = dist) %>% 
  dplyr::select(sort(everything()))

dist_mat <- dplyr::select(dist, -focalID) %>% 
  dplyr::select(sort(everything())) %>% 
  as.matrix()

dist_mat <- dist_mat[,order(colnames(dist_mat))]
rownames(dist_mat) <- dist$focalID

grow_mat <- neigh1 %>% 
  pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
  na.omit() 

mat <- round(rcorr(as.matrix(grow_mat[,5:ncol(grow_mat)]), type='pearson')$r,2)
pr <- as.dist(round(cor(grow_mat[,5:ncol(grow_mat)], method = 'pearson'),2))


mat <- correlate(grow_mat[,5:ncol(grow_mat)])
network_plot(mat)

as.dist(dist_mat)
network_plot(dist_mat)


pearson <- read.csv("Processed Data/synchrony_dat.csv") %>% 
  na.omit() %>% 
  filter(Site.x !="ic") %>% 
  mutate(pair=paste(Species.x, Species.y, sep="-"),
         pair=ifelse(pair=="PJ-PP", "PP-PJ", ifelse(pair=="PL-PP", "PP-PL", pair)))

ggplot(pearson, aes(dist, pearson_r)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x, 2), se=FALSE) +
  facet_grid(Site.x~Neighborhood.x)

test_mod <- brm(pearson_r~dist + (1|pair) + (1|Site.x), data=pearson, cores=4, control = list(adapt_delta=0.99))
marginal_effects(test_mod)







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
