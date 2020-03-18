dist_mat <- read.csv("Raw data/dist mat 20m.csv")
precip <- read.csv("Processed Data/precip_temp_spei.csv")
rings_field <- read.csv("Processed Data/matched_rings_field.csv") %>% 
  filter(siteno==11) %>% 
  filter(Neighborhood==1) %>% 
  left_join(precip) %>% 
  group_by(Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID))))

grow_mat <- rings_field %>% 
  pivot_wider(id_cols=c("Year", "spei12", "mean_temp_C", "total_ppt_mm"), names_from="tree.uniqueID", values_from="spline_growth", names_prefix="t") %>% 
  na.omit() %>% 
  mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)

names(grow_mat) <- c("year", "spei12", "mean_temp_C","total_ppt_mm","a1", "a2", "a3", "a4", "b5", "b6", "b7", "b8", "c9", "c10", "c11", "c12")
a <- c(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12)
fit <- brm(bf(mvbind(a) ~ 0, sigma~0), data=grow_mat, cores=4)
fit2 <- brm(bf(mvbind(a1, a2, a3, a4, b5, b6, b7, b8, c9, c10, c11, c12) ~ 0 + spei12 + mean_temp_C), data=grow_mat, cores=4)

vars <- fit %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )

test <- vars %>% 
  group_by(pair) %>% 
  summarize_all(median)
ggplot(test,aes(x=.value, y=pair)) +
  geom_point()

vars2 <- fit2 %>% 
  gather_draws(`rescor.*`, regex=T) %>% 
  median_qi() %>% 
  mutate(
    temp=str_split_fixed(.variable, "_", 3)[,3],
    tree1=str_split_fixed(temp, "__", 2)[,1],
    tree2=str_split_fixed(temp, "__", 2)[,2],
    pair=paste0(str_sub(tree1, 1,1), str_sub(tree2, 1,1))
  )


test2 <- vars2 %>% 
  group_by(pair) %>% 
  summarize_all(median)
ggplot(test2,aes(x=.value, y=pair)) +
  geom_point()
  

plot(vars$.value, vars2$.value)
abline(0,1)
  
ggplot(vars, aes(x=.value, y=.variable, upper=.upper, lower=.lower)) +
  geom_pointintervalh() +
  facet_wrap(~pair)

ggplot(vars2, aes(x=.value, y=.variable, upper=.upper, lower=.lower)) +
  geom_pointintervalh() +
  facet_wrap(~pair)


vars2 <-  arrange(vars2, tree1, tree2)

ggplot(vars2, aes(tree1, tree2, fill= .value)) + 
  geom_tile()

ggplot(vars, aes(tree1, tree2, fill= .value)) + 
  geom_tile()




### Species responses ----
rings_field <- read.csv("Processed Data/matched_rings_field.csv") %>% 
  filter(siteno==11) %>% 
  left_join(precip) %>% 
  group_by(Neighborhood, Species) %>% 
  mutate(treeID=as.numeric(droplevels(as.factor(tree.uniqueID)))) %>% 

comp_spp <- rings_field %>% 
  select(c(Year, spline_growth, Site, Species, Neighborhood, tree.uniqueID, treeID, mean_temp_C, spei12)) %>% 
  pivot_wider(id_cols=-c(spline_growth, Species, tree.uniqueID),
              names_from="Species", values_from=c("tree.uniqueID","spline_growth")) %>% 
  mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)

bfa <- bf(spline_growth_AC | mi() ~ 0 + (1|Neighborhood:treeID))
bfb <- bf(spline_growth_PJ | mi() ~ 0 + (1|Neighborhood:treeID))
bfc <- bf(spline_growth_PL | mi() ~ 0 + (1|Neighborhood:treeID))

fit3 <- brm(data = comp_spp,
                mvbf(bfa, bfb, bfc, rescor = TRUE), cores=4)

fit3 <- brm(bf(mvbind(spline_growth_AC, spline_growth_PJ, spline_growth_PL | mi()) ~ 0 ), data=comp_spp, cores=4)
  

comp_long <- rings_field %>% 
  select(c(Year, spline_growth, Site, Neighborhood, tree.uniqueID, mean_temp_C, spei12)) %>% 
  filter()
  mutate_at(.vars = vars(-Year, spei12, mean_temp_C, total_ppt_mm), .funs = scale)



### Compare with pearson----
grow_dat <- dplyr::select(rings_field, c("Year", "spline_growth", "tree.uniqueID", "Species"))


sync_data <- read.csv("Processed Data/synchrony_dat.csv") %>% 
  mutate(spp_pair=ifelse(Species.x==Species.y, "same", "different"),
         int_group=paste0(Species.x, Species.y),
         s.no1=str_sub(tree1, 1,2),
         s.no2=str_sub(tree2, 1,2),
         nb1=str_sub(tree1, 5,5),
         nb2=str_sub(tree2, 5,5)) %>% 
  filter(s.no1==s.no2&nb1==nb2) %>% 
  filter(Site.x=='bm') %>% 
  filter(nb1==1) 

test2 <- sync_data %>% 
  group_by(int_group) %>% 
  summarize_at(.vars=vars(pearson_r), median)

ggplot(test2, aes(pearson_r, int_group)) +
  geom_point()


library(Hmisc)
library(reshape2)
tree_sub <- filter(rings_field, Year>=1952) %>% 
    select(tree.uniqueID, Year, spline_growth) %>% 
    pivot_wider(names_from="tree.uniqueID", values_from="spline_growth") 
  
  if(dim(tree_sub)[1]==0){    
    cordf <- NA
  } else {
    treemat <-tree_sub %>% 
      select(-Year) %>%  as.matrix()
    
    cormat <- rcorr(treemat, type='pearson')$r
    cormat[upper.tri(cormat)] <- NA
    
    cordf <- melt(cormat) %>% filter(value!=1&!is.na(value)) %>% 
      rename(tree1=Var1, tree2=Var2, pearson_r=value) 
  }

sub <- select(rings_field, c(tree.uniqueID, Site, Species, Region, hilo, hegyi, DBH)) %>% 
  group_by(tree.uniqueID) %>% 
  summarise_all(first)

cordf2 <- cordf %>% 
  na.omit() %>% 
  mutate_at("tree1", as.character) %>% 
  left_join(sub, by=c("tree1"="tree.uniqueID")) %>% 
  left_join(select(sub, c(tree.uniqueID, Site, Species)), by=c("tree2"="tree.uniqueID")) %>% 
  mutate(int_group=paste0(Species.x, Species.y))

test2 <- cordf2 %>% 
  group_by(int_group) %>% 
  summarize_at(.vars=vars(pearson_r), median)

ggplot(test2, aes(pearson_r, int_group)) +
  geom_point()

uni <- brm(a1~1+a2, data=grow_mat, cores=4)
cor(grow_mat$a1,grow_mat$a2)

ggplot(grow_mat, aes(a1, a2)) +
  geom_point()


### Covariance nonsense----
sds <- fit2 %>% 
  gather_draws(`sigma.*`, regex=T) %>% 
  mutate_at(vars(.value), exp) %>% 
  median_qi() %>% 
  mutate(
    tree1=str_split_fixed(.variable, "_", 2)[,2]
  )


vars2_wide <- vars2 %>%
  pivot_wider(id_cols="tree1", names_from="tree2", values_from=".value") %>% 
  left_join(dplyr::select(sds, c(tree1,.value)), by="tree1")


test <- vars2_wide %>% 
  mutate(a2=a2*sds$.value[which(sds$tree1=="a2")]*.value,
         a3=a3*sds$.value[which(sds$tree1=="a3")]*.value,
         a4=a4*sds$.value[which(sds$tree1=="a4")]*.value,
         b5=b5*sds$.value[which(sds$tree1=="b5")]*.value,
         b6=b6*sds$.value[which(sds$tree1=="b6")]*.value,
         b7=b7*sds$.value[which(sds$tree1=="b7")]*.value,
         b8=b8*sds$.value[which(sds$tree1=="b8")]*.value,
         c10=c10*sds$.value[which(sds$tree1=="c10")]*.value,
         c11=c11*sds$.value[which(sds$tree1=="c11")]*.value,
         c12=c12*sds$.value[which(sds$tree1=="c12")]*.value,
         c9=c9*sds$.value[which(sds$tree1=="c9")]*.value,
         .value=NULL) %>% 
  pivot_longer(cols=-"tree1", names_to="tree2", values_to="covar", values_drop_na=T)

ggplot(test, aes(tree1, tree2, fill= covar)) + 
  geom_tile()

cor_cov <- left_join(vars2, test)

ggplot(cor_cov,(aes(.value, covar))) +
  geom_point() +
  geom_abline(intercept = 0, slope=1)

## Long form intra vs inter nonsense----
str(rings_field)

spptable <- dplyr::select(rings_field, c("tree.uniqueID", "Species", "treeID")) %>% 
  group_by(tree.uniqueID) %>% 
  summarize_all(first)

comp_long <- rings_field %>% 
  select(c(Year, spline_growth, Site, Neighborhood, tree.uniqueID, mean_temp_C, spei12)) %>% 
  pivot_wider(id_cols=-c("spline_growth", "tree.uniqueID"), 
              names_from="tree.uniqueID", values_from=c("spline_growth")) %>% 
  na.omit() %>% 
  pivot_longer(cols=-c(Year, Site, Neighborhood, mean_temp_C, spei12), 
               names_to="tree.uniqueID", values_to="spline_growth") %>% 
  left_join(spptable) 

int_long <- comp_long %>% 
  rename(comp.ID=tree.uniqueID, comp_growth=spline_growth, comp_sp=Species) %>% 
  full_join(comp_long) %>% 
  mutate(comp_type=ifelse(Species==comp_sp, "intra", "inter")) %>% 
  pivot_wider(id_cols=-c(comp_growth, comp_sp, comp_type), 
              names_from="comp_type", 
              values_from=c("comp_growth", "comp_sp"))
