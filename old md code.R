## old markdown code:



dat_post <- pair_mod10 %>% 
  spread_draws(b_Intercept, `r_Site:pair`[comp,]) %>% 
  mutate(comp_med=b_Intercept + `r_Site:pair`) %>% 
  dplyr::select(-c(b_Intercept, `r_Site:pair`)) %>% 
  median_qi() %>% 
  mutate(Site=str_sub(comp, 1,2),
         pair=str_sub(comp, 4,8),
         spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5)) %>% 
  left_join(region, by="Site")


ggplot(dat_post, aes(y=interaction(Region, hilo, lex.order=TRUE), x=comp_med, color=Region)) +
  geom_pointintervalh(position = position_dodgev(-0.4)) +
  scale_color_viridis(discrete=T, option='plasma', direction = -1) +
  facet_grid(spp1~spp2) +
  theme_test() +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  xlab("posterior median intercept") +
  ylab("Site") +
  theme(legend.position = 'none')


pred_dat <- sync_mod_dat %>% 
  group_by(Site, pair) %>% 
  data_grid(cor, dist=round(quantile(sync_mod_dat$dist, probs=c(0.5)), digits=2), 
            sizeratio=round(quantile(sync_mod_dat$sizeratio, probs=c(0.5)), digits=2)) %>% 
  add_fitted_draws(pair_mod10, n=100, re_formula = NULL, allow_new_levels=TRUE) %>% 
  mutate(spp1=str_sub(pair, 1,2),
         spp2=str_sub(pair, 4,5),
         .value=.value) %>% 
  left_join(region, by="Site") 

## check out predicted data 
ggplot(pred_dat, aes(y = interaction(Region, hilo, lex.order=TRUE), x = .value, color=Region)) +
  stat_pointintervalh(.width = 0.95) +
  # geom_violinh(aes(x = .value)) +
  geom_vline(xintercept=0, linetype="dashed")+
  facet_grid(spp1~spp2) +
  scale_color_viridis(discrete=TRUE, option='plasma', direction = -1) +
  geom_hline(yintercept = 0, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  xlab("conditional effect")  +
  ylab("Site") +
  theme(legend.position = 'none')

dat_post <- pair_mod10 %>% 
  gather_draws(b_dist, b_sizeratio)

params <- ggplot(dat_post, aes(x=.value, y=.variable)) +
  geom_halfeyeh() +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_test() +
  ylab("parameter") +
  xlab("estimate")

pred_dat <- sync_mod_dat %>% 
  group_by(Site, pair) %>% 
  data_grid(cor, dist=round(quantile(sync_mod_dat$dist, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), digits=2), 
            sizeratio=round(quantile(sync_mod_dat$sizeratio, probs=c(0.5)), digits=2)) %>% 
  add_fitted_draws(pair_mod10, n=100, re_formula = NULL, allow_new_levels=TRUE) %>% 
  mutate(.value=.value)

dist <- ggplot(pred_dat, aes(x = dist, y = .value)) +
  stat_lineribbon(aes(y = .value), .width = c(.95), alpha=0.2) +
  stat_lineribbon(aes(y = .value), .width = c(.01)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_viridis(discrete=TRUE) +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted correlation (Rij)") +
  theme(legend.position = 'none') +
  xlab("Distance")

pred_dat <- sync_mod_dat %>% 
  group_by(Site, pair) %>% 
  data_grid(cor, dist=round(quantile(sync_mod_dat$dist, probs=c(0.5)), digits=2), 
            sizeratio=round(quantile(sync_mod_dat$sizeratio, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), digits=2)) %>% 
  add_fitted_draws(pair_mod10, n=100, re_formula = NULL, allow_new_levels=TRUE)  %>% 
  mutate(.value=.value)

ratio <- ggplot(pred_dat, aes(x = sizeratio, y = .value)) +
  stat_lineribbon(aes(y = .value), .width = c(.95), alpha=0.2) +
  stat_lineribbon(aes(y = .value), .width = c(.01)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_viridis(discrete=TRUE) +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted correlation (Rij)")  +
  theme(legend.position = 'none') +
  xlab("Size Ratio")

plot_grid(params, plot_grid(dist, ratio, nrow=1), nrow=2)

bern_mod9 <- readRDS("saved models/intra_greater_mods/bern_mod9.rds")
head(dplyr::select(all_cor, c(tree1.x, tree2.x, comp.x, .value.x, tree1.y, tree2.y, comp.y, .value.y, intra_greater)))

dat_post <- bern_mod9 %>% 
  spread_draws(b_Intercept, `r_Site:intra_spp:pair2`[comp,]) %>% 
  mutate(comp_med=inv_logit_scaled(b_Intercept + `r_Site:intra_spp:pair2`)) %>% 
  dplyr::select(-c(b_Intercept, `r_Site:intra_spp:pair2`)) %>% 
  median_qi() %>% 
  # mutate(Region=str_split_fixed(comp, "_", 3)[,1],
  #        intra_spp=str_sub(str_split_fixed(comp, "_", 3)[,3], 1,2),
  #        inter_spp=str_sub(str_split_fixed(comp, "_", 3)[,3], 4,5),
  #        Region=str_replace(Region, "\\.", " ")) %>% 
  mutate(Site=str_sub(comp, 1,2),
         intra_spp=str_sub(str_split_fixed(comp, "_", 3)[,3], 1,2),
         inter_spp=str_sub(str_split_fixed(comp, "_", 3)[,3], 4,5)) %>% 
  left_join(region, by="Site")


ggplot(dat_post, aes(y=reorder(intra_spp, desc(intra_spp)), x=comp_med, 
                     color=interaction(Region, hilo, lex.order=TRUE))) +
  geom_pointintervalh(position = position_dodgev(-0.4)) +
  scale_color_viridis(discrete=T, option='plasma', direction = -1) +
  geom_vline(xintercept=0.5, linetype="dashed") +
  facet_grid(~inter_spp) +
  theme_test() +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab("posterior median proportion of intra > inter correlations") +
  ylab("Intraspecific species pair") +
  theme(legend.position = 'none')

pred_dat <- all_cor %>% 
  group_by(Site, intra_spp, pair2) %>% 
  data_grid(intra_greater) %>% 
  add_fitted_draws(bern_mod9, re_formula = ~(1|Site/intra_spp/pair2), allow_new_levels=TRUE) %>% 
  mutate(.value=inv_logit_scaled(.value),
         spp1=str_sub(pair2, 1,2),
         spp2=str_sub(pair2, 4,5))

## check out predicted data 
ggplot(pred_dat, aes(x = spp2, y = .value, fill=spp2)) +
  geom_violin(aes(y = .value)) +
  facet_grid(~spp1) +
  scale_fill_viridis(discrete=TRUE) +
  # scale_alpha_ordinal(range=c(0.3,1)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  # ylim(-.4, .4) +
  theme_test() +
  ylab("predicted proportion of intra > inter correlations") 

all_cordf_raw <- all_cordf %>% 
  na.omit() %>% 
  mutate_at("tree1", as.character) %>% 
  left_join(sub, by=c("tree1"="tree.uniqueID")) %>% 
  left_join(dplyr::select(sub, c(tree.uniqueID, Site, Neighborhood, Species)), by=c("tree2"="tree.uniqueID")) %>% 
  left_join(dist_mat2) %>% 
  filter(Neighborhood.x==Neighborhood.y, Site.x==Site.y)

write.csv(all_cordf_raw, "Processed Data/synchrony_raw.csv", row.names = F)



ggplot(bigvars_t2, aes(x=.value, y=as.factor(time))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~Site) +
  theme_test() +
  ggtitle("A. concolor") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("20 year block (year +/- 10)")

## subset species----
acdat <- spp_subsetter(bigvars_t2, "ac")
pjdat <- spp_subsetter(bigvars_t2, "pj")
pldat <- spp_subsetter(bigvars_t2, "pl")
ppdat <- spp_subsetter(bigvars_t2, "pp")

cols <- c("self"="#440154FF", "pj"="#3B528BFF", "pl"="#21908CFF", "pp"="#5DC863FF", "ac"="#FDE725FF") ## keep species cols consistent

ggplot(acdat, aes(x=.value, y=as.factor(time), fill=comp_sp)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~comp_sp) +
  theme_test() +
  ggtitle("A. concolor") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("20 year block (year +/- 10)")

ggplot(pjdat, aes(x=.value, y=as.factor(time), fill=comp_sp)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~comp_sp) +
  theme_test() +
  ggtitle("P. jeffreyi") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("20 year block (year +/- 10)")

ggplot(pldat, aes(x=.value, y=as.factor(time), fill=comp_sp)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~comp_sp) +
  theme_test() +
  ggtitle("P. lambertiana") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("20 year block (year +/- 10)")

ggplot(ppdat, aes(x=.value, y=as.factor(time), fill=comp_sp)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,  linetype="dashed") +
  stat_density_ridges(alpha=0) +
  facet_grid(~comp_sp) +
  theme_test() +
  ggtitle("P. ponderosa") +
  scale_fill_manual(values=cols) +
  theme(legend.position = 'none') +
  xlab("Correlation") +
  ylab("20 year block (year +/- 10)")