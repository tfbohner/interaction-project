#### dendro clim analyses
library(tidyverse)
library(dplR)
library(treeclim)

bm<- read.rwl("Processed Data/detrended rwl/bm_spl_rwl.rwl")

bm2 <- bm[which(rownames(bm)>1895),]
chron <- chron(bm2)


climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
  filter(Site=="bm") %>% 
  dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
  group_by(Year, Month) %>% 
  summarize_all(mean) %>% 
  as.matrix



corr <- dcc(chron, climate, selection = -6:9, method = "response",
    dynamic = "static")
plot(corr)

cor_df <- as.data.frame(corr$coef) %>% 
  filter(significant=="TRUE") %>% 
  mutate(abs_coef=abs(coef),
         chron=i) %>% 
  arrange(desc(abs_coef)) 



dc2 <-dcc(chrono = chron, climate = climate, method="correlation", selection = -9:8, 
          dynamic="moving", win_offset=5)
plot(dc2)

g_bm <-g_test(dc2, sb = FALSE)
g_bm

for(i in c("sp", "bm", "lc", "pp", "pr", "cm", "sl")){ 
  
  rings<- read.rwl(paste0("Processed Data/detrended rwl/", i, "_spl_rwl.rwl"))
  
  rings2 <- rings[which(rownames(rings)>1895),]
  rings2$year <- rownames(rings2)
  

  
  if(i=="sp"){
    allrings <- rings2
  }else{
    allrings <- full_join(allrings, rings2)
  }
  
}

rownames(allrings) <- allrings$year
allrings$year <- NA


## Species level----
for(s in list("AC", "PJ", "PL", "PP")) {
  spp <- allrings[,which(str_detect(colnames(allrings), s)==TRUE)]
  chrono <- chron(spp)
  
  sites <- unique(str_sub(colnames(spp), 1,2))
  
  sites <- str_replace(sites, "4.", "4")
  sites <- str_replace(sites, "8.", "8")
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(siteno %in% sites) %>% 
    mutate(spei1=ifelse(is.infinite(spei1), NA, spei1)) %>% 
    dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
    group_by(Year, Month) %>% 
    summarize_all(mean, na.rm=T) %>% 
    as.matrix
  
  resp <- dcc(chrono, climate, selection = -6:9, method = "response",
              dynamic = "static")
  resp1 <- resp$coef %>% 
    mutate(model="response",
           Species=s)
  
  corr <- dcc(chrono, climate, selection = -6:9, method = "correlation",
              dynamic = "static")
  
  coefs <- corr$coef %>% 
    mutate(model="correlation",
           Species=s) %>% 
    bind_rows(resp1)
  
  if(s=="AC") {
    all_coefs <- coefs
  } else { all_coefs <- bind_rows(all_coefs, coefs)}
}


all_coefs <- all_coefs %>% 
  mutate(monthno=ifelse(varname=="spei1", id, id-16))

write.csv(all_coefs, "Processed Data/treeclim_spp.csv", row.names = F)

ggplot(filter(all_coefs, model=="correlation"), aes(reorder(month, monthno), coef,  alpha=significant)) +
  geom_linerange(stat="identity", aes(ymin=0, ymax=coef, color=Species), size=3, position = position_dodge(width=0.5)) +
  geom_point(data=filter(all_coefs, model=="response"),
             aes(reorder(month, monthno), coef, fill=Species),
             color="black", shape=21, size=2, position = position_dodge(width=0.5)) +
  theme_test() +
  facet_wrap(~varname, ncol = 1, scales='free') +
  geom_hline(yintercept = 0, size=0.2)

  


plot(resp$coef$coef, corr$coef$coef)
abline(0,1)

for(s in list("AC", "PJ", "PL", "PP")) {
  spp <- allrings[,which(str_detect(colnames(allrings), s)==TRUE)]
  chrono <- chron(spp)
  
  sites <- unique(str_sub(colnames(spp), 1,2))
  
  sites <- str_replace(sites, "4.", "4")
  sites <- str_replace(sites, "8.", "8")
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(siteno %in% sites) %>% 
    mutate(spei1=ifelse(is.infinite(spei1), NA, spei1)) %>% 
    dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
    group_by(Year, Month) %>% 
    summarize_all(mean, na.rm=T) %>% 
    as.matrix
  
  dc2 <-dcc(chrono = chrono, climate = climate, method="correlation", selection = -10:9, 
            dynamic="moving")
  plot(dc2)
  
  g_bm <-g_test(dc2,  sb = FALSE)
  g_bm <- g_bm %>% 
    filter(p<0.05) %>% 
    mutate(Species=s)
  
  if(s=="AC"){
    dyn_df <- g_bm
  }else{
    dyn_df <- bind_rows(dyn_df, g_bm)
  }
  
}



## Site level----
for(i in c("sp", "bm", "lc", "pp", "pr", "cm", "sl")){ 
  
  rings<- read.rwl(paste0("Processed Data/detrended rwl/", i, "_spl_rwl.rwl"))
  
  rings2 <- rings[which(rownames(rings)>1895),]
  chrono <- chron(rings2)
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(Site==i) %>% 
    dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
    group_by(Year, Month) %>% 
    summarize_all(mean) %>% 
    as.matrix
  
  corr <- dcc(chrono, climate, selection = -6:9, method = "response",
              dynamic = "static")
  plot(corr)
  
  cor_df <- as.data.frame(corr$coef) %>% 
    filter(significant=="TRUE") %>% 
    mutate(abs_coef=abs(coef),
           chron=i) %>% 
    arrange(desc(abs_coef)) 
  
  if(i=="sp"){
    static <- cor_df
  }else{
    static <- bind_rows(static, cor_df)
  }
  
}

for(i in c("sp", "bm", "lc", "pp", "pr", "cm", "sl")){ 
  
  rings<- read.rwl(paste0("Processed Data/detrended rwl/", i, "_spl_rwl.rwl"))
  
  rings2 <- rings[which(rownames(rings)>1895),]
  chrono <- chron(rings2)
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(Site==i) %>% 
    dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
    group_by(Year, Month) %>% 
    summarize_all(mean) %>% 
    as.matrix
  
  dc2 <-dcc(chrono = chrono, climate = climate, method="correlation", selection = -10:9, 
            dynamic="moving", win_offset=5)
  plot(dc2)
  
  g_bm <-g_test(dc2,  sb = FALSE)
  g_bm <- g_bm %>% 
    filter(p<0.05) %>% 
    mutate(chron=i)
  
  if(i=="sp"){
    dyn_df <- g_bm
  }else{
    dyn_df <- bind_rows(dyn_df, g_bm)
  }
  
}
