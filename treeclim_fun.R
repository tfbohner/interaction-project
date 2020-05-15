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
  chron <- chron(rings2)
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(Site==i) %>% 
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
  
  if(i=="sp"){
    static <- cor_df
  }else{
    static <- bind_rows(static, cor_df)
  }
  
}

for(i in c("sp", "bm", "lc", "pp", "pr", "cm", "sl")){ 
  
  rings<- read.rwl(paste0("Processed Data/detrended rwl/", i, "_spl_rwl.rwl"))
  
  rings2 <- rings[which(rownames(rings)>1895),]
  chron <- chron(rings2)
  
  climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
    filter(Site==i) %>% 
    dplyr::select(c(Year, Month, spei1, tmean_C)) %>% 
    group_by(Year, Month) %>% 
    summarize_all(mean) %>% 
    as.matrix
  
  dc2 <-dcc(chrono = chron, climate = climate, method="correlation", selection = -10:9, 
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
