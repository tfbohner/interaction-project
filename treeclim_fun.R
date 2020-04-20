#### dendro clim analyses
library(tidyverse)
library(dplR)
library(treeclim)

bm <- read.rwl("Raw Data/RWL copies/Black_Mountain.rwl")
names(bm)[str_detect(names(bm), "11PL1_3A")==T] <- "11PL1_1C"
names(bm)[str_detect(names(bm), "11PL1_3B")==T] <- "11PL1_1D"

bm2 <- bm[which(rownames(bm)>1895),]
bm3 <- detrend(bm2, method="Spline")
bm_chron <- chron(bm3)

climate <- read.csv("Raw Data/monthly_precip_temp_spei.csv") %>% 
  filter(Site=="bm") %>% 
  dplyr::select(c(Year, Month, ppt_mm, tmean_C)) %>% 
  group_by(Year, Month) %>% 
  summarize_all(mean) %>% 
  as.matrix


corr <- dcc(bm_chron, climate, selection = -6:9, method = "response",
    dynamic = "static", win_size = 25)
plot(corr)

dc2 <-dcc(chrono = bm_chron, climate = climate, method="correlation", selection = -9:8, 
          dynamic="moving", win_offset=5)
plot(dc2)

g_test(dc2, sb = FALSE)
