# Prepare data for (G)LMM analysis


rm(list=ls(all=TRUE))
#dev.off()

# libraries
library(tidyverse)
library(nlme)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
#library(cowplot)
# library(patchwork)
library(ggeffects) # to get predictions
library(sjPlot) # 'plot_model' function

# data
k.dat<-read_csv ("Data/headers_PCA.csv") 
str(k.dat) 
names (k.dat)

# the dataset combines all diversity measures
# alpha diversity measures (SR and ENSPIE) include doubled 10 m2 plots, thus "series" (i.e. 100m2 plots) should be random factor
# gamma diversity measures (SR and ENSPIE)include 100m2 or 20 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
# beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha
## SR - species richness
## ENSPIE - evenness measure calculated as inverse Simpson using species cover

# Thus we need two separate detasets for "gamma & beta" and for "alpha" (which is double sample size of that in "gamma & beta")

# 1) Entire dataset ----
## Dataset beta and gamma diversity measures for the entire dataset ----

alpha_beta_gamma<-read_csv ("Data/alpha_beta_gamma.csv") %>% 
  unite("plotID", c(series, subplot), sep="", remove=F) %>% 
  # rename(series=plot_id, subplot=subplot_id) %>% 
  filter(!(series%in%c("NFD21_44", "NFD21_46")))  # To do: do we need to remove this line?


alpha <- alpha_beta_gamma %>% 
  filter(type=="alpha")%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(k.dat, by=c("dataset", "plotID", "series", "subplot") )


# check if there are NAs after the uniting data
alpha %>% 
  filter(is.na(lat)) %>% 
  dplyr::select(plotID, series, subplot)

k.dat_mean <- k.dat %>% 
  select(c(series,eunis_group,zonality, habitat_broad, 
           where(is.numeric))) %>% 
  group_by(series, eunis_group, zonality, habitat_broad) %>% 
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  %>% 
  ungroup()



beta_gamma <- alpha_beta_gamma %>% 
  filter(type=="gamma" | type=="beta" )%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(k.dat_mean, by=c("dataset", "series") )

str(beta_gamma)


write_csv(alpha, "Data/alpha_GLM.csv")
write_csv(beta_gamma, "Data/beta_gamma_GLM.csv")


# filter header data:

library(tidyverse)
library(factoextra)

names(header_data)
header_data <- read_csv("data/headers.csv") %>% 
  select(dataset, plotID, series, subplot, 
         habitat_broad, zonality,
         BIO1, BIO12, # climate data: temperature, precipitation
         BIO7, # climate data: Temperature range
         BIO15, # climate data: precipitation CV
         grazing_intencity, mowing, 
         cover_litter,
         pH, Corg_percent)

write_csv(header_data, "data/Environm_data.csv")

read_csv("data/Environm_data.csv")
