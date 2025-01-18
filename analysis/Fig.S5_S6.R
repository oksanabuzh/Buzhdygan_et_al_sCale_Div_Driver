# Purpose: produce Fig S5, Fig.S6

# libraries----
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(performance)
library(patchwork)


# data----
# series are 100m2 plots
# subplot is one of the two corners (two 10 m2 plots) within same series
# dataset is sampling campaign

## climate PCA data----
climate_PCA <- read_csv("data/climate_PCA.csv")

## environmental variables data----
header <- read_csv("data/Environm_variabl.csv") %>% 
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header) 
names (header)

## data for 100m2 (gamma) ----

# mean per series (per 100m2 plots)
header_mean <- header %>% 
  select(c(series,zonality, habitat_broad, 
           where(is.numeric))) %>% 
  group_by(series, zonality, habitat_broad) %>% 
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  %>% 
  ungroup()

# join with dataset for diversity measures
beta_gamma <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>% 
  filter(type=="gamma" | type=="beta" )%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header_mean, by=c("dataset", "series") )

str(beta_gamma) 
names (beta_gamma)


beta_gamma$dataset <- factor(beta_gamma$dataset)# dataset is a separate vegetation survey campaign


# selected variables
gamma_data <- beta_gamma %>% 
  dplyr::select(gamma_100_div, gamma_100_ENSPIE, 
                pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,BIO1, BIO12,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         Temprt = BIO1,
         Precipt=BIO12,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, 
                              c("saline", "complex", "dry", 
                                "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na



str(gamma_data)


## data for 10 m2 (alpha) ----
alpha <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type=="alpha")%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header, 
            by=c("dataset", "plotID", "series", "subplot")
  )
str(alpha) 
names (alpha)


alpha$dataset <- factor(alpha$dataset)# dataset is a separate vegetation survey campaign

# selected variables 
alpha_data <- alpha %>% 
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,  
                pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad,
                subplot) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, 
                              c("saline", "complex", "dry",
                                "wet" , "mesic", "fringe", "alpine"))) %>% 
  
  drop_na

# Take mean across two subplots for the purpose of plotting:
# plot on a mean alpha per series to omit pseudo-replication of the plots on the figure:
alpha_mean <- alpha_data %>% 
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,  
                pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

str(alpha_mean)


## data for beta ----

beta_data <- beta_gamma %>% 
  dplyr::select(beta_100_div, beta_100_ENSPIE, 
                pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, 
                              c("saline", "complex", "dry",
                                "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

str(beta_data)

# 1) Get the predictions for the final selected models ----  


# alpha SR ----
m1_3 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

m1_1 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Corg_percent,2)+
                 poly(pH,2) +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

m2_1 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Prec_Varieb, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

## Predictions----
clima_pred_10m <-get_model_data(m1_3,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m <-get_model_data(m1_3,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_10m <-get_model_data(m1_3,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_10m <-get_model_data(m1_1,type = "pred", terms="pH[3.8:9, by=.001]")
grazing_pred_10m <-get_model_data(m1_3,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m <-get_model_data(m2_1,type = "pred", terms="Prec_Varieb")


# alpha ENSPIE----
m1_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

m2_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

## Predictions----

clima_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="pH[3.8:9, by=.001]")
grazing_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m_Ensp <-get_model_data(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb")

# gamma SR ----
m1_1 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

m2_1 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Prec_Varieb, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

## Predictions----

clima_pred_100m <- get_model_data(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m <- get_model_data(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_100m <- get_model_data(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_100m <- get_model_data(m1_1,type = "pred", terms="pH[3.7:9, by=.001]")
grazing_pred_100m <- get_model_data(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m <- get_model_data(m2_1,type = "pred", terms="Prec_Varieb")

# gamma ENSPIE ----

m1_3_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)

m2_1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)

## Predictions----

clima_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE,type = "pred", terms="pH[3.7:9, by=.001]")
grazing_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m_Ensp <- get_model_data(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb")


# beta SR ----

m1_1 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)


m2_1 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Prec_Varieb, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)

## Predictions----

clima_pred_beta <- get_model_data(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_beta <- get_model_data(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_beta <- get_model_data(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_beta <- get_model_data(m1_1,type = "pred", terms="pH[3.7:9, by=.001]")
grazing_pred_beta <- get_model_data(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_beta <- get_model_data(m2_1,type = "pred", terms="Prec_Varieb")

# beta ENSPIE ----

m1_1_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(pH, 2) +
                      poly(Corg_percent,2) +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)
m2_1_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)

## Predictions----

clima_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Litter_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]")
pH_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="pH[3.7:9, by=.001]")
grazing_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_beta_ENSPIE <- get_model_data(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb")


# 2) Plot the predictions  ----  

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


# alpha-gamma----

## SR----

### Clima ----

Fig_Clima_10  <- ggplot(clima_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(pca1_clima, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Climate gradient')+
  geom_line(linetype=1, size=1, col="#50A0C8") 



Fig_Clima_100  <- ggplot(clima_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=clima_pred_100m, aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(pca1_clima, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Climate gradient')+
  geom_line(data=clima_pred_100m, linetype=1, size=1, col="#D6604D") 

### Soil C ----
Fig_Humus_10  <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(Corg_percent, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Soil C')+
  geom_line(linetype=1, size=1, col="#50A0C8")


Fig_Humus_100  <- ggplot(Humus_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=Humus_pred_100m  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(Corg_percent, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Soil C')+
  geom_line(data=Humus_pred_100m, linetype=1, size=1, col="#D6604D") 


### Litter % ----

Fig_Litter_10  <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(cover_litter, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Litter cover')+
  geom_line(linetype=1, size=1, col="#50A0C8")

Fig_Litter_100  <- ggplot(Litter_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=Litter_pred_100m  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(cover_litter, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Litter cover')+
  geom_line(data=Litter_pred_100m, linetype=1, size=1, col="#D6604D") 

### Soil pH ----

Fig_pH_10  <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
    geom_point(data=alpha_mean, 
             aes(pH, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Soil pH')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") 

Fig_pH_100  <- ggplot(pH_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=pH_pred_100m  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(pH, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Soil pH')+
  geom_line(data=pH_pred_100m, linetype=1, size=0.5, col="#D6604D") 


### Grazing ----

Fig_grazing_10  <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(grazing_intencity, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21, position=position_jitter(w=0.2))+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") 


Fig_grazing_100  <- ggplot(grazing_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=grazing_pred_100m  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(grazing_intencity, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21, position=position_jitter(w=0.2))+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Grazing intencity')+
  geom_line(data=grazing_pred_100m, linetype=5, size=0.5, col="#D6604D") 


### Prec_Varieb ----

Fig_precipCV_10  <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(Prec_Varieb, alpha_10_div, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Precipitation CV')+
  geom_line(linetype=1, size=1, col="#50A0C8") 

Fig_precipCV_100  <- ggplot(precipCV_pred_100m, aes(x, predicted)) +
  geom_ribbon(data=precipCV_pred_100m, aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(Prec_Varieb, gamma_100_div, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Precipitation CV')+
  geom_line(data=precipCV_pred_100m, linetype=1, size=1, col="#D6604D") 


## ENSPIE ----

### Clima----

Fig_Clima_10_ENSPIE  <- ggplot(clima_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(pca1_clima, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Climate gradient')+ 
  geom_line(linetype=1, size=1, col="#50A0C8") 


Fig_Clima_100_ENSPIE  <- ggplot(clima_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=clima_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(pca1_clima, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
 scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Climate gradient')+ 
  geom_line(data=clima_pred_100m_Ensp, linetype=1, size=1, col="#D6604D") 



### Soil C ----
Fig_Humus_10_ENSPIE  <- ggplot(Humus_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(Corg_percent, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil C')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") 

Fig_Humus_100_ENSPIE  <- ggplot(Humus_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=Humus_pred_100m_Ensp  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(Corg_percent, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil C')+
  geom_line(data=Humus_pred_100m_Ensp, linetype=5, size=0.5, col="#D6604D") 

### Litter % ----

Fig_Litter_10_ENSPIE  <- ggplot(Litter_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(cover_litter, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Litter cover')+
  geom_line(linetype=1, size=0.5, col="#50A0C8") 


Fig_Litter_100_ENSPIE  <- ggplot(Litter_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=Litter_pred_100m_Ensp  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(cover_litter, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
   scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Litter cover')+
  geom_line(data=Litter_pred_100m_Ensp, linetype=1, size=1, col="#D6604D") 



### Soil pH ----

Fig_pH_10_ENSPIE  <- ggplot(pH_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(pH, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil pH')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") 

Fig_pH_100_ENSPIE  <- ggplot(pH_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=pH_pred_100m_Ensp  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
 geom_point(data=gamma_data, 
             aes(pH, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil pH')+
  geom_line(data=pH_pred_100m_Ensp, linetype=5, size=0.5, col="#D6604D") 



### Grazing ----

Fig_grazing_10_ENSPIE  <- ggplot(grazing_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(grazing_intencity, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21, position=position_jitter(w=0.2))+
  scale_fill_manual(values = col) +  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") 

Fig_grazing_100_ENSPIE  <- ggplot(grazing_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=grazing_pred_100m_Ensp  , aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(grazing_intencity, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21, position=position_jitter(w=0.2))+
  scale_fill_manual(values = col) +  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Grazing intencity')+
  geom_line(data=grazing_pred_100m_Ensp, linetype=1, size=0.5, col="#D6604D") 



### Prec_Varieb ----

Fig_precipCV_10_ENSPIE  <- ggplot(precipCV_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4,  fill="#A6CEE3")+
  geom_point(data=alpha_mean, 
             aes(Prec_Varieb, alpha_10_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8,
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Precipitation CV')+
  geom_line(linetype=1, size=1, col="#50A0C8") 



Fig_precipCV_100_ENSPIE  <- ggplot(precipCV_pred_100m_Ensp, aes(x, predicted)) +
  geom_ribbon(data=precipCV_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.1, fill="#D6604D")+
  geom_point(data=gamma_data, 
             aes(Prec_Varieb, gamma_100_ENSPIE, fill = habitat, col = habitat), size=1, alpha=0.8, 
             pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y=expression(paste("ENS"[PIE])), x='Precipitation CV')+
  geom_line(data=precipCV_pred_100m_Ensp, linetype=1, size=1, col="#D6604D") 


## Merge plots----

library(patchwork)

set_theme(base = theme_bw(), axis.textsize.x = 0.6, axis.textsize.y = 0.6, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 0.65,legend.pos = "None", 
          geom.linetype = 2)
FigS5 <- 
  Fig_Clima_10 +  Fig_Clima_100 + 
      Fig_Clima_10_ENSPIE +  Fig_Clima_100_ENSPIE +   
  Fig_precipCV_10 +   Fig_precipCV_100 + 
     Fig_precipCV_10_ENSPIE + Fig_precipCV_100_ENSPIE +   
  Fig_Humus_10 +  Fig_Humus_100 + 
     Fig_Humus_10_ENSPIE +  Fig_Humus_100_ENSPIE + 
  Fig_pH_10 +  Fig_pH_100 + 
     Fig_pH_10_ENSPIE + Fig_pH_100_ENSPIE +   
  Fig_Litter_10 +  Fig_Litter_100 +  
     Fig_Litter_10_ENSPIE +  Fig_Litter_100_ENSPIE +  
  Fig_grazing_10 + Fig_grazing_100 + 
    Fig_grazing_10_ENSPIE +  Fig_grazing_100_ENSPIE +    
  #
  #
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol=4) & ylab(NULL) & theme(plot.margin = margin(3, 1, 3, 20), 
                                           plot.tag = element_text(size = 6, face = 'bold'), 
                                           plot.tag.position = c(0.15, 1.06))

FigS5



## Fig S6 

Fig.alphaSR_mowing <-ggplot(alpha_mean, aes(mowing, alpha_10_div)) + 
  geom_boxplot(color="#64ABCE")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
  labs(y="Species richness", x='Mowing')

Fig.alphaSR_mowing 


Fig.alphaENSPIE_mowing <-ggplot(alpha_mean, aes(mowing, alpha_10_ENSPIE)) + 
  geom_boxplot(color="#64ABCE")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
  labs(y=expression(paste("ENS"[PIE])), x='Mowing')+ labs(color='Habitat type')

Fig.alphaENSPIE_mowing 

Fig.gammaSR_mowing <-ggplot(gamma_data, aes(mowing, gamma_100_div)) + 
  geom_boxplot(color="#D6604D")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
  labs(y="Species richness", x='Mowing')

Fig.gammaSR_mowing 
 

Fig.gammaENSPIE_mowing <-ggplot(gamma_data, aes(mowing, gamma_100_ENSPIE)) + 
  geom_boxplot(color="#D6604D")+
  labs(y=expression(paste("gamma ENS"[PIE])), x='Mowing')+ 
  geom_point(aes(color = habitat), position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+
  labs(color='Habitat type')

Fig.gammaENSPIE_mowing 


Fig.betaSR_mowing <-ggplot(beta_data, aes(mowing, beta_100_div)) + 
  geom_boxplot(color="#00AC7F")+
  labs(y="beta SR", x='Mowing')+ 
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) 

Fig.betaSR_mowing 

Fig.betaENSPIE_mowing <-ggplot(beta_data, aes(mowing, beta_100_ENSPIE)) + 
  geom_boxplot(color="#00AC7F")+
  labs(y=expression(paste("ENS"[PIE])), x='Mowing')+ 
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=1, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) +
  labs(color='Habitat type')

Fig.betaENSPIE_mowing 



Fig.alphaSR_mowing + Fig.gammaSR_mowing + Fig.betaSR_mowing +
Fig.alphaENSPIE_mowing + Fig.gammaENSPIE_mowing + Fig.betaENSPIE_mowing +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol=3) & ylab(NULL) & theme(plot.margin = margin(3, 1, 3, 20), 
                                           plot.tag = element_text(size = 6, face = 'bold'), 
                                           plot.tag.position = c(0.15, 1.06))






