# Purpose: Analysis for among-subplot dissimilarity in species composition
# and spatial aggregation

# libraries----
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(performance)
library(patchwork)


# prepare data----

# "data/climate_PCA.csv" contains scores for the compound climate variable, 
# derived from the PCA analysis in "1_prepare_data/ PCA_environment.R"

# "data/Environm_variabl.csv" contains all environmental data

bray_curtis<- read_csv("data/aggregation.csv") 

climate_PCA <- read.csv("data/climate_PCA.csv")

header <- read_csv("data/Environm_variabl.csv") %>% 
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header) 
names (header)


header_mean <- header %>% 
  select(c(series, lat, lon, zonality, habitat_broad, 
           where(is.numeric))) %>% 
  group_by(series,  zonality, habitat_broad) %>% 
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  %>% 
  ungroup()


beta_gamma <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>% 
  filter(type=="gamma" | type=="beta" )%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header_mean, by=c("dataset", "series") )%>% 
  mutate(dataset=factor(dataset))

str(beta_gamma) 
names (beta_gamma)

# dataset is a separate vegetation survey campaign
beta_gamma$dataset 

## data beta----
beta_data <- beta_gamma %>% 
  dplyr::select(dataset, series, habitat_broad, zonality,
                gamma_100_div, gamma_100_ENSPIE, gamma_100_cover,
                beta_100_div, beta_100_ENSPIE, 
                lat, lon, pca1_clima, 
                grazing_intencity, mowing, 
                # cover_shrub_total,     inclination, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                ) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na %>%
  left_join(aggr, by = join_by(series)) %>%
  left_join(bray_curtis, by = join_by(series)) %>%
  mutate(aggreg=mean_corner_diff*sd_corner_diff,
         BRAY_tn = beta.BRAY.BAL, 
         BRAY_nest = beta.BRAY.GRA, 
         BRAY_dissiml = beta.BRAY)

str(beta_data)
summary(beta_data)


beta_data  %>% 
  select(gamma_100_div, gamma_100_ENSPIE, gamma_100_cover,,
         beta_100_div, beta_100_ENSPIE,
         mean_corner_diff, # mean value of between-subplot differences in cover among species 
         sd_corner_diff,   # sd value of between-subplot differences in cover among species
         aggreg,
         beta.BRAY.BAL, beta.BRAY.GRA, beta.BRAY ) %>% 
 # rename(SR=sowndiv, "FG richness"=numfg) %>% 
  cor() |>
  ggcorrplot::ggcorrplot(
    lab = TRUE, type = "lower",
    outline.color = "white",
    colors = c("red", "white", "blue")
  ) 








### Simulate data to check indices --------

sp1 <-c(2, 14)
sp2 <-c(15, 1)
sp3 <-c(0, 2)

sim_dat <- tibble(sp1, sp2, sp3) 
sim_dat

betapart::beta.pair(sim_dat %>% 
                      mutate_at(1:3, ~ replace(., . > 0, 1)),
                    index.family = "jaccard")|>
  as_tibble() %>% 
 rename(Jaccard_tn = beta.jtu,
        Jaccard_nest = beta.jne,
        Jaccard_dissiml = beta.jac)

betapart::beta.multi.abund(sim_dat, index.family = "bray")|>
  as_tibble() %>% 
  rename(BRAY_tn = beta.BRAY.BAL,  # balanced variation in abundance, 
# whereby the individuals of some species in one site 
# are substituted by the same number of individuals of different species in another site
BRAY_nest = beta.BRAY.GRA,  # some individuals are lost from one site to the other
BRAY_dissiml = beta.BRAY)

sim_dat


# 2) 
sp1 <-c(1, 13)
sp2 <-c(13, 1)
sp3 <-c(0, 14)
sp4 <-c(14, 0)

sim_dat <- tibble(sp1, sp2, sp3, sp4) 
sim_dat

betapart::beta.pair(sim_dat %>% 
                      mutate_at(1:4, ~ replace(., . > 0, 1)),
                    index.family = "jaccard")|>
  as_tibble() %>% 
  rename(Jaccard_tn = beta.jtu,
         Jaccard_nest = beta.jne,
         Jaccard_dissiml = beta.jac)

betapart::beta.multi.abund(sim_dat, index.family = "bray")|>
  as_tibble() %>% 
  rename(BRAY_tn = beta.BRAY.BAL,  # balanced variation in abundance, 
         # whereby the individuals of some species in one site 
         # are substituted by the same number of individuals of different species in another site
         BRAY_nest = beta.BRAY.GRA,  # some individuals are lost from one site to the other
         BRAY_dissiml = beta.BRAY)

sim_dat


# 3) 
sp1 <-c(1, 13)
sp2 <-c(13, 1)
sp3 <-c(0, 15)
sp4 <-c(15, 0)
sp5 <-c(0, 1)
sim_dat <- tibble(sp1, sp2, sp3, sp4, sp5) 
sim_dat

betapart::beta.pair(sim_dat %>% 
                      mutate_at(1:5, ~ replace(., . > 0, 1)),
                    index.family = "jaccard")|>
  as_tibble() %>% 
  rename(Jaccard_tn = beta.jtu,
         Jaccard_nest = beta.jne,
         Jaccard_dissiml = beta.jac)

betapart::beta.multi.abund(sim_dat, index.family = "bray")|>
  as_tibble() %>% 
  rename(BRAY_tn = beta.BRAY.BAL,  # balanced variation in abundance, 
         # whereby the individuals of some species in one site 
         # are substituted by the same number of individuals of different species in another site
         BRAY_nest = beta.BRAY.GRA,  # some individuals are lost from one site to the other
         BRAY_dissiml = beta.BRAY)

sim_dat




# (1) aggregation -------

# GLLM----

## Exploration----

# Check the model:
m <- lmer(BRAY_tn ~ 
                   poly(pca1_clima, 2) +
                   poly(pH, 2) +
                   poly(Corg_percent,2)+
                   poly(cover_litter,2) +
                   grazing_intencity + mowing +
                   (1|dataset),  data = beta_data)

# check model
plot(m) # heteroscedasticity
qqnorm(resid(m))
qqline(resid(m))



Anova(m)

# Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 
# Model selection ----

# poly(pca1_clima, 2) 
# poly(pH, 2) is marginal

m1_1 <- lmer(BRAY_tn ~  poly(pca1_clima, 2) +
                      pH +
                     Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)

m1_2 <- lmer(BRAY_tn ~ pca1_clima +
               pH +
               Corg_percent +
               poly(cover_litter,2) +
               grazing_intencity + mowing +
               (1|dataset),  data = beta_data)


m1_3 <- lmer(BRAY_tn ~  poly(pca1_clima, 2) +
                      pH +
                      Corg_percent +
                      cover_litter +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)


m1_4 <- lmer(BRAY_tn ~ 
               pca1_clima +
               pH +
               Corg_percent +
               cover_litter +
               grazing_intencity + mowing +
               (1|dataset),  data = beta_data)


# calculate and compare AIC 
AIC(m1_1 , m1_2 , m1_3 , m1_4 ) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m1_1 )
# Anova(m1_4_ENSPIE)


# Model 2: Add Prec_Varieb ----

### Model selection -----
m2_1 <- lmer(BRAY_tn ~ 
               poly(pca1_clima, 2) +
               poly(Prec_Varieb,2) +
               pH +
               Corg_percent +
               poly(cover_litter,2) +
               grazing_intencity + mowing +
               (1|dataset),  data = beta_data)


m2_2 <- lmer(BRAY_tn ~ 
               poly(pca1_clima, 2) +
               Prec_Varieb +
               pH +
               Corg_percent +
               poly(cover_litter,2) +
               grazing_intencity + mowing +
               (1|dataset),  data = beta_data)


# calculate and compare AIC 
AIC(m2_1, m2_2) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))


Anova(m2_1)



## -> Additional analysis (requested by reviewer) -----
# Comment "I am unsure if you can conclude that the precipitation variability 
# effect is shown. Perhaps you just repeat the hump-shape pattern of the environmental gradient. 
# Maybe you can compare if models, where you have both linear and quadratic terms of climate, 
# differ from those where the quadratic term of climate is replaced by precipitation variability. 
# If the latter model is better (e.g. AIC is smaller than two units), you have more support for this claim."

m3_1  <- lmer(BRAY_tn ~ 
                      poly(pca1_clima, 2) +
                      (1|dataset),  data = beta_data)

m3_2  <- lmer(BRAY_tn ~ 
                      pca1_clima +
                      Prec_Varieb +
                      (1|dataset),  data = beta_data)

# calculate and compare AIC 
AIC(m3_1 , m3_2 ) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m3_1 )
Anova(m3_2 )

# Plots----
#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")

set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,legend.pos = "None", geom.linetype = 2)

# final models:
Anova(m1_1)
Anova(m2_1)


## Clima ----

clima_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_beta_aggr

Fig.betaaggr_clima <- ggplot(clima_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(pca1_clima, BRAY_tn, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Climate gradient')+ 
  geom_line(linetype=1, size=1, col="black") 

Fig.betaaggr_clima 


plot_model(m1_1_aggr,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]", #  show.data=F
           title = "", line.size=0.5) + aes(linetype="solid") +
  labs(y=expression(paste("beta ENS"[PIE])), x='Climate gradient (PC)')+ 
  geom_point(data = beta_data, aes(pca1_clima, BRAY_tn, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 



## Soil C ----
Humus_pred_beta_aggr <- get_model_data(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_beta_aggr

Fig.betaaggr_soilC <- ggplot(Humus_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Corg_percent, BRAY_tn, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Soil C')+ 
  geom_line(linetype=5, size=0.5, col="black") 

Fig.betaaggr_soilC 


## Litter % ----
Litter_pred_beta_aggr <- get_model_data(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_beta_aggr

Fig.betaaggr_Litter <- ggplot(Litter_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(cover_litter, BRAY_tn, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Litter cover')+ 
  geom_line(linetype=1, size=1, col="black") 

Fig.betaaggr_Litter 

## Soil pH ----
pH_pred_beta_aggr <- get_model_data(m1_1,type = "pred", terms="pH[3.7:9, by=.001]")
pH_pred_beta_aggr

Fig.betaaggr_soil.pH <- ggplot(pH_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(pH, BRAY_tn, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Soil pH')+ 
  geom_line(linetype=5, size=0.5, col="black") 


Fig.betaaggr_soil.pH 


## Grazing ----
grazing_pred_beta_aggr <- get_model_data(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_beta_aggr

Fig.betaaggr_grazing <- ggplot(grazing_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(grazing_intencity, BRAY_tn, 
                                 fill = habitat, col=habitat), 
             size=1, alpha=0.8, pch=21,
             position=position_jitter(w=0.2))+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Grazing intencity')+ 
  geom_line(linetype=5, size=0.5, col="black") 

Fig.betaaggr_grazing


## Mowing ----

plot_model(m1_1,type = "pred", terms="mowing",  show.data=T,
           title = "", line.size=0.5)

Fig.betaaggr_mowing <-ggplot(beta_data, aes(mowing, BRAY_tn)) + 
  geom_boxplot(color="black")+
  labs(y="Species aggregation", x='Mowing')+ 
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) +
  labs(color='Habitat type')

Fig.betaaggr_mowing 


## Prec_Varieb ----
precipCV_pred_beta_aggr <- get_model_data(m2_1,type = "pred", terms="Prec_Varieb")
precipCV_pred_beta_aggr

Fig.betaaggr_precip.CV <- ggplot(precipCV_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Prec_Varieb, BRAY_tn, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species aggregation", x='Precipitation CV')+ 
  geom_line(linetype=5, size=0.5, col="black") 

Fig.betaaggr_precip.CV

## Merge plots----

library(patchwork)

set_theme(base = theme_bw(), axis.textsize.x = 0.8, axis.textsize.y = 0.8, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 0.9,legend.pos = "None", 
          geom.linetype = 2)

Fig.betaaggr_clima + 
  Fig.betaaggr_precip.CV +
  Fig.betaaggr_soilC + 
  Fig.betaaggr_soil.pH + 
  Fig.betaaggr_Litter + 
  Fig.betaaggr_grazing +
  plot_annotation(tag_levels = 'a') +
  plot_layout(ncol=3)  & #ylab(NULL) & 
  theme( #plot.margin = margin(3, 1, 3, 20), 
                                           plot.tag = element_text(size = 11, face = 'bold'), 
                                           plot.tag.position = c(0.25, 1.06))
