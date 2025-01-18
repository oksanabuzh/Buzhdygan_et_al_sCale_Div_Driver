# Purpose: GLMM analysis for beta diversity


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

# "data/alpha_beta_gamma_community_variabl.csv" combines all diversity measures and plant cover
# alpha diversity measures (SR and ENSPIE) include doubled 10 m2 plots, 
## thus "series" (i.e. 100m2 plots), nested in dataset (separate vegetation survey campaign) 
## are fitted as a random effect
# gamma diversity measures (SR and ENSPIE)include 100m2 plots (i.e. the sample size is half of what we have for the 10m2 plots)
# beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

## SR - species richness
## ENSPIE - evenness measure calculated as inverse Simpson using species cover
## cover - is cumulative plant cover


# "data/climate_PCA.csv" contains scores for the compound climate variable, 
# derived from the PCA analysis in "1_prepare_data/ PCA_environment.R"

# "data/Environm_variabl.csv" contains all environmental data


climate_PCA <- read.csv("data/climate_PCA.csv")

header <- read_csv("data/Environm_variabl.csv") %>% 
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header) 
names (header)


header_mean <- header %>% 
  select(c(series, zonality, habitat_broad, 
           where(is.numeric))) %>% 
  group_by(series,  zonality, habitat_broad) %>% 
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  %>% 
  ungroup()


# prepare subset of data for beta-diversity analysis
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

# Select variables, remove NAs 
beta_data <- beta_gamma %>% 
  dplyr::select(beta_100_div, beta_100_ENSPIE, 
                pca1_clima, 
                grazing_intencity, mowing, 
                # cover_shrub_total,     inclination, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na





# (1) beta SR -----


# Data Exploration----

str(beta_data)

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(beta_data, aes(pca1_clima, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(Prec_Varieb, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(cover_litter, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(pH, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(Corg_percent, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(grazing_intencity, beta_100_div)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(mowing, beta_100_div)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
 scale_color_manual(values = col)+
  labs(color='Habitat type')



# GLLM----

## Exploration----

m <- lmer (beta_100_div ~ 
              poly(pca1_clima, 2) +
              poly(Corg_percent,2)+
              poly(pH, 2) +
              poly(cover_litter,2) +
              grazing_intencity + mowing +
              (1|dataset),  data = beta_data)


check_convergence(m)
  
# check model
plot(m)
qqnorm(resid(m))
qqline(resid(m))

# check multicolinearity
check_collinearity(m)
  

 
Anova(m)
summary(m)

# Model selection ----

## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 

# poly(pca1_clima, 2) nonsignificant
# poly(Corg_percent,2)nonsignificant

m1_1 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)

m1_2 <- lmer (beta_100_div ~ 
                pca1_clima +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)

m1_3 <- lmer (beta_100_div ~ 
                pca1_clima +
                Corg_percent +
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)


m1_4 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                Corg_percent+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)

# calculate and compare AIC 
AIC(m1_1, m1_2, m1_3, m1_4) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m1_1)

## Model 2: Add Prec_Varieb ----
m2_1 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Prec_Varieb, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)

m2_2 <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                Prec_Varieb +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)


# calculate and compare AIC 
AIC(m2_1, m2_2) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m2_1)
# Anova(m2_2)


# Plots----
#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")

set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,legend.pos = "None", geom.linetype = 2)

# final models:
Anova(m1_1)
Anova(m2_1)


## Clima ----

clima_pred_beta <- get_model_data(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_beta

Fig.betaSR_clima <- ggplot(clima_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(pca1_clima, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="beta SR", x='Climate gradient')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaSR_clima 

plot_model(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]", #  show.data=F
           title = "", line.size=0.5) + 
  labs(y="beta SR", x='Climate gradient (PC)')+ 
  geom_point(data = beta_data, aes(pca1_clima, beta_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 




## Soil C ----

Humus_pred_beta <- get_model_data(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_beta

Fig.betaSR_soilC <- ggplot(Humus_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Corg_percent, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="beta SR", x='Soil C')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaSR_soilC 

plot_model(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]",# show.data=F,
                                title = "", line.size=0.5)  +
  labs(y="beta SR", x='Soil C')+ 
  geom_point(data = beta_data, aes(Corg_percent, beta_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

 


## Litter % ----

Litter_pred_beta <- get_model_data(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_beta

Fig.betaSR_Litter <- ggplot(Litter_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(cover_litter, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="beta SR", x='Litter cover')+ 
  geom_line(linetype=1, size=1, col="#00AC7F") 

Fig.betaSR_Litter 

 plot_model(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]",# show.data=F,
                                  title = "", line.size=1) + aes(linetype="solid") +
  labs(y="beta SR", x='Litter cover')+ 
  geom_point(data = beta_data, aes(cover_litter, beta_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 




## Soil pH ----
pH_pred_beta <- get_model_data(m1_1,type = "pred", terms="pH[3.7:9, by=.001]")
pH_pred_beta
 
Fig.betaSR_soil.pH <- ggplot(pH_pred_beta, aes(x, predicted)) +
   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
   geom_point(data=beta_data, aes(pH, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
   scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
   labs(y="beta SR", x='Soil pH')+ 
   geom_line(linetype=1, size=1, col="#00AC7F") 
 
Fig.betaSR_soil.pH 
 
plot_model(m1_1,type = "pred", terms="pH[3.7:9, by=.001]",# show.data=F,
                                title = "", line.size=1) +  aes(linetype="solid") +
  labs(y="beta SR", x='Soil pH')+ 
  geom_point(data = beta_data, aes(pH, beta_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

 

## Grazing ----
grazing_pred_beta <- get_model_data(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_beta

Fig.betaSR_grazing <- ggplot(grazing_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(grazing_intencity, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="beta SR", x='Grazing intencity')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaSR_grazing 


plot_model(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]",# show.data=F,
                                  title = "", line.size=0.5)  +  
  labs(y="beta SR", x='Grazing intencity')+ 
  geom_point(data = beta_data, aes(grazing_intencity, beta_100_div, color = habitat), size=3, alpha=0.8,
               position=position_jitter(w=0.2))+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

Fig.betaSR_grazing


## Mowing ----

plot_model(m1_1,type = "pred", terms="mowing",  show.data=T,
                                  title = "", line.size=0.5)

Fig.betaSR_mowing <-ggplot(beta_data, aes(mowing, beta_100_div)) + 
  geom_boxplot(color="#00AC7F")+
  labs(y="beta SR", x='Mowing')+ 
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) 

Fig.betaSR_mowing 



## Prec_Varieb ----

precipCV_pred_beta <- get_model_data(m2_1,type = "pred", terms="Prec_Varieb")
precipCV_pred_beta

Fig.betaSR_precip.CV <- ggplot(precipCV_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Prec_Varieb, beta_100_div, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="beta SR", x='Precipitation variability')+ 
  geom_line(linetype=1, size=1, col="#00AC7F") 

Fig.betaSR_precip.CV

plot_model(m2_1,type = "pred", terms="Prec_Varieb", #  show.data=F
                                title = "", line.size=1) + aes(linetype="solid") +
  labs(y="beta SR", x='Precipitation variability')+ 
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

 



# Plots combined ----

# for Gamma SR (100 m^2)

Fig.betaSR_clima + 
  Fig.betaSR_precip.CV +
  Fig.betaSR_soilC + 
  Fig.betaSR_soil.pH + 
  Fig.betaSR_Litter + 
  Fig.betaSR_grazing +
  plot_layout(ncol=3)


# final statistics----
Anova(m1_1)
Anova(m2_1) # for Prec_Varieb

# write.csv(Anova(m1_1),  file = "results/glmer_beta_SR.csv")
# write.csv(coef(summary(m1_1)),  file = "results/summary_beta_SR.csv")
# write.csv(Anova(m2_1),  file = "results/glmer_beta_SR_2.csv")
# write.csv(coef(summary(m2_1)),  file = "results/summary_beta_SR_2.csv")

# R2 for the entire model---------
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination

# write.csv(MuMIn::r.squaredGLMM(m1_1),  file = "results/Mod1_R2_beta_SR.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1),  file = "results/Mod2_R2_beta_SR.csv")

# Partial R2 for fixed effects
R1_SR <- r2glmm::r2beta(m1_1,  partial = T, data = beta_data, method = 'sgv')
R1_SR

R2_SR <- r2glmm::r2beta(m2_1,  partial = T, data = beta_data, method = 'sgv')
R2_SR

R_SR <- R2_SR %>%
  filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2" | 
           Effect=="Model") %>%  
  bind_rows(R1_SR %>% filter(!Effect=="Model"))

write.csv(R_SR,  file = "results/R2_beta_SR.csv")

# (2) beta ENSPIE -----


# Data Exploration----

str(beta_data)

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(beta_data, aes(pca1_clima, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(Prec_Varieb, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(cover_litter, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(pH, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(Corg_percent, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(grazing_intencity, beta_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(beta_data, aes(mowing, beta_100_ENSPIE)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
  scale_color_manual(values = col)+
  labs(color='Habitat type')



# GLLM----

## Exploration----

# Check the model:
m_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                   poly(pca1_clima, 2) +
                   poly(pH, 2) +
                   poly(Corg_percent,2)+
                   poly(cover_litter,2) +
                   grazing_intencity + mowing +
                   (1|dataset),  data = beta_data)

# check model
plot(m_ENSPIE) # heteroscedasticity
qqnorm(resid(m_ENSPIE))
qqline(resid(m_ENSPIE))



Anova(m_ENSPIE)


# Model selection ----

## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 

# poly(pca1_clima, 2) 
# poly(pH, 2) is marginal

m1_1_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(pH, 2) +
                     poly(Corg_percent,2) +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)

m1_2_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)


m1_3_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      pca1_clima +
                      pH +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)


m1_4_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      pca1_clima +
                      pH +
                      Corg_percent+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)


# calculate and compare AIC 
AIC(m1_1_ENSPIE, m1_2_ENSPIE, m1_3_ENSPIE, m1_4_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m1_1_ENSPIE)
# Anova(m1_4_ENSPIE)


## Model 2: Add Prec_Varieb ----

m2_1_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)


m2_2_ENSPIE <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      Prec_Varieb +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)



# calculate and compare AIC 
AIC(m2_1_ENSPIE, m2_2_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))


Anova(m2_1_ENSPIE)



# Plots----
#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")

set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,legend.pos = "None", geom.linetype = 2)

# final models:
Anova(m1_1_ENSPIE)
Anova(m2_1_ENSPIE)


## Clima ----

clima_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_beta_ENSPIE

Fig.betaENSPIE_clima <- ggplot(clima_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(pca1_clima, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Climate gradient')+ 
  geom_line(linetype=1, size=0.5, col="#00AC7F") 

Fig.betaENSPIE_clima 


plot_model(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]", #  show.data=F
                                    title = "", line.size=0.5) + aes(linetype="solid") +
  labs(y=expression(paste("beta ENS"[PIE])), x='Climate gradient (PC)')+ 
  geom_point(data = beta_data, aes(pca1_clima, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 



## Soil C ----
Humus_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_beta_ENSPIE

Fig.betaENSPIE_soilC <- ggplot(Humus_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Corg_percent, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Soil C')+ 
  geom_line(linetype=1, size=1, col="#00AC7F") 

Fig.betaENSPIE_soilC 

plot_model(m1_1_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]",# show.data=F,
                                    title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("beta ENS"[PIE])), x='Soil C')+ 
  geom_point(data = beta_data, aes(Corg_percent, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 
# xlim(-1.9,5)



## Litter % ----
Litter_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_beta_ENSPIE

Fig.betaENSPIE_Litter <- ggplot(Litter_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(cover_litter, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Litter cover')+ 
  geom_line(linetype=1, size=1, col="#00AC7F") 

Fig.betaENSPIE_Litter 

plot_model(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]",# show.data=F,
                                     title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("beta ENS"[PIE])), x='Litter cover')+ 
  geom_point(data = beta_data, aes(cover_litter, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

 


## Soil pH ----
pH_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="pH[3.7:9, by=.001]")
pH_pred_beta_ENSPIE

Fig.betaENSPIE_soil.pH <- ggplot(pH_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(pH, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Soil pH')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaENSPIE_soil.pH 

plot_model(m1_1_ENSPIE,type = "pred", terms="pH[3.7:9, by=.001]",# show.data=F,
                                      title = "", line.size=0.5)  +
  labs(y=expression(paste("beta ENS"[PIE])), x='Soil pH')+ 
  geom_point(data = beta_data, aes(pH, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 


## Grazing ----
grazing_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_beta_ENSPIE

Fig.betaENSPIE_grazing <- ggplot(grazing_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(grazing_intencity, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Grazing intencity')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaENSPIE_grazing

plot_model(m1_1_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]",# show.data=F,
                                      title = "", line.size=0.5)  +   
  labs(y=expression(paste("beta ENS"[PIE])), x='Grazing intencity')+ 
  geom_point(data = beta_data, aes(grazing_intencity, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8,
             position=position_jitter(w=0.2))+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 



## Mowing ----

plot_model(m1_1_ENSPIE,type = "pred", terms="mowing",  show.data=T,
           title = "", line.size=0.5)

Fig.betaENSPIE_mowing <-ggplot(beta_data, aes(mowing, beta_100_ENSPIE)) + 
  geom_boxplot(color="#00AC7F")+
  labs(y=expression(paste("ENS"[PIE])), x='Mowing')+ 
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) +
  labs(color='Habitat type')

Fig.betaENSPIE_mowing 


## Prec_Varieb ----
precipCV_pred_beta_ENSPIE <- get_model_data(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb")
precipCV_pred_beta_ENSPIE

Fig.betaENSPIE_precip.CV <- ggplot(precipCV_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=beta_data, aes(Prec_Varieb, beta_100_ENSPIE, fill = habitat, col=habitat), size=1, alpha=0.8, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Precipitation CV')+ 
  geom_line(linetype=5, size=0.5, col="#00AC7F") 

Fig.betaENSPIE_precip.CV

plot_model(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb", #  show.data=F
                                        title = "", line.size=0.5) + 
  labs(y=expression(paste("beta ENS"[PIE])), x='Precipitation variability')+ 
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type')


# Plots combined ----

# for Gamma SR (100 m^2)

Fig.betaENSPIE_clima + 
  Fig.betaENSPIE_precip.CV +
  Fig.betaENSPIE_soilC + 
  Fig.betaENSPIE_soil.pH + 
  Fig.betaENSPIE_Litter + 
  Fig.betaENSPIE_grazing +
  plot_layout(ncol=3)


# combine with SR


Fig.betaSR_clima + Fig.betaENSPIE_clima +
  Fig.betaSR_precip.CV + Fig.betaENSPIE_precip.CV +
  Fig.betaSR_soilC + Fig.betaENSPIE_soilC +
  Fig.betaSR_soil.pH + Fig.betaENSPIE_soil.pH +
  Fig.betaSR_Litter + Fig.betaENSPIE_Litter +
  Fig.betaSR_grazing + Fig.betaENSPIE_grazing +
  plot_layout(ncol=2)+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face = 'bold'), plot.tag.position = c(0.22, 0.95),
        plot.margin = unit(c(0,0,0,0), "pt"))



  
  
  
# final statistics -----
Anova(m1_1_ENSPIE)
Anova(m2_1_ENSPIE)

# write.csv(Anova(m1_1_ENSPIE),  file = "results/glmer_beta_ENSPIE.csv")
# write.csv(coef(summary(m1_1_ENSPIE)),  file = "results/summary_beta_ENSPIE.csv")
# write.csv(Anova(m2_1_ENSPIE),  file = "results/glmer_beta_ENSPIE_2.csv")
# write.csv(coef(summary(m2_1_ENSPIE)),  file = "results/summary_beta_ENSPIE_2.csv")


# R2 for the entire model---------

# write.csv(MuMIn::r.squaredGLMM(m1_1_ENSPIE),  file = "results/Mod1_R2_beta_ENSPIE.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1_ENSPIE),  file = "results/Mod2_R2_beta_ENSPIE.csv")



R1_ENSPIE <- r2glmm::r2beta(m1_1_ENSPIE,  partial = T, data = beta_data, method = 'sgv')
R1_ENSPIE

R2_ENSPIE <- r2glmm::r2beta(m2_1_ENSPIE,  partial = T, data = beta_data, method = 'sgv')
R2_ENSPIE


R_ENSPIE <- R2_ENSPIE %>%
  filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2" | 
           Effect=="Model") %>% 
  bind_rows(R1_ENSPIE %>% filter(!Effect=="Model")) 

write.csv(R_ENSPIE,  file = "results/R2_beta_ENSPIE.csv")


