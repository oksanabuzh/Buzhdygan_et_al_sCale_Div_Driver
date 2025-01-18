# Purpose: GLMM analysis for the 10-m2 plots (alpha diversity)


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

# prepare subset of data for alpha scale (10-m2 plots)

alpha <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type=="alpha")%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header, 
            by=c("dataset", "plotID", "series", "subplot")
            ) %>% 
  mutate(dataset=factor(dataset))
                   
str(alpha) 
names (alpha)

# dataset is a separate vegetation survey campaign
alpha$dataset 


# selected variables, removed NAs
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

str(alpha_data)
alpha_data$dataset


# Analysis -----------------

#---------------------------#
# (1) Species richness -----
# --------------------------#

##  Data Exploration----
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

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(alpha_mean, aes(pca1_clima, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(Prec_Varieb, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(cover_litter, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(pH, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(Corg_percent, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(grazing_intencity, alpha_10_div)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(mowing, alpha_10_div)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
 scale_color_manual(values = col)+
  labs(color='Habitat type')



## GLLM----

### Exploration----

m <- glmer (alpha_10_div ~ 
              poly(pca1_clima, 2) +
              poly(Corg_percent,2)+
              poly(pH, 2) +
              poly(cover_litter,2) +
              grazing_intencity + mowing +
              (1|dataset/series),  
            family = "poisson", data = alpha_data)


check_convergence(m)
  
# check model
plot(m)
qqnorm(resid(m))
qqline(resid(m))

# check multicolinearity
check_collinearity(m)
  

# check overdispersion
sum(residuals(m, type = "pearson")^2) / df.residual(m)
  

Anova(m)
summary(m)

# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m)
  
# Partial R2 for fixed effects
r2glmm::r2beta(m,  partial = T, data=alpha_data)
  
  
  
### Test random effects ----
  
# check random effects
ranef(m) # not zeros
hist(ranef(m)$`series:dataset`[,1])
  
RE_m1 <- glmer (alpha_10_div ~ 
                  poly(pca1_clima, 2) +
                  poly(Corg_percent,2)+
                  poly(pH, 2) +
                  poly(cover_litter,2) +
                  grazing_intencity + mowing +
                  (1|dataset/series),  
                family = "poisson", data = alpha_data)

RE_m2 <- glmer (alpha_10_div ~ 
                  poly(pca1_clima, 2) +
                  poly(Corg_percent,2)+
                  poly(pH, 2) +
                  poly(cover_litter,2) +
                  grazing_intencity + mowing +
                  (1|dataset),  
                family = "poisson", data = alpha_data)
  
anova(RE_m1, RE_m2)
# RE_m1 significantly differ from RE_m2, thus we should add series nested in dataset 
# as a random effect  

# Model selection-----
## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 

# poly(pca1_clima, 2) is marginal
# poly(pH, 2) is non significant

m1_1 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Corg_percent,2)+
                 poly(pH, 2) +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

m1_2 <- glmer (alpha_10_div ~ 
                 pca1_clima +
                 poly(Corg_percent,2)+
                 poly(pH, 2) +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

m1_3 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)


m1_4 <- glmer (alpha_10_div ~ 
                 pca1_clima +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)


# calculate and compare AIC 
AIC(m1_1, m1_2, m1_3, m1_4) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

# best models (within 2 units of the model with lowest AIC)
Anova(m1_3)
# Anova(m1_1) # lowest AIC

## Model 2: Add Prec_Varieb ----

m2_1 <- glmer (alpha_10_div ~ 
                poly(pca1_clima, 2) +
                poly(Prec_Varieb, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

m2_2 <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 Prec_Varieb +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)

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


# Final models
Anova(m1_3) # for all variables
Anova(m2_1) # for precipitation variability

## Clima ----

# Predictions
clima_pred_10m <-get_model_data(m1_3,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_10m

clima_pred_10m$x
clima_pred_10m$predicted
clima_pred_10m$conf.low
clima_pred_10m$conf.high

Fig.alphaSR_clima <- ggplot(clima_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(pca1_clima, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Climate gradient (PC)')+
  geom_line(linetype=1, size=1, col="#64ABCE") 

Fig.alphaSR_clima 


## Soil C ----

Humus_pred_10m <-get_model_data(m1_3,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_10m

Fig.alphaSR_soilC <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, 
             aes(Corg_percent, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Soil C')+ 
  geom_line(linetype=1, size=1, col="#64ABCE")

Fig.alphaSR_soilC 


## Litter % ----

Litter_pred_10m <-get_model_data(m1_3,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_10m

Fig.alphaSR_Litter <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(cover_litter, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Litter cover')+
  geom_line(linetype=1, size=1, col="#64ABCE")

Fig.alphaSR_Litter


## Soil pH ----

pH_pred_10m <-get_model_data(m1_1,type = "pred", terms="pH[3.8:9, by=.001]")
pH_pred_10m

Fig.alphaSR_soil.pH <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(pH, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Soil pH')+
  geom_line(linetype=5, size=0.5, col="#64ABCE")

Fig.alphaSR_soil.pH

plot_model(m1_1,type = "pred", terms="pH[3.7:9, by=.001]",# show.data=F,
           title = "", line.size=0.5) +  aes(linetype="solid") +
  labs(y="gamma SR", x='Soil pH')+ 
  geom_point(data = alpha_data, aes(pH, alpha_10_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 



## Grazing ----

grazing_pred_10m <-get_model_data(m1_3,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_10m

Fig.alphaSR_grazing <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(grazing_intencity, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col) +  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#64ABCE")

Fig.alphaSR_grazing


## Mowing ----

plot_model(m1_3,type = "pred", terms="mowing",  show.data=T,
                                  title = "", line.size=0.5)

Fig.alphaSR_mowing <-ggplot(alpha_mean, aes(mowing, alpha_10_div)) + 
  geom_boxplot(color="#64ABCE")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
labs(y="Species richness", x='Mowing')

Fig.alphaSR_mowing 


## Prec_Varieb ----
precipCV_pred_10m <-get_model_data(m2_1,type = "pred", terms="Prec_Varieb")
precipCV_pred_10m

Fig.alphaSR_precip.CV <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(Prec_Varieb, alpha_10_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Precipitation variability')+
  geom_line(linetype=1, size=1, col="#64ABCE")

Fig.alphaSR_precip.CV

plot_model(m2_1,type = "pred", terms="Prec_Varieb", #  show.data=F
                                title = "", line.size=0.5) + aes(linetype="solid") +
  labs(y="alpha SR ", x='Precipitation variability')+ 
  geom_point(data = alpha_mean, aes(Prec_Varieb, alpha_10_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

# Plots combined ----
# for alpha SR (10 m^2)

Fig.alphaSR_clima + 
  Fig.alphaSR_precip.CV +
  Fig.alphaSR_soilC + 
  Fig.alphaSR_soil.pH + 
  Fig.alphaSR_Litter + 
  Fig.alphaSR_grazing +
  plot_layout(ncol=3)


# R2 for the entire model---------

# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1_3)
MuMIn::r.squaredGLMM(m2_1)
# write.csv(MuMIn::r.squaredGLMM(m1_3),  file = "results/Mod1_R2_alpha_SR.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1),  file = "results/Mod2_R2_alpha_SR.csv")


# Partial R2 for fixed effects
Anova(m1_3) # fo all variables
Anova(m2_1) # for Prec_Varieb

R1 <- r2glmm::r2beta(m1_3,  partial = T, data = alpha_data, method = 'sgv')
R1

R2 <- r2glmm::r2beta(m2_1,  partial = T, data = alpha_data, method = 'sgv')
R2


R <- R2 %>%
  filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2"| Effect=="Model")%>% 
    bind_rows(R1 %>% filter(!Effect=="Model")) 

write.csv(R,  file = "results/R2_alpha_SR.csv")

# write.csv(MuMIn::r.squaredGLMM(m1_3),  file = "results/Mod_R2_alpha_SR.csv")
# write.csv(Anova(m1_3),  file = "results/glmer_alpha_SR.csv")
# write.csv(coef(summary(m1_3)),  file = "results/summary_alpha_SR.csv")
# write.csv(Anova(m2_1),  file = "results/glmer_alpha_SR_2.csv")
# write.csv(coef(summary(m2_1)),  file = "results/summary_alpha_SR_2.csv")



# (2) ENSPIE -----
# Data Exploration----

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(alpha_mean, aes(pca1_clima, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(Prec_Varieb, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(cover_litter, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(pH, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(Corg_percent, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(grazing_intencity, alpha_10_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(alpha_mean, aes(mowing, alpha_10_ENSPIE)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
  scale_color_manual(values = col)+
  labs(color='Habitat type')



# GLLM----

## Exploration----

m_ENSPIE <- lmer(alpha_10_ENSPIE ~ 
                   poly(pca1_clima, 2) +
                   poly(pH, 2) +
                   poly(Corg_percent,2)+
                   poly(cover_litter,2) +
                   grazing_intencity + mowing +
                   (1|dataset/series),  data = alpha_data)

# check model
plot(m_ENSPIE) # heteroscedasticity
qqnorm(resid(m_ENSPIE))
qqline(resid(m_ENSPIE))

# log-transform response
m_ENSPIE_b <- lmer(log(alpha_10_ENSPIE) ~ 
                     poly(pca1_clima, 2) +
                     poly(pH, 2) +
                     poly(Corg_percent,2)+
                     poly(cover_litter,2) +
                     grazing_intencity + mowing +
                     (1|dataset/series),  data = alpha_data)

plot(m_ENSPIE_b) # better
qqnorm(resid(m_ENSPIE_b))
qqline(resid(m_ENSPIE_b))

check_collinearity(m_ENSPIE_b)

Anova(m_ENSPIE_b)

# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_ENSPIE_b)

# Partial R2 for fixed effects
r2glmm::r2beta(m_ENSPIE_b,  partial = T, data=alpha_data)

# Model selection----

## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 

Anova(m_ENSPIE_b)

# poly(cover_litter, 2) is marginal

m1_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

m1_2_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      Corg_percent +
                      cover_litter +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)



m1_3_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(pH,2) +
                      Corg_percent +
                      cover_litter +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

m1_4_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      cover_litter +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

 
# calculate and compare AIC 
AIC(m1_1_ENSPIE, m1_2_ENSPIE, m1_3_ENSPIE, m1_4_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m1_1_ENSPIE )

## Model 2: Add Prec_Varieb ----

### Model selection -----

m2_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)


m2_2_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      Prec_Varieb +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)



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

# Final models
Anova(m1_1_ENSPIE) # for all variables
Anova(m2_1_ENSPIE) # for precipitation variability


get_model_data(m1_1_ENSPIE, type = "pred")

## Clima ----

# Predictions
clima_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_10m_Ensp

Fig.alphaENSPIE_clima <- ggplot(clima_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(pca1_clima, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Climate gradient (PC)')+ 
  geom_line(linetype=1, size=1, col="#64ABCE") 

Fig.alphaENSPIE_clima 

plot_model(m1_1_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]",  # show.data=T,
                                    title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("alpha ENS"[PIE])), x='Climate gradient (PC)')+ 
  geom_point(data = alpha_mean, aes(pca1_clima, alpha_10_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 


## Soil C ----

Humus_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_10m_Ensp

Fig.alphaENSPIE_soilC <- ggplot(Humus_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(Corg_percent, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Soil C')+ 
  geom_line(linetype=5, size=0.5, col="#64ABCE") 

Fig.alphaENSPIE_soilC 



## Litter % ----
Litter_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_10m_Ensp

Fig.alphaENSPIE_Litter <- ggplot(Litter_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(cover_litter, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Litter cover')+ 
  geom_line(linetype=1, size=0.5, col="#64ABCE") 

Fig.alphaENSPIE_Litter 

plot_model(m1_1_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]",# show.data=F,
                                     title = "", line.size=0.5) + aes(linetype="solid") +
  labs(y=expression(paste("alpha ENS"[PIE])), x='Litter cover')+ 
  geom_point(data = alpha_mean, aes(cover_litter, alpha_10_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 



## Soil pH ----
pH_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="pH[3.8:9, by=.001]")
pH_pred_10m_Ensp

Fig.alphaENSPIE_soil.pH <- ggplot(pH_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(pH, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Soil pH')+ 
  geom_line(linetype=5, size=0.5, col="#64ABCE") 

Fig.alphaENSPIE_soil.pH 


## Grazing ----
grazing_pred_10m_Ensp <-get_model_data(m1_1_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_10m_Ensp

Fig.alphaENSPIE_grazing <- ggplot(grazing_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(grazing_intencity, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Grazing intencity')+ 
  geom_line(linetype=5, size=0.5, col="#64ABCE") 

Fig.alphaENSPIE_grazing 


## Mowing ----

plot_model(m1_1_ENSPIE,type = "pred", terms="mowing",  show.data=T,
           title = "", line.size=0.5)

Fig.alphaENSPIE_mowing <-ggplot(alpha_mean, aes(mowing, alpha_10_ENSPIE)) + 
  geom_boxplot(color="#64ABCE")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
  labs(y=expression(paste("ENS"[PIE])), x='Mowing')+ labs(color='Habitat type')

Fig.alphaENSPIE_mowing 


## Prec_Varieb ----
precipCV_pred_10m_Ensp <-get_model_data(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb")
precipCV_pred_10m_Ensp

Fig.alphaENSPIE_precip.CV <- ggplot(precipCV_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=alpha_mean, aes(Prec_Varieb, alpha_10_ENSPIE, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y=expression(paste("ENS"[PIE])), x='Precipitation variability')+ 
  geom_line(linetype=1, size=1, col="#64ABCE") 

Fig.alphaENSPIE_precip.CV 


## Plots combined ----
# for alpha ENSPIE (10 m^2)

Fig.alphaENSPIE_clima + 
  Fig.alphaENSPIE_precip.CV +
  Fig.alphaENSPIE_soilC + 
  Fig.alphaENSPIE_soil.pH + 
  Fig.alphaENSPIE_Litter + 
  Fig.alphaENSPIE_grazing +
  plot_layout(ncol=3)


# R2 for the entire model---------

# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1_1_ENSPIE)
MuMIn::r.squaredGLMM(m2_1_ENSPIE)

# write.csv(MuMIn::r.squaredGLMM(m1_1_ENSPIE),  file = "results/Mod1_R2_alpha_ENSPIE.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1_ENSPIE),  file = "results/Mod2_R2_alpha_ENSPIE.csv")

# Partial R2 for fixed effects
Anova(m1_1_ENSPIE) # fo all variables
Anova(m2_1_ENSPIE) # for Prec_Varieb


R1_ENSPIE <- r2glmm::r2beta(m1_1_ENSPIE,  partial = T, data = alpha_data, method = 'sgv')
R1_ENSPIE

R2_ENSPIE <- r2glmm::r2beta(m2_1_ENSPIE,  partial = T, data = alpha_data, method = 'sgv')
R2_ENSPIE


R_ENSPIE <- R2_ENSPIE %>%
              filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2" |
                       Effect=="Model") %>% 
  bind_rows(R1_ENSPIE %>% 
  filter(!Effect=="Model")
  ) 

write.csv(R_ENSPIE,  file = "results/R2_alpha_ENSPIE.csv")

# write.csv(Anova(m1_1_ENSPIE),  file = "results/glmer_alpha_ENSPIE.csv")
# write.csv(coef(summary(m1_1_ENSPIE)),  file = "results/summary_alpha_ENSPIE.csv")
# write.csv(Anova(m2_1_ENSPIE),  file = "results/glmer_alpha_ENSPIE_2.csv")
# write.csv(coef(summary(m2_1_ENSPIE)),  file = "results/summary_alpha_ENSPIE_2.csv")

