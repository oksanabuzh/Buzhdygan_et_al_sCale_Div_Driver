# Purpose: GLMM analysis for the 100 m^2 plots (gamma diversity)


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


# prepare subset of data for gamma scale (100-m2 plots)

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
gamma_data <- beta_gamma %>% 
  dplyr::select(gamma_100_div, gamma_100_ENSPIE, 
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
  mutate(habitat =fct_relevel(habitat_broad, 
                              c("saline", "complex", "dry",
                                "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na


# (1) Species Richness -----


# Data Exploration----

str(gamma_data)

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(gamma_data, aes(pca1_clima, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(Prec_Varieb, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(cover_litter, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(pH, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(Corg_percent, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(grazing_intencity, gamma_100_div)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(mowing, gamma_100_div)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
 scale_color_manual(values = col)+
  labs(color='Habitat type')



# GLLM----

## Exploration----

m <- glmer (gamma_100_div ~ 
              poly(pca1_clima, 2) +
              poly(Corg_percent,2)+
              poly(pH, 2) +
              poly(cover_litter,2) +
              grazing_intencity + mowing +
              (1|dataset),  
            family = "poisson", data = gamma_data)


check_convergence(m)
  
# check model
plot(m)
qqnorm(resid(m))
qqline(resid(m))

# check multicolinearity
check_collinearity(m)
  

# check overdispersion
sum(residuals(m, type = "pearson")^2) / df.residual(m)
  

# high overdispersion, use negative binomial

m_b <- glmer.nb (gamma_100_div ~ 
                   poly(pca1_clima, 2) +
                   poly(Corg_percent,2)+
                   poly(pH, 2) +
                   poly(cover_litter,2) +
                   grazing_intencity + mowing + 
                   (1|dataset), data = gamma_data)

sum(residuals(m_b, type = "pearson")^2) / df.residual(m_b)
# corrected for the overdispersion

check_collinearity(m_b)

Anova(m_b)
summary(m_b)




# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_b)
  
# Partial R2 for fixed effects
r2glmm::r2beta(m_b,  partial = T, data=gamma_data)
  
  
  
## Random effects ----
ranef(m_b) # not zeros
hist(ranef(m_b)$`dataset`[,1])


Anova(m_b)
summary(m_b)

# Model selection ----
## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 
# poly(pca1_clima, 2) 
# poly(pH, 2) is marginal

m1_1 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

m1_2 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Corg_percent,2)+
                    pH +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)


# calculate and compare AIC 
AIC(m1_1, m1_2) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))


Anova(m1_1)
# Anova(m1_2)

## Model 2: Add Prec_Varieb ----
m2_1 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Prec_Varieb, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

m2_2 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    Prec_Varieb +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)


# calculate and compare AIC 
AIC(m2_1, m2_2) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

# Anova(m2_2)
Anova(m2_1)


## -> Additional analysis for precipitation CV effects  -----
# Analysis whether precipitation variability adds explanatory power beyond the 
# nonlinear effect of the climate gradient:
# Model m3_1 includes both linear and quadratic terms for the climate gradient. 
# Model m3_2 includes the linear effects of both climate gradient and 
# precipitation variability (i.e., the quadratic term for the climate gradient was replaced by precipitation variability).
# If the latter model is better (e.g. AIC is smaller than two units), we have more support for claim that the precipitation variability effect is shown.


m3_1 <- glmer.nb (gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    (1|dataset), data = gamma_data)


m3_2 <- glmer.nb (gamma_100_div ~ 
                    pca1_clima + Prec_Varieb +
                    (1|dataset), data = gamma_data)

# calculate and compare AIC 
AIC(m3_1, m3_2) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m3_1)
Anova(m3_2)

# Plots----
#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")

set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,legend.pos = "None", geom.linetype = 2)

# final models:
Anova(m1_1)
Anova(m2_1)

## Clima ----

clima_pred_100m <- get_model_data(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]")
clima_pred_100m

clima_pred_10m$x
clima_pred_10m$predicted
clima_pred_10m$conf.low
clima_pred_10m$conf.high

Fig.gammaSR_clima <- ggplot(clima_pred_100m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=gamma_data, aes(pca1_clima, gamma_100_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Climate gradient (PC)')+
  geom_line(linetype=1, size=1, col="#D6604D") 

Fig.gammaSR_clima 



plot_model(m1_1,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]", #  show.data=F
           title = "", line.size=1) + aes(linetype="solid") +
  labs(y="gamma SR", x='Climate gradient (PC)')+ 
  geom_point(data = gamma_data, aes(pca1_clima, gamma_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 


## Soil C ----
Humus_pred_100m <- get_model_data(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]")
Humus_pred_100m

Fig.gammaSR_soilC <- ggplot(Humus_pred_100m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=gamma_data, aes(Corg_percent, gamma_100_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Soil C')+
  geom_line(linetype=1, size=1, col="#D6604D") 

Fig.gammaSR_soilC

plot_model(m1_1,type = "pred", terms="Corg_percent[0:9.5, by=.001]",# show.data=F,
                                title = "", line.size=1) + aes(linetype="solid") +
  labs(y="gamma SR", x='Soil C')+ 
  geom_point(data = gamma_data, aes(Corg_percent, gamma_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 



## Litter % ----
Litter_pred_100m <- get_model_data(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]")
Litter_pred_100m

Fig.gammaSR_Litter <- ggplot(Litter_pred_100m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=gamma_data, aes(cover_litter, gamma_100_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Litter cover')+
  geom_line(linetype=1, size=1, col="#D6604D") 

Fig.gammaSR_Litter

plot_model(m1_1,type = "pred", terms="cover_litter[0:100, by=0.01]",# show.data=F,
                                  title = "", line.size=1) + aes(linetype="solid") +
  labs(y="gamma SR", x='Litter cover')+ 
  geom_point(data = gamma_data, aes(cover_litter, gamma_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 


## Soil pH ----
pH_pred_100m <- get_model_data(m1_1,type = "pred", terms="pH[3.7:9, by=.001]")
pH_pred_100m

Fig.gammaSR_soil.pH <- ggplot(pH_pred_100m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=gamma_data, aes(pH, gamma_100_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Soil pH')+
  geom_line(linetype=1, size=0.5, col="#D6604D") 

Fig.gammaSR_soil.pH

plot_model(m1_1,type = "pred", terms="pH[3.7:9, by=.001]",# show.data=F,
                                title = "", line.size=0.5) +  aes(linetype="solid") +
  labs(y="gamma SR", x='Soil pH')+ 
  geom_point(data = gamma_data, aes(pH, gamma_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 



## Grazing ----
grazing_pred_100m <- get_model_data(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]")
grazing_pred_100m

Fig.gammaSR_grazing <- ggplot(grazing_pred_100m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=gamma_data, aes(grazing_intencity, gamma_100_div, fill = habitat, col=habitat), size=3, alpha=0.7, pch=21)+
  scale_fill_manual(values = col)+  scale_color_manual(values = col)  +  
  labs(y="Species richness", x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#D6604D") 

Fig.gammaSR_grazing

plot_model(m1_1,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]",# show.data=F,
                                  title = "", line.size=0.5)  +  
  labs(y="gamma SR", x='Grazing intencity')+ 
  geom_point(data = gamma_data, aes(grazing_intencity, gamma_100_div, color = habitat), size=3, alpha=0.8,
               position=position_jitter(w=0.2))+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

## Mowing ----

plot_model(m1_1,type = "pred", terms="mowing",  show.data=T,
                                  title = "", line.size=0.5)

Fig.gammaSR_mowing <-ggplot(gamma_data, aes(mowing, gamma_100_div)) + 
  geom_boxplot(color="#D6604D")+
  geom_point(aes(color = habitat, fill = habitat), pch=21, position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+ scale_fill_manual(values = col) + 
  labs(y="Species richness", x='Mowing')

Fig.gammaSR_mowing 


## Prec_Varieb ----
Fig.gammaSR_precip.CV <- plot_model(m2_1,type = "pred", terms="Prec_Varieb", #  show.data=F
                                title = "", line.size=1) + aes(linetype="solid") +
  labs(y="gamma SR", x='Precipitation variability')+ 
  geom_point(data = gamma_data, aes(Prec_Varieb, gamma_100_div, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

Fig.gammaSR_precip.CV 



# Plots combined ----

# for Gamma SR (100 m^2)

Fig.gammaSR_clima + 
  Fig.gammaSR_precip.CV +
  Fig.gammaSR_soilC + 
  Fig.gammaSR_soil.pH + 
  Fig.gammaSR_Litter + 
  Fig.gammaSR_grazing +
  plot_layout(ncol=3)



# combine with alpha
# from "analysis/02_GLMM_alpha_10_div.R"
Fig.alphaSR_clima + Fig.gammaSR_clima +
  Fig.alphaSR_precip.CV + Fig.gammaSR_precip.CV +
  Fig.alphaSR_soilC + Fig.gammaSR_soilC +
  Fig.alphaSR_soil.pH + Fig.gammaSR_soil.pH +
  Fig.alphaSR_Litter + Fig.gammaSR_Litter +
  Fig.alphaSR_grazing + Fig.gammaSR_grazing +
  plot_layout(ncol=2)+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face = 'bold'), plot.tag.position = c(0.22, 0.95),
        plot.margin = unit(c(0,0,0,0), "pt"))

# R2 for the entire model---------

# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1_1)
MuMIn::r.squaredGLMM(m2_1)

# write.csv(MuMIn::r.squaredGLMM(m1_1),  file = "results/Mod1_R2_gamma_SR.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1),  file = "results/Mod2_R2_gamma_SR.csv")


# Partial R2 for fixed effects

Anova(m1_1)
Anova(m2_1) # for Prec_Varieb


R1 <- r2glmm::r2beta(m1_1,  partial = T, data = gamma_data, method = 'sgv')
R1

R2 <- r2glmm::r2beta(m2_1,  partial = T, data = gamma_data, method = 'sgv')
R2




R <- R2 %>%
  filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2" | 
  Effect=="Model") %>%  
    bind_rows(R1 %>% filter(!Effect=="Model"))

write.csv(R,  file = "results/R2_gamma_SR.csv")

# write.csv(Anova(m1_1),  file = "results/glmer_gamma_SR.csv")
# write.csv(coef(summary(m1_1)),  file = "results/summary_gamma_SR.csv")
# write.csv(Anova(m2_1),  file = "results/glmer_gamma_SR_2.csv")
# write.csv(coef(summary(m2_1)),  file = "results/summary_gamma_SR_2.csv")



# (2) ENSPIE -----


# Data Exploration----

str(gamma_data)

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")


ggplot(gamma_data, aes(pca1_clima, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(Prec_Varieb, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(pca1_clima, Prec_Varieb)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(cover_litter, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(pH, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(Corg_percent, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat)) + scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(grazing_intencity, gamma_100_ENSPIE)) + 
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.4)) + 
  scale_color_manual(values = col)+
  labs(color='Habitat type')

ggplot(gamma_data, aes(mowing, gamma_100_ENSPIE)) + 
  geom_boxplot()+
  geom_point(size=2, aes(color = habitat), position=position_jitter(w=0.1))+
  scale_color_manual(values = col)+
  labs(color='Habitat type')



# GLLM----

## Exploration----

# Check the model:
m_ENSPIE <- lmer(gamma_100_ENSPIE ~ 
                   poly(pca1_clima, 2) +
                   poly(pH, 2) +
                   poly(Corg_percent,2)+
                   poly(cover_litter,2) +
                   grazing_intencity + mowing +
                   (1|dataset),  data = gamma_data)

# check model
plot(m_ENSPIE) # heteroscedasticity
qqnorm(resid(m_ENSPIE))
qqline(resid(m_ENSPIE))

m_ENSPIE_b <- lmer(log(gamma_100_ENSPIE) ~ 
                     poly(pca1_clima, 2) +
                     poly(pH, 2) +
                     poly(Corg_percent,2)+
                     poly(cover_litter,2) +
                     grazing_intencity + mowing +
                     (1|dataset),  data = gamma_data)

plot(m_ENSPIE_b) # better
qqnorm(resid(m_ENSPIE_b))
qqline(resid(m_ENSPIE_b))

check_collinearity(m_ENSPIE_b)

Anova(m_ENSPIE_b)


plot_model(m_ENSPIE_b,type = "pred", terms="grazing_intencity", show.data=T,
           title = "", line.size=1) + aes(linetype="solid")

# Model selection ----

## Model 1: all predictors (except precipitation CV)----
# test quadratic effects of climate, soil C, pH, and litter 

Anova(m_ENSPIE_b)

# poly(pca1_clima, 2) 
# poly(pH, 2) is marginal

m1_1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      Corg_percent+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)

m1_2_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)


m1_3_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)

# calculate and compare AIC 
AIC(m1_1_ENSPIE, m1_2_ENSPIE, m1_3_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))


Anova(m1_3_ENSPIE)


## Model 2: Add Prec_Varieb ----

### Model selection -----
m2_1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)


m2_2_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      Prec_Varieb +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)



# calculate and compare AIC 
AIC(m2_1_ENSPIE, m2_2_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))


# Anova(m2_2)
Anova(m2_1_ENSPIE)


## -> Additional analysis for precipitation CV effects  -----
# Analysis whether precipitation variability adds explanatory power beyond the 
# nonlinear effect of the climate gradient:
# Model m3_1_ENSPIE includes both linear and quadratic terms for the climate gradient. 
# Model m3_2_ENSPIE includes the linear effects of both climate gradient and 
# precipitation variability (i.e., the quadratic term for the climate gradient was replaced by precipitation variability).
# If the latter model is better (e.g. AIC is smaller than two units), we have more support for claim that the precipitation variability effect is shown.


m3_1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      (1|dataset),  data = gamma_data)

m3_2_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ 
                      pca1_clima +
                      Prec_Varieb +
                      (1|dataset),  data = gamma_data)

# calculate and compare AIC 
AIC(m3_1_ENSPIE, m3_2_ENSPIE) %>% 
  arrange(+AIC) %>% 
  mutate(delta_AIC=AIC-min(AIC))

Anova(m3_1_ENSPIE)
Anova(m3_2_ENSPIE)

# Plots----
#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600",  "#e3c28b", "#CC99FF", "#0066FF" ,  "#00B200",  "#006600")

set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,legend.pos = "None", geom.linetype = 2)

# final models:
Anova(m1_3_ENSPIE)
Anova(m2_1_ENSPIE)

## Clima ----

Fig.gammaENSPIE_clima <- plot_model(m1_3_ENSPIE,type = "pred", terms="pca1_clima[-1.2:4.8, by=.001]", #  show.data=F
                                    title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("gamma ENS"[PIE])), x='Climate gradient (PC)')+ 
  geom_point(data = gamma_data, aes(pca1_clima, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

Fig.gammaENSPIE_clima 


## Soil C ----
Fig.gammaENSPIE_soilC <- plot_model(m1_3_ENSPIE,type = "pred", terms="Corg_percent[0:9.5, by=.001]",# show.data=F,
                                    title = "", line.size=0.5)  +
  labs(y=expression(paste("gamma ENS"[PIE])), x='Soil C')+ 
  geom_point(data = gamma_data, aes(Corg_percent, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 
# xlim(-1.9,5)

Fig.gammaENSPIE_soilC 


## Litter % ----

Fig.gammaENSPIE_Litter <- plot_model(m1_3_ENSPIE,type = "pred", terms="cover_litter[0:100, by=0.01]",# show.data=F,
                                     title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("gamma ENS"[PIE])), x='Litter cover')+ 
  geom_point(data = gamma_data, aes(cover_litter, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

Fig.gammaENSPIE_Litter 


## Soil pH ----

Fig.gammaENSPIE_soil.pH <- plot_model(m1_3_ENSPIE,type = "pred", terms="pH[3.7:9, by=.001]",# show.data=F,
                                      title = "", line.size=0.5)  +
  labs(y=expression(paste("gamma ENS"[PIE])), x='Soil pH')+ 
  geom_point(data = gamma_data, aes(pH, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

Fig.gammaENSPIE_soil.pH 

## Grazing ----

Fig.gammaENSPIE_grazing <- plot_model(m1_3_ENSPIE,type = "pred", terms="grazing_intencity[-0.2:3.2, by=0.01]",# show.data=F,
                                      title = "", line.size=0.5)  +   aes(linetype="solid") + 
  labs(y=expression(paste("gamma ENS"[PIE])), x='Grazing intencity')+ 
  geom_point(data = gamma_data, aes(grazing_intencity, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8,
             position=position_jitter(w=0.2))+
  scale_color_manual(values = col) # +  labs(color='Habitat type') + 

Fig.gammaENSPIE_grazing


## Mowing ----

plot_model(m1_3_ENSPIE,type = "pred", terms="mowing",  show.data=T,
           title = "", line.size=0.5)

Fig.gammaENSPIE_mowing <-ggplot(gamma_data, aes(mowing, gamma_100_ENSPIE)) + 
  geom_boxplot(color="grey")+
  labs(y=expression(paste("gamma ENS"[PIE])), x='Mowing')+ 
  geom_point(aes(color = habitat), position=position_jitter(w=0.1), size=3, alpha=0.8)+
  scale_color_manual(values = col)+
  labs(color='Habitat type')

Fig.gammaENSPIE_mowing 


## Prec_Varieb ----
Fig.gammaENSPIE_precip.CV <- plot_model(m2_1_ENSPIE,type = "pred", terms="Prec_Varieb", #  show.data=F
                                        title = "", line.size=1) + aes(linetype="solid") +
  labs(y=expression(paste("gamma ENS"[PIE])), x='Precipitation variability')+ 
  geom_point(data = gamma_data, aes(Prec_Varieb, gamma_100_ENSPIE, color = habitat), size=3, alpha=0.8)+
  scale_color_manual(values = col) # +  labs(color='Habitat type') 

Fig.gammaENSPIE_precip.CV 



# Plots combined ----

# for Gamma SR (100 m^2)

Fig.gammaENSPIE_clima + 
  Fig.gammaENSPIE_precip.CV +
  Fig.gammaENSPIE_soilC + 
  Fig.gammaENSPIE_soil.pH + 
  Fig.gammaENSPIE_Litter + 
  Fig.gammaENSPIE_grazing +
  plot_layout(ncol=3)


# combine with alpha

Fig.alphaENSPIE_clima + Fig.gammaENSPIE_clima +  
Fig.alphaENSPIE_precip.CV + Fig.gammaENSPIE_precip.CV +
Fig.alphaENSPIE_soilC + Fig.gammaENSPIE_soilC + 
Fig.alphaENSPIE_soil.pH + Fig.gammaENSPIE_soil.pH + 
Fig.alphaENSPIE_Litter + Fig.gammaENSPIE_Litter + 
Fig.alphaENSPIE_grazing + Fig.gammaENSPIE_grazing + 
  plot_layout(ncol=2)+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face = 'bold'), plot.tag.position = c(0.22, 0.95),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "pt"))



Fig.alphaSR_clima + Fig.gammaSR_clima +
  Fig.alphaSR_precip.CV + Fig.gammaSR_precip.CV +
  Fig.alphaSR_soilC + Fig.gammaSR_soilC +
  Fig.alphaSR_soil.pH + Fig.gammaSR_soil.pH +
  Fig.alphaSR_Litter + Fig.gammaSR_Litter +
  Fig.alphaSR_grazing + Fig.gammaSR_grazing +
  plot_layout(ncol=2)+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face = 'bold'), plot.tag.position = c(0.22, 0.95),
        plot.margin = unit(c(0,0,0,0), "pt"))



  

# R2 for the entire model---------

Anova(m1_3_ENSPIE)
Anova(m2_1_ENSPIE)# for Prec_Varieb


# write.csv(Anova(m1_3_ENSPIE),  file = "results/glmer_gamma_ENSPIE.csv")
# write.csv(coef(summary(m1_3_ENSPIE)),  file = "results/summary_gamma_ENSPIE.csv")
# write.csv(Anova(m2_1_ENSPIE),  file = "results/glmer_gamma_ENSPIE_2.csv")
# write.csv(coef(summary(m2_1_ENSPIE)),  file = "results/summary_gamma_ENSPIE_2.csv")


# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1_3_ENSPIE)
MuMIn::r.squaredGLMM(m2_1_ENSPIE)


# write.csv(MuMIn::r.squaredGLMM(m1_3_ENSPIE),  file = "results/Mod1_R2_gamma_ENSPIE.csv")
# write.csv(MuMIn::r.squaredGLMM(m2_1_ENSPIE),  file = "results/Mod2_R2_gamma_ENSPIE.csv")


# Partial R2 for fixed effects

R1_ENSPIE <- r2glmm::r2beta(m1_3_ENSPIE,  partial = T, data = gamma_data, method = 'sgv')
R1_ENSPIE

R2_ENSPIE <- r2glmm::r2beta(m2_1_ENSPIE,  partial = T, data = gamma_data, method = 'sgv')
R2_ENSPIE



R_ENSPIE <- R2_ENSPIE %>%
  filter(Effect=="poly(Prec_Varieb, 2)1" | Effect=="poly(Prec_Varieb, 2)2" | 
           Effect=="Model") %>% 
  bind_rows(R1_ENSPIE %>% filter(!Effect=="Model")) 

write.csv(R_ENSPIE,  file = "results/R2_gamma_ENSPIE.csv")




