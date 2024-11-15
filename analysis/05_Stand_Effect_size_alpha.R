# Purpose: Calculate standardized effect size for each predictor on alpha species richness

# Steps:
# 1) Create composites to captures the collective effect of x and x^2 on y 
# method: https://jslefche.github.io/sem_book/composite-variables.html#what-is-a-composite-variable
#         https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/atoms/files/SEM_10_2_Composites_for_Endogenous_Nonlinearities.pdf
# 2) Use the composites to predict y
# 3) Obtain the standardized coefficients of predictors
#         Standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model


dev.off()

#
library(tidyverse)

# data----
alpha <-read_csv ("Data/alpha_GLM.csv") 
str(alpha) 
names (alpha)

# dataset is a separate vegetation survey campaign
alpha$dataset <- factor(alpha$dataset)


# Remove NAs  
alpha_data <- alpha %>% 
  select(alpha_10_div, alpha_10_ENSPIE,  
         pca1_clima, 
         grazing_intencity, mowing, 
         # cover_shrub_total,     inclination, 
         cover_litter,
         BIO7, BIO15,
         pH, Corg_percent,
         dataset, series, eunis_group, habitat_broad, zonality) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

str(alpha_data)

alpha_data$dataset



alpha_data$pca1_clima_2<- (as.vector(scale(alpha_data$pca1_clima, center = TRUE, scale = F)))^2
alpha_data$pH_2<- (as.vector(scale(alpha_data$pH,  center = T, scale=F)))^2 
#alpha_$pH_2<-alpha_$pH^2 
cor(alpha_data$pH, alpha_data$pH_2)

alpha_data$cover_litter_2<- (as.vector(scale(alpha_data$cover_litter,  center = T, scale=F)))^2 


# alpha SR -----


# Check the model:
m <- glmer (alpha_10_div ~ 
               poly(pca1_clima, 2) +
               poly(pH, 2) +
               poly(Corg_percent,2)+
               poly(cover_litter,2) +
               grazing_intencity + mowing +
               (1|dataset/series),  
             family = "poisson", data = alpha_data)

# check overdispersion
sum(residuals(m, type = "pearson")^2) / df.residual(m)

check_collinearity(m)


plot_model(m,type = "pred", terms = c("pca1_clima"), show.data=T)
plot_model(m,type = "pred", terms = c("pH"), show.data=T)
plot_model(m,type = "pred", terms = c("Corg_percent"), show.data=T)
plot_model(m,type = "pred", terms = c("cover_litter"), show.data=T)
plot_model(m,type = "pred", terms = c("grazing_intencity"), show.data=T)
plot_model(m,type = "pred", terms = c("mowing"), show.data=T)


## 1) Create composites to captures the collective effect of x and x^2 on y ----

m1 <- glmer (alpha_10_div ~ pca1_clima + I(pca1_clima^2) + 
               (1|dataset/series), family = "poisson", data = alpha_data)

m2 <- glmer (alpha_10_div ~ pH + I(pH^2) + 
               (1|dataset/series), family = "poisson", data = alpha_data)

m3 <- glmer (alpha_10_div ~ Corg_percent + I(Corg_percent^2) + 
               (1|dataset/series), family = "poisson", data = alpha_data)
check_convergence(m3)

m4 <- glmer (alpha_10_div ~ cover_litter + cover_litter_2 + 
               (1|dataset/series), family = "poisson", data = alpha_data)
check_convergence(m4)

# extract the coefficients, use them to generate the factor scores for the composits 
alpha_dt <-alpha_data %>% 
  mutate(clima_Comp = summary(m1)$coefficients[2, 1] * pca1_clima + 
                      summary(m1)$coefficients[3, 1] * pca1_clima^2) %>% 
  mutate(pH_Comp = summary(m2)$coefficients[2, 1] * pH + 
                   summary(m2)$coefficients[3, 1] * pH^2) %>% 
  mutate(Corg_Comp = summary(m3)$coefficients[2, 1] * Corg_percent + 
                     summary(m3)$coefficients[3, 1] * Corg_percent^2) %>%
  mutate(Litter_Comp = summary(m4)$coefficients[2, 1] * cover_litter + 
                       summary(m4)$coefficients[3, 1] * cover_litter_2)

## 2) Use the composites to predict y ----

m_alpha <- glmer (alpha_10_div ~ 
                    clima_Comp + 
                    pH_Comp + # pH+
                    Corg_Comp +
                    Litter_Comp +
                    grazing_intencity + mowing +
                    (1|dataset/series), family = "poisson", data = alpha_dt)

check_convergence(m_alpha)
Anova(m_alpha)
check_collinearity(m_alpha)

## 3) Obtain  the standardized coefficients of predictors ----
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
alpha.SR_Std.Estimate <- coefs(m_alpha, standardize = "scale", standardize.type = "Menard.OE")[,c(1,2,7,8,9)]
alpha.SR_Std.Estimate



# alpha ENSPIE -----


# Check the model:
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

plot_model(m_ENSPIE_b,type = "pred", terms = c("pca1_clima"), show.data=T)
plot_model(m_ENSPIE_b,type = "pred", terms = c("pH"), show.data=T)
plot_model(m_ENSPIE_b,type = "pred", terms = c("Corg_percent"), show.data=T)
plot_model(m_ENSPIE_b,type = "pred", terms = c("cover_litter"), show.data=T)
plot_model(m_ENSPIE_b,type = "pred", terms = c("grazing_intencity"), show.data=T)
plot_model(m_ENSPIE_b,type = "pred", terms = c("mowing"), show.data=T)


## 1) Create composites to captures the collective effect of x and x^2 on y ----

m1_ENSPIE <- lmer (log(alpha_10_ENSPIE) ~ pca1_clima + I(pca1_clima^2) + 
               (1|dataset/series),  data = alpha_data)

m2_ENSPIE <- lmer (log(alpha_10_ENSPIE) ~ pH + I(pH^2) + 
               (1|dataset/series), data = alpha_data)

m3_ENSPIE <- lmer (log(alpha_10_ENSPIE) ~ Corg_percent + I(Corg_percent^2) + 
               (1|dataset/series), data = alpha_data)

m4_ENSPIE <- lmer (log(alpha_10_ENSPIE) ~ cover_litter + cover_litter_2 + 
               (1|dataset/series), data = alpha_data)
check_convergence(m4_ENSPIE)

# extract the coefficients, use them to generate the factor scores for the composits 
alpha_dt_ENSPIE <-alpha_data %>% 
  mutate(clima_Comp = summary(m1_ENSPIE)$coefficients[2, 1] * pca1_clima + 
           summary(m1_ENSPIE)$coefficients[3, 1] * pca1_clima^2) %>% 
  mutate(pH_Comp = summary(m2_ENSPIE)$coefficients[2, 1] * pH + 
           summary(m2_ENSPIE)$coefficients[3, 1] * pH^2) %>% 
  mutate(Corg_Comp = summary(m3_ENSPIE)$coefficients[2, 1] * Corg_percent + 
           summary(m3_ENSPIE)$coefficients[3, 1] * Corg_percent^2) %>%
  mutate(Litter_Comp = summary(m4_ENSPIE)$coefficients[2, 1] * cover_litter + 
           summary(m4_ENSPIE)$coefficients[3, 1] * cover_litter_2)

## 2) Use the composites to predict y ----

m_alpha_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                    clima_Comp + 
                    pH_Comp + # pH+
                    Corg_Comp + # Corg_percent +
                    Litter_Comp +
                    grazing_intencity + mowing +
                    (1|dataset/series), data = alpha_dt_ENSPIE)

check_convergence(m_alpha_ENSPIE)
Anova(m_alpha_ENSPIE)
check_collinearity(m_alpha_ENSPIE)

## 3) Obtain  the standardized coefficients of predictors ----
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
alpha.ENSPIE_Std.Estimate <- coefs(m_alpha_ENSPIE, standardize = "scale", standardize.type = "Menard.OE")[,c(1,2,7,8,9)]
alpha.ENSPIE_Std.Estimate
