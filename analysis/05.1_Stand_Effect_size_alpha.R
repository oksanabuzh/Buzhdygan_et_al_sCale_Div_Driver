# Calculate standardized effect size for each predictor on
## alpha species richness and alpha ENSPIE

# Steps:
# 1) Create composites to captures the collective effect of x and x^2 on y
# method: https://jslefche.github.io/sem_book/composite-variables.html#what-is-a-composite-variable
#         https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/atoms/files/SEM_10_2_Composites_for_Endogenous_Nonlinearities.pdf
# 2) Use the composites to predict y
# 3) Obtain the standardized coefficients of predictors
#         Standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model


dev.off()

# libraries----
library(tidyverse)
library(performance)
library(car)
library(lme4)
library(lmerTest)
# library(ggplot2)
library(sjPlot)
library(piecewiseSEM)


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
names(header)

# prepare subset of data for alpha scale (10-m2 plots)

alpha <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
    by = c("dataset", "plotID", "series", "subplot")
  ) %>%
  mutate(dataset = factor(dataset))

str(alpha)
names(alpha)

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
    Prec_Varieb = BIO15) %>%
  mutate(habitat = fct_relevel(habitat_broad,
    c("saline", "complex", "dry",
      "wet", "mesic", "fringe", "alpine"))) %>%

  drop_na

str(alpha_data)
alpha_data$dataset

# alpha SR -----


## 1) Create composites to captures the collective effect of x and x^2 on y ----

alpha_data$pca1_clima_2 <- (as.vector(scale(alpha_data$pca1_clima, center = T, scale = F)))^2

m1 <- glmer(alpha_10_div ~ pca1_clima + pca1_clima_2 +
  (1 | dataset / series), family = "poisson", data = alpha_data)


alpha_data$pH_2 <- (as.vector(scale(alpha_data$pH, center = T, scale = F)))^2

m2 <- glmer(alpha_10_div ~ pH + pH_2 +
  (1 | dataset / series), family = "poisson", data = alpha_data)


Corg_percent_2 <- (as.vector(scale(alpha_data$Corg_percent, center = T, scale = F)))^2

m3 <- glmer(alpha_10_div ~ Corg_percent + Corg_percent_2 +
  (1 | dataset / series), family = "poisson", data = alpha_data)
check_convergence(m3)


alpha_data$cover_litter_10 <- alpha_data$cover_litter / 10
alpha_data$cover_litter_2 <- (as.vector(scale(alpha_data$cover_litter_10, center = T, scale = T)))^2

m4 <- glmer(alpha_10_div ~ cover_litter_10 + cover_litter_2 +
  (1 | dataset / series), family = "poisson", data = alpha_data)

check_convergence(m4)
summary(m4)


alpha_data$Prec_Varieb_10 <- alpha_data$Prec_Varieb / 10
alpha_data$Prec_Varieb_2 <- (as.vector(scale(alpha_data$Prec_Varieb_10, center = T, scale = F)))^2

m5 <- glmer(alpha_10_div ~ Prec_Varieb_10 + Prec_Varieb_2 +
  (1 | dataset / series), family = "poisson", data = alpha_data)

# extract the coefficients, use them to generate the factor scores for the composits
alpha_dt <- alpha_data %>%
  mutate(clima_Comp = summary(m1)$coefficients[2, 1] * pca1_clima +
    summary(m1)$coefficients[3, 1] * pca1_clima_2) %>%
  mutate(pH_Comp = summary(m2)$coefficients[2, 1] * pH +
    summary(m2)$coefficients[3, 1] * pH_2) %>%
  mutate(Corg_Comp = summary(m3)$coefficients[2, 1] * Corg_percent +
    summary(m3)$coefficients[3, 1] * Corg_percent_2) %>%
  mutate(Litter_Comp = summary(m4)$coefficients[2, 1] * cover_litter_10 +
    summary(m4)$coefficients[3, 1] * cover_litter_2) %>%
  mutate(PrecipCV_Comp = summary(m5)$coefficients[2, 1] * Prec_Varieb_10 +
    summary(m5)$coefficients[3, 1] * Prec_Varieb_2)

alpha_dt

## 2) Use the composites to predict y ----
# use the final model (from the selection procedure)

m_alpha <- glmer(alpha_10_div ~
  clima_Comp +
  Corg_Comp +
  pH +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset / series), family = "poisson", data = alpha_dt)

check_convergence(m_alpha)
Anova(m_alpha)
check_collinearity(m_alpha)

## 3) Obtain  the standardized coefficients of predictors ----
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
alpha.SR_Std.Estimate <- coefs(m_alpha, standardize = "scale", standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
alpha.SR_Std.Estimate

### 3.1 Precipitation CV  ----

m_alpha <- glmer(alpha_10_div ~
  clima_Comp +
  PrecipCV_Comp +
  Corg_Comp +
  pH +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset / series), family = "poisson", data = alpha_dt)

Anova(m_alpha)
alpha.SR_Std.Estimate_PrecipCV <- coefs(m_alpha, standardize = "scale", standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
alpha.SR_Std.Estimate_PrecipCV





# alpha ENSPIE -----


# Check the model:

## 1) Create composites to captures the collective effect of x and x^2 on y ----

m1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ pca1_clima + pca1_clima_2 +
  (1 | dataset / series), data = alpha_data)

m2_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ pH + pH_2 +
  (1 | dataset / series), data = alpha_data)

m3_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ Corg_percent + Corg_percent_2 +
  (1 | dataset / series), data = alpha_data)

m4_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ cover_litter_10 + cover_litter_2 +
  (1 | dataset / series), data = alpha_data)


m5_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ Prec_Varieb_10 + Prec_Varieb_2 +
  (1 | dataset / series), data = alpha_data)


# extract the coefficients, use them to generate the factor scores for the composits
alpha_dt_ENSPIE <- alpha_data %>%
  mutate(clima_Comp = summary(m1_ENSPIE)$coefficients[2, 1] * pca1_clima +
    summary(m1_ENSPIE)$coefficients[3, 1] * pca1_clima_2) %>%
  mutate(pH_Comp = summary(m2_ENSPIE)$coefficients[2, 1] * pH +
    summary(m2_ENSPIE)$coefficients[3, 1] * pH_2) %>%
  mutate(Corg_Comp = summary(m3_ENSPIE)$coefficients[2, 1] * Corg_percent +
    summary(m3_ENSPIE)$coefficients[3, 1] * Corg_percent_2) %>%
  mutate(Litter_Comp = summary(m4_ENSPIE)$coefficients[2, 1] * cover_litter_10 +
    summary(m4_ENSPIE)$coefficients[3, 1] * cover_litter_2) %>%
  mutate(PrecipCV_Comp = summary(m5_ENSPIE)$coefficients[2, 1] * Prec_Varieb_10 +
    summary(m5_ENSPIE)$coefficients[3, 1] * Prec_Varieb_2)

## 2) Use the composites to predict y ----
# use the final model (from the selection procedure)

m_alpha_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~
  clima_Comp +
  pH +
  Corg_percent +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_dt_ENSPIE)

check_convergence(m_alpha_ENSPIE)
Anova(m_alpha_ENSPIE)
check_collinearity(m_alpha_ENSPIE)

## 3) Obtain  the standardized coefficients of predictors ----
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
alpha.ENSPIE_Std.Estimate <- coefs(m_alpha_ENSPIE, standardize = "scale", standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
alpha.ENSPIE_Std.Estimate

### 3.1 Precipitation CV  ----

m_alpha_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~
  clima_Comp +
  PrecipCV_Comp +
  pH +
  Corg_percent +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_dt_ENSPIE)

alpha.ENSPIE_Std.Estimate_PrecipCV <- coefs(m_alpha_ENSPIE, standardize = "scale", standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
alpha.ENSPIE_Std.Estimate_PrecipCV


# rebind tables----

alpha_st.eff <- rbind(alpha.SR_Std.Estimate,
  alpha.SR_Std.Estimate_PrecipCV[2, ],
  alpha.ENSPIE_Std.Estimate,
  alpha.ENSPIE_Std.Estimate_PrecipCV[2, ])

alpha_st.eff

write.csv(alpha_st.eff, "results/alpha_st.eff.csv", row.names = TRUE)
