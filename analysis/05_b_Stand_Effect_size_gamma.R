# Purpose: Calculate standardized effect size for each predictor on
## gamma species richness and alpha ENSPIE

# Steps:
# 1) Create composites to captures the collective effect of x and x^2 on y
# method: https://jslefche.github.io/sem_book/composite-variables.html#what-is-a-composite-variable
#         https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/atoms/files/SEM_10_2_Composites_for_Endogenous_Nonlinearities.pdf
# 2) Use the composites to predict y
# 3) Obtain the standardized coefficients of predictors
#         Standardized coefficient of composite is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model


# Load libraries -----------------------------------------------------------
library(tidyverse)
library(performance)
library(car)
library(lme4)
library(lmerTest)
library(piecewiseSEM)


# Read and prepare data -------------------------------------------------------
# This script prepared the following necessary data for all analyses
# - alpha_data: diversity, ENSPIE, cover and environmental variables for the 10m2 plots
# - gamma_data: diversity, ENSPIE, cover and environmental variables for the 
#    100m2 plots
# - beta_data: diversity, ENSPIE, cover and environmental variables for the
#    beta scale
source("analysis/helper_scripts/prepare_data.R")

# Check how the dataset looks like
gamma_data

# gamma SR -------------------------------------------------------------------

# 1) Create composites to captures the collective effect of x and x^2 on y -------------------------

# Climate
gamma_data$pca1_clima_2 <- (as.vector(scale(gamma_data$pca1_clima,
  center = TRUE, scale = FALSE)))^2

m1 <- glmer.nb(gamma_100_div ~ pca1_clima + pca1_clima_2 +
  (1 | dataset), data = gamma_data)

# pH
gamma_data$pH_2 <- (as.vector(scale(gamma_data$pH, center = TRUE, scale = FALSE)))^2

m2 <- glmer.nb(gamma_100_div ~ pH + pH_2 +
  (1 | dataset), data = gamma_data)
check_convergence(m2)

# C organic
Corg_percent_2 <- (as.vector(scale(gamma_data$Corg_percent, center = TRUE,
  scale = FALSE)))^2

m3 <- glmer.nb(gamma_100_div ~ Corg_percent + Corg_percent_2 +
  (1 | dataset), data = gamma_data)

# Litter cover
gamma_data$cover_litter_10 <- gamma_data$cover_litter / 10
gamma_data$cover_litter_2 <- (as.vector(scale(gamma_data$cover_litter_10,
  center = TRUE, scale = TRUE)))^2

m4 <- glmer.nb(gamma_100_div ~ cover_litter_10 + cover_litter_2 +
  (1 | dataset), data = gamma_data)

# Precipitation CV
gamma_data$Prec_Varieb_10 <- gamma_data$Prec_Varieb / 10
gamma_data$Prec_Varieb_2 <- (as.vector(scale(gamma_data$Prec_Varieb_10,
  center = TRUE, scale = FALSE)))^2

m5 <- glmer.nb(gamma_100_div ~ Prec_Varieb_10 + Prec_Varieb_2 +
  (1 | dataset), data = gamma_data)

# extract the coefficients, use them to generate the factor scores for the composits
gamma_dt <- gamma_data %>%
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

# 2) Use the composites to predict y -----------------------------------------
m_gamma <- glmer.nb(gamma_100_div ~
  clima_Comp +
  Corg_Comp +
  pH_Comp +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_dt)

check_convergence(m_gamma)

Anova(m_gamma)
check_collinearity(m_gamma)

# 3) Obtain  the standardized coefficients of predictors ---------------------
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
gamma.SR_Std.Estimate <- coefs(m_gamma, standardize = "scale",
  standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
gamma.SR_Std.Estimate


# 3.1 Precipitation CV  -----------------------------------------------------
m_gamma <- glmer.nb(gamma_100_div ~
  clima_Comp +
  PrecipCV_Comp +
  Corg_Comp +
  pH_Comp +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_dt)

Anova(m_gamma)

gamma.SR_Std.Estimate_PrecipCV <- coefs(m_gamma, standardize = "scale",
  standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
gamma.SR_Std.Estimate_PrecipCV

# gamma ENSPIE ----------------------------------------------------------------

# 1) Create composites to captures the collective effect of x and x^2 on y -------------------------

m1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ pca1_clima + pca1_clima_2 +
  (1 | dataset), data = gamma_data)

m2_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ pH + pH_2 +
  (1 | dataset), data = gamma_data)

m3_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ Corg_percent + Corg_percent_2 +
  (1 | dataset), data = gamma_data)

m4_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ cover_litter_10 + cover_litter_2 +
  (1 | dataset), data = gamma_data)

m5_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~ Prec_Varieb_10 + Prec_Varieb_2 +
  (1 | dataset), data = gamma_data)

# extract the coefficients, use them to generate the factor scores for the composits
gamma_dt_ENSPIE <- gamma_data %>%
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

# 2) Use the composites to predict y -----------------------------------------

m_gamma_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~
  clima_Comp +
  pH_Comp +
  Corg_Comp +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_dt_ENSPIE)

check_convergence(m_gamma_ENSPIE)
Anova(m_gamma_ENSPIE)
check_collinearity(m_gamma_ENSPIE)

# 3) Obtain  the standardized coefficients of predictors ---------------------
# standardized coefficient of composit is the total non-linear effect of predictor(i.e., x and x^2), controlling for other predictors in the model
# use the coefs function from piecewiseSEM to obtain the standardized coefficients
gamma.ENSPIE_Std.Estimate <- coefs(m_gamma_ENSPIE, standardize = "scale",
  standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
gamma.ENSPIE_Std.Estimate

# 3.1 Precipitation CV  -----------------------------------------------------
m_gamma_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~
  clima_Comp +
  PrecipCV_Comp +
  pH_Comp +
  Corg_Comp +
  Litter_Comp +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_dt_ENSPIE)

gamma.ENSPIE_Std.Estimate_PrecipCV <- coefs(m_gamma_ENSPIE, standardize = "scale",
  standardize.type = "Menard.OE")[, c(1, 2, 7, 8, 9)]
gamma.ENSPIE_Std.Estimate_PrecipCV

# Combine tables -----------------------------------------------------------
gamma_st.eff <- rbind(gamma.SR_Std.Estimate,
  gamma.SR_Std.Estimate_PrecipCV[2, ],
  gamma.ENSPIE_Std.Estimate,
  gamma.ENSPIE_Std.Estimate_PrecipCV[2, ])

gamma_st.eff

write.csv(gamma_st.eff, "results/gamma_st.eff.csv", row.names = TRUE)
