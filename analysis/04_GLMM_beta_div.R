# Purpose: GLMM analysis for beta diversity


# load libraries --------------------------------------------
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(sjPlot)
library(performance)
library(patchwork)

# Define habitat colors to be used in all plots to distinguish habitats
habitat_colors = c(
  saline = "#4e3910",
  complex = "#CC6600",
  dry = "#e3c28b",
  wet = "#CC99FF",
  mesic = "#0066FF",
  fringe = "#00B200",
  alpine = "#006600")

# Read and prepare data --------------------------------------------

# SR - species richness
# ENSPIE - evenness measure calculated as inverse Simpson using species cover
# cover - is cumulative plant cover

# "data/alpha_beta_gamma_community_variabl.csv" combines all diversity measures and plant cover
# alpha diversity measures (SR and ENSPIE) include doubled 10 m2 plots,
## thus "series" (i.e. 100m2 plots), nested in dataset (separate vegetation survey campaign)
## are fitted as a random effect
# gamma diversity measures (SR and ENSPIE)include 100m2 plots (i.e. the sample size is half of what we have for the 10m2 plots)
# beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

# Read climate data and compund climate variable from PCA analysis in "1_prepare_data/ PCA_environment.R"
climate_PCA <- read.csv("data/climate_PCA.csv")

# Read all environmental data
header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

# Calulate mean header data to combine two corners within each plot into one value
header_mean <- header %>%
  select(c(series, zonality, habitat_broad,
    where(is.numeric))) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# Prepare subset of data for beta scale -------------------------

beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series")) %>%
  mutate(dataset = factor(dataset))

# Remove NAs and select only needed variables
beta_data <- beta_gamma %>%
  dplyr::select(beta_100_div, beta_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
    cover_litter,
    Tem_range, Prec_Varieb,
    pH, Corg_percent,
    dataset, series, habitat_broad, zonality) %>%
  mutate(mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()


# Start Analysis -------------------------------------------------------------

#-----------------------------------------------------------------------------#
# (1) Species richness -------------------------------------------------------
# ----------------------------------------------------------------------------#

# GLLM analyses ----------------------------------------------------------------

# Exploration ------------------------------------------------------------------

m <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

check_convergence(m)

# check model assumptions
plot(m)
qqnorm(resid(m))
qqline(resid(m))

# check multicolinearity
check_collinearity(m)

Anova(m)
summary(m)

# Model selection ---------------------------------------------------------------

## Model 1: all predictors (except precipitation CV) ----------------------------
# test quadratic effects of climate, soil C, pH, and litter

# poly(pca1_clima, 2) nonsignificant
# poly(Corg_percent,2)nonsignificant

m1_1 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_2 <- lmer(beta_100_div ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_3 <- lmer(beta_100_div ~
  pca1_clima +
  Corg_percent +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_4 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  Corg_percent +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m1_1, m1_2, m1_3, m1_4) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

Anova(m1_1)

# Model 2: Add precipitation variablity ---------------------------------------
m2_1 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m2_2 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m2_1, m2_2) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

# Check the best model
Anova(m2_1)
# Anova(m2_2)

# Additional analysis for precipitation CV effects  -----------------------------
# Analysis whether precipitation variability adds explanatory power beyond the
# nonlinear effect of the climate gradient:
# Model m3_1 includes both linear and quadratic terms for the climate gradient.
# Model m3_2 includes the linear effects of both climate gradient and
# precipitation variability (i.e., the quadratic term for the climate gradient was replaced by precipitation variability).
# If the latter model is better (e.g. AIC is smaller than two units), we have more support for claim that the precipitation variability effect is shown.

m3_1 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  (1 | dataset), data = beta_data)

m3_2 <- lmer(beta_100_div ~
  pca1_clima + Prec_Varieb +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m3_1, m3_2) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

Anova(m3_1)
Anova(m3_2)

# Make model plots ------------------------------------------------------------

# final models:
Anova(m1_1)
Anova(m2_1)

# Climate plots ----------------------------------------------------------------

clima_pred_beta <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.betaSR_clima <- ggplot(clima_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(pca1_clima, beta_100_div, col = habitat),
    size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Climate gradient') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaSR_clima

# Soil C plot -----------------------------------------------------------------

Humus_pred_beta <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")

Fig.betaSR_soilC <- ggplot(Humus_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(Corg_percent, beta_100_div, col = habitat),
    size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Soil C') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaSR_soilC

# Litter % plot --------------------------------------------------------------

Litter_pred_beta <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")

Fig.betaSR_Litter <- ggplot(Litter_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(cover_litter, beta_100_div, col = habitat),
    size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Litter cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaSR_Litter

# Soil pH plot ----------------------------------------------------------------
pH_pred_beta <- get_model_data(m1_1, type = "pred", terms = "pH[3.7:9, by=.001]")

Fig.betaSR_soil.pH <- ggplot(pH_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(pH, beta_100_div, col = habitat),
    size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Soil pH') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaSR_soil.pH

# Grazing plot ---------------------------------------------------------------
grazing_pred_beta <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")

Fig.betaSR_grazing <- ggplot(grazing_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(grazing_intencity, beta_100_div, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Grazing intencity') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaSR_grazing

# Mowing plot ----------------------------------------------------------------

Fig.betaSR_mowing <- ggplot(beta_data, aes(mowing, beta_100_div)) +
  geom_boxplot(color = "#00AC7F") +
  labs(y = "beta SR", x = 'Mowing') +
  geom_point(aes(color = habitat), pch = 19, 
             position = position_jitter(width = 0.1), size = 3, alpha = 0.8) +
  scale_color_manual(values = habitat_colors)

Fig.betaSR_mowing

# Precipitation variability plot ----------------------------------------------

precipCV_pred_beta <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")

Fig.betaSR_precip.CV <- ggplot(precipCV_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_div, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = "beta SR", x = 'Precipitation variability') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaSR_precip.CV

# Final statistics ------------------------------------------------------------
Anova(m1_1)
Anova(m2_1) # for Prec_Varieb

# R2 for the entire model -----------------------------------------------------
# R2m and R2c are marginal (for fixed predictors) and
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1_1)
MuMIn::r.squaredGLMM(m2_1)

# Partial R2 for fixed effects
R1_SR <- r2glmm::r2beta(m1_1, partial = T, data = beta_data, method = 'sgv')
R1_SR

R2_SR <- r2glmm::r2beta(m2_1, partial = T, data = beta_data, method = 'sgv')
R2_SR

R_SR <- R2_SR %>%
  filter(Effect == "poly(Prec_Varieb, 2)1" | Effect == "poly(Prec_Varieb, 2)2" |
    Effect == "Model") %>%
  bind_rows(R1_SR %>% filter(!Effect == "Model"))

write.csv(R_SR, file = "results/R2_beta_SR.csv")

#-----------------------------------------------------------------------------#
# (1) ENSPIE -------------------------------------------------------
# ----------------------------------------------------------------------------#

# GLLM ------------------------------------------------------------------------

# Exploration -----------------------------------------------------------------

# Check the model:
m_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# check model
plot(m_ENSPIE) # heteroscedasticity
qqnorm(resid(m_ENSPIE))
qqline(resid(m_ENSPIE))

Anova(m_ENSPIE)

# Model 1: all predictors (except precipitation CV) ---------------------------
# test quadratic effects of climate, soil C, pH, and litter

# poly(pca1_clima, 2)
# poly(pH, 2) is marginal

m1_1_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_2_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  pH +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_3_ENSPIE <- lmer(beta_100_ENSPIE ~
  pca1_clima +
  pH +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_4_ENSPIE <- lmer(beta_100_ENSPIE ~
  pca1_clima +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m1_1_ENSPIE, m1_2_ENSPIE, m1_3_ENSPIE, m1_4_ENSPIE) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

# Check the best model
Anova(m1_1_ENSPIE)
# Anova(m1_4_ENSPIE)

## Model 2: Add precipitation variability -------------------------------------

m2_1_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m2_2_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  Prec_Varieb +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m2_1_ENSPIE, m2_2_ENSPIE) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

# Check the best model
Anova(m2_1_ENSPIE)

# Additional analysis for precipitation CV effects  -----------------------------
# Analysis whether precipitation variability adds explanatory power beyond the
# nonlinear effect of the climate gradient:
# Model m3_1_ENSPIE includes both linear and quadratic terms for the climate gradient.
# Model m3_2_ENSPIE includes the linear effects of both climate gradient and
# precipitation variability (i.e., the quadratic term for the climate gradient was replaced by precipitation variability).
# If the latter model is better (e.g. AIC is smaller than two units), we have more support for claim that the precipitation variability effect is shown.

m3_1_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  (1 | dataset), data = beta_data)

m3_2_ENSPIE <- lmer(beta_100_ENSPIE ~
  pca1_clima +
  Prec_Varieb +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m3_1_ENSPIE, m3_2_ENSPIE) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

Anova(m3_1_ENSPIE)
Anova(m3_2_ENSPIE)

# Make model plots ------------------------------------------------------------

# final models:
Anova(m1_1_ENSPIE)
Anova(m2_1_ENSPIE)

# Climate plot ----------------------------------------------------------------

clima_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.betaENSPIE_clima <- ggplot(clima_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data,
    aes(pca1_clima, beta_100_ENSPIE, col = habitat),
    size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Climate gradient') +
  geom_line(linetype = 1, linewidth = 0.5, col = "#00AC7F")

Fig.betaENSPIE_clima

# Soil C plot -----------------------------------------------------------------

Humus_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")

Fig.betaENSPIE_soilC <- ggplot(Humus_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Corg_percent, beta_100_ENSPIE, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil C') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaENSPIE_soilC

# Litter % plot --------------------------------------------------------------

Litter_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "cover_litter[0:100, by=0.01]")

Fig.betaENSPIE_Litter <- ggplot(Litter_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(cover_litter, beta_100_ENSPIE, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Litter cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaENSPIE_Litter

# Soil pH plot ----------------------------------------------------------------
pH_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pH[3.7:9, by=.001]")

Fig.betaENSPIE_soil.pH <- ggplot(pH_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(pH, beta_100_ENSPIE, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil pH') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaENSPIE_soil.pH

# Grazing plot ---------------------------------------------------------------
grazing_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")

Fig.betaENSPIE_grazing <- ggplot(grazing_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(grazing_intencity, beta_100_ENSPIE, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Grazing intencity') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaENSPIE_grazing

# Mowing plot ----------------------------------------------------------------

Fig.betaENSPIE_mowing <- ggplot(beta_data, aes(mowing, beta_100_ENSPIE)) +
  geom_boxplot(color = "#00AC7F") +
  labs(y = expression(paste("ENS"[PIE])), x = 'Mowing') +
  geom_point(aes(color = habitat), pch = 19, position = position_jitter(width = 0.1), size = 3, alpha = 0.8) +
  scale_color_manual(values = habitat_colors)

Fig.betaENSPIE_mowing

# Precipitation variability plot ----------------------------------------------
precipCV_pred_beta_ENSPIE <- get_model_data(m2_1_ENSPIE, type = "pred", terms = "Prec_Varieb")

Fig.betaENSPIE_precip.CV <- ggplot(precipCV_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_ENSPIE, col = habitat), size = 1, alpha = 0.8, pch = 19) +
  scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Precipitation CV') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00AC7F")

Fig.betaENSPIE_precip.CV
# final statistics -----------------------------------------------------------
Anova(m1_1_ENSPIE)
Anova(m2_1_ENSPIE)

# R2 for the entire model -----------------------------------------------------
MuMIn::r.squaredGLMM(m1_1_ENSPIE)
MuMIn::r.squaredGLMM(m2_1_ENSPIE)

R1_ENSPIE <- r2glmm::r2beta(m1_1_ENSPIE, partial = T, data = beta_data, method = 'sgv')
R1_ENSPIE

R2_ENSPIE <- r2glmm::r2beta(m2_1_ENSPIE, partial = T, data = beta_data, method = 'sgv')
R2_ENSPIE

R_ENSPIE <- R2_ENSPIE %>%
  filter(Effect == "poly(Prec_Varieb, 2)1" | Effect == "poly(Prec_Varieb, 2)2" |
    Effect == "Model") %>%
  bind_rows(R1_ENSPIE %>% filter(!Effect == "Model"))

write.csv(R_ENSPIE, file = "results/R2_beta_ENSPIE.csv")
