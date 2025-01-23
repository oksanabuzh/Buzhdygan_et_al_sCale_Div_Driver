# Purpose: produce Figure 2

# Load libraries --------------------------------------------------------------
library(tidyverse)
library(sjPlot)
library(lme4)
library(lmerTest)
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

# Set theme for the model plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1,
          axis.textsize.y = 1,
          axis.textcolor = "black",
          axis.title.color = "black",
          axis.title.size = 1.4,
          legend.pos = "None",
          geom.linetype = 2)


# Read and prepare data -------------------------------------------------------
# This script prepared the following necessary data for all analyses
# - alpha_data: diversity, ENSPIE, cover and environmental variables for the 10m2 plots
# - gamma_data: diversity, ENSPIE, cover and environmental variables for the
#    100m2 plots
# - beta_data: diversity, ENSPIE, cover and environmental variables for the
#    beta scale
source("analysis/helper_scripts/prepare_data.R")

# Turn mowing variable into a factor
alpha_data <- mutate(alpha_data, mowing = factor(mowing))
beta_data <- mutate(beta_data, mowing = factor(mowing))
gamma_data <- mutate(gamma_data, mowing = factor(mowing))

# Check how the dataset looks like
alpha_data
beta_data
gamma_data

# 1) Get the predictions for the final selected models ------------------------

# alpha SR --------------------------------------------------------------------
m1_3 <- glmer(alpha_10_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  pH +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), family = "poisson", data = alpha_data)

m1_1 <- glmer(alpha_10_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), family = "poisson", data = alpha_data)

m2_1 <- glmer(alpha_10_div ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  pH +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), family = "poisson", data = alpha_data)

# Extract predictions from models ---------------------------------------------
clima_pred_10m <- get_model_data(m1_3, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m <- get_model_data(m1_3, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_10m <- get_model_data(m1_3, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_10m <- get_model_data(m1_1, type = "pred", terms = "pH[3.8:9, by=.001]")
grazing_pred_10m <- get_model_data(m1_3, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")


# alpha ENSPIE-----------------------------------------------------------------
m1_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~
  poly(pca1_clima, 2) +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

m2_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

# Extract predctions from models ----------------------------------------------

clima_pred_10m_Ensp <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m_Ensp <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_10m_Ensp <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_10m_Ensp <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pH[3.8:9, by=.001]")
grazing_pred_10m_Ensp <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m_Ensp <- get_model_data(m2_1_ENSPIE, type = "pred", terms = "Prec_Varieb")

# gamma SR --------------------------------------------------------------------
m1_1 <- glmer.nb(gamma_100_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

m2_1 <- glmer.nb(gamma_100_div ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

# Extract predictions from models ----------------------------------------------

clima_pred_100m <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_100m <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_100m <- get_model_data(m1_1, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_100m <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")

# gamma ENSPIE ----------------------------------------------------------------

m1_3_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~
  poly(pca1_clima, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

m2_1_ENSPIE <- lmer(log(gamma_100_ENSPIE) ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

# Extract predictions from models ----------------------------------------------

clima_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_100m_Ensp <- get_model_data(m1_3_ENSPIE, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m_Ensp <- get_model_data(m2_1_ENSPIE, type = "pred", terms = "Prec_Varieb")


# beta SR --------------------------------------------------------------------

m1_1 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)


m2_1 <- lmer(beta_100_div ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# Extract predictions from models ----------------------------------------------
clima_pred_beta <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_beta <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_beta <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_beta <- get_model_data(m1_1, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_beta <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_beta <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")

# beta ENSPIE ----------------------------------------------------------------

m1_1_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m2_1_ENSPIE <- lmer(beta_100_ENSPIE ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# Extract predictions from models ----------------------------------------------
clima_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_beta_ENSPIE <- get_model_data(m1_1_ENSPIE, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_beta_ENSPIE <- get_model_data(m2_1_ENSPIE, type = "pred", terms = "Prec_Varieb")


# 2) Plot the predictions  ----------------------------------------------------


# alpha-gamma--------------------------------------------------------------

# SR-----------------------------------------------------------------------

# Climate
Fig_Clima_10_100 <- ggplot(clima_pred_10m, aes(x, predicted)) +
  geom_point(data = gamma_data,
    aes(pca1_clima, gamma_100_div, fill = habitat, col = habitat), 
    size = 1, alpha = 0.8, pch = 8) +
  geom_ribbon(data = clima_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(data = clima_pred_10m, aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data,
    aes(pca1_clima, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Climate gradient') +
  geom_line(linetype = 1, linewidth = 1, col = "#50A0C8") +
  geom_line(data = clima_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")


## Soil C
Fig_Humus_10_100 <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(data = Humus_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(Corg_percent, gamma_100_div, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(Corg_percent, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Soil C') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8") +
  geom_line(data = Humus_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Litter % 
Fig_Litter_10_100 <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(data = Litter_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(cover_litter, gamma_100_div, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(cover_litter, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Litter cover') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8") +
  geom_line(data = Litter_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Soil pH
Fig_pH_10_100 <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(data = pH_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(pH, gamma_100_div, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(pH, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Soil pH') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8") +
  geom_line(data = pH_pred_100m, linetype = 1, linewidth = 0.5, col = "#D6604D")


# Grazing
Fig_grazing_10_100 <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(data = grazing_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(grazing_intencity, gamma_100_div, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8, position = position_jitter(width = 0.2)) +
  geom_point(data = alpha_data,
    aes(grazing_intencity, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Grazing intencity') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8") +
  geom_line(data = grazing_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")


# Prec_Varieb
Fig_precipCV_10_100 <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(data = precipCV_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(Prec_Varieb, gamma_100_div, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(Prec_Varieb, alpha_10_div, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Precipitation CV') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8") +
  geom_line(data = precipCV_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")


# ENSPIE ---------------------------------------------------------------------

# Climate
Fig_Clima_10_100_ENSPIE <- ggplot(clima_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = clima_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(pca1_clima, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(pca1_clima, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Climate gradient') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8") +
  geom_line(data = clima_pred_100m_Ensp, linetype = 1, linewidth = 1, col = "#D6604D")


# Soil C
Fig_Humus_10_100_ENSPIE <- ggplot(Humus_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = Humus_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(Corg_percent, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(Corg_percent, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil C') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8") +
  geom_line(data = Humus_pred_100m_Ensp, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Litter %
Fig_Litter_10_100_ENSPIE <- ggplot(Litter_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = Litter_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(cover_litter, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(cover_litter, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Litter cover') +
  geom_line(linetype = 1, size = 0.5, col = "#50A0C8") +
  geom_line(data = Litter_pred_100m_Ensp, linetype = 1, linewidth = 1, col = "#D6604D")

# Soil pH
Fig_pH_10_100_ENSPIE <- ggplot(pH_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = pH_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(pH, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(pH, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil pH') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8") +
  geom_line(data = pH_pred_100m_Ensp, linetype = 5, linewidth = 0.5, col = "#D6604D")


# Grazing
Fig_grazing_10_100_ENSPIE <- ggplot(grazing_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = grazing_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(grazing_intencity, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8, position = position_jitter(width = 0.2)) +
  geom_point(data = alpha_data,
    aes(grazing_intencity, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Grazing intencity') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8") +
  geom_line(data = grazing_pred_100m_Ensp, linetype = 1, linewidth = 0.5, col = "#D6604D")


# Prec_Varieb
Fig_precipCV_10_100_ENSPIE <- ggplot(precipCV_pred_10m_Ensp, aes(x, predicted)) +
  geom_ribbon(data = precipCV_pred_100m_Ensp, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = gamma_data,
    aes(Prec_Varieb, gamma_100_ENSPIE, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 8) +
  geom_point(data = alpha_data,
    aes(Prec_Varieb, alpha_10_ENSPIE, fill = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Precipitation CV') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8") +
  geom_line(data = precipCV_pred_100m_Ensp, linetype = 1, linewidth = 1, col = "#D6604D")

# beta scale -----------------------------------------------------------------

# SR -------------------------------------------------------------------------
# Climate

Fig.betaSR_clima <- ggplot(clima_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(pca1_clima, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Climate gradient') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")

# Soil C
Fig.betaSR_soilC <- ggplot(Humus_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Corg_percent, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Soil C') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")

# Litter %
Fig.betaSR_Litter <- ggplot(Litter_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(cover_litter, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Litter cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00966F")

# Soil pH
Fig.betaSR_soil.pH <- ggplot(pH_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(pH, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Soil pH') +
  geom_line(linetype = 1, linewidth = 1, col = "#00966F")

# Grazing
Fig.betaSR_grazing <- ggplot(grazing_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(grazing_intencity, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24,
    position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Grazing intencity') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")


# Prec_Varieb
Fig.betaSR_precip.CV <- ggplot(precipCV_pred_beta, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_div, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Precipitation CV') +
  geom_line(linetype = 1, linewidth = 1, col = "#00966F")

# ENSPIE ---------------------------------------------------------------------

# Clima
Fig.betaENSPIE_clima <- ggplot(clima_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(pca1_clima, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Climate gradient') +
  geom_line(linetype = 1, linewidth = 0.5, col = "#00966F")

# Soil C
Fig.betaENSPIE_soilC <- ggplot(Humus_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Corg_percent, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil C') +
  geom_line(linetype = 1, size = 1, col = "#00966F")

# Litter %
Fig.betaENSPIE_Litter <- ggplot(Litter_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(cover_litter, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Litter cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00966F")

# Soil pH
Fig.betaENSPIE_soil.pH <- ggplot(pH_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(pH, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Soil pH') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")

# Grazing
Fig.betaENSPIE_grazing <- ggplot(grazing_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(grazing_intencity, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24,
    position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Grazing intencity') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")

# Prec_Varieb
Fig.betaENSPIE_precip.CV <- ggplot(precipCV_pred_beta_ENSPIE, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(Prec_Varieb, beta_100_ENSPIE, fill = habitat, col = habitat), size = 0.7, alpha = 0.8, pch = 24) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = expression(paste("ENS"[PIE])), x = 'Precipitation CV') +
  geom_line(linetype = 5, linewidth = 0.5, col = "#00966F")

# Combine all plots -----------------------------------------------------------

Fig2 <-
  Fig_Clima_10_100 + Fig_Clima_10_100_ENSPIE + Fig.betaSR_clima + Fig.betaENSPIE_clima +
    Fig_precipCV_10_100 + Fig_precipCV_10_100_ENSPIE + Fig.betaSR_precip.CV + Fig.betaENSPIE_precip.CV +
    Fig_Humus_10_100 + Fig_Humus_10_100_ENSPIE + Fig.betaSR_soilC + Fig.betaENSPIE_soilC +
    Fig_pH_10_100 + Fig_pH_10_100_ENSPIE + Fig.betaSR_soil.pH + Fig.betaENSPIE_soil.pH +
    Fig_Litter_10_100 + Fig_Litter_10_100_ENSPIE + Fig.betaSR_Litter + Fig.betaENSPIE_Litter +
    Fig_grazing_10_100 + Fig_grazing_10_100_ENSPIE + Fig.betaSR_grazing + Fig.betaENSPIE_grazing +
    plot_annotation(tag_levels = 'a') +
    plot_layout(ncol = 4) & ylab(NULL) & 
    theme(plot.margin = margin(3, 1, 3, 20),
      plot.tag = element_text(size = 6, face = 'bold'),
      plot.tag.position = c(0.15, 1.06)
      )

Fig2
