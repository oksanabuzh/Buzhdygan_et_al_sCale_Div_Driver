# Purpose: Make Figure S8

# Load libraries -------------------------------------------------------------
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(sjPlot)
library(performance)
library(patchwork)


# Define habitat colors to be used in all plots to distinguish habitats
habitat_colors <- c(
  saline = "#4e3910",
  complex = "#CC6600",
  dry = "#e3c28b",
  wet = "#CC99FF",
  mesic = "#0066FF",
  fringe = "#00B200",
  alpine = "#006600"
)

# Set theme for the plots
set_theme(
  base = theme_bw(),
  axis.textsize.x = 0.8,
  axis.textsize.y = 0.8,
  axis.textcolor = "black",
  axis.title.color = "black",
  axis.title.size = 0.9,
  legend.pos = "None",
  geom.linetype = 2
)

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

# GLMM models ----------------------------------------------------------------

# alpha scale ----------------------------------------------------------------

m1_1 <- lmer(sqrt(alpha_10_cover) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

# check model
plot(m1_1) # heteroscedasticity
qqnorm(resid(m1_1))
qqline(resid(m1_1))

summary(m1_1)
Anova(m1_1)

m1_2 <- lmer(sqrt(alpha_10_cover) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

summary(m1_2)
Anova(m1_2)

plot_model(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)
plot_model(m1_2, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


# calculate AIC to compare models
AIC(m1_1, m1_2)

Anova(m1_1)

# Add precipiation variability -----------------------------------------------
# Select model
m2_1 <- lmer(sqrt(alpha_10_cover) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

m2_2 <- lmer(sqrt(alpha_10_cover) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data)

# calculate AIC
AIC(m2_1, m2_2)

Anova(m2_1)

# Extract predictions alpha scale ---------------------------------------------
clima_pred_10m <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_10m <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_10m <- get_model_data(m1_1, type = "pred", terms = "pH[3.8:9, by=.001]")
grazing_pred_10m <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")

# gamma models ---------------------------------------------------------------
m3_1 <- lmer(sqrt(gamma_100_cover) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

# check model
plot(m3_1) # heteroscedasticity
qqnorm(resid(m3_1))
qqline(resid(m3_1))

summary(m3_1)
Anova(m3_1)

m3_2 <- lmer(sqrt(gamma_100_cover) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

summary(m3_2)
Anova(m3_2)

m3_3 <- lmer(sqrt(gamma_100_cover) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  cover_litter +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

Anova(m3_3)

# calculate AIC
AIC(m3_1, m3_2, m3_3)

# Checkout the best model
Anova(m3_1)
plot_model(m3_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


# Add precipiation variability -----------------------------------------------
# select model
m4_1 <- lmer(sqrt(gamma_100_cover) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

m4_2 <- lmer(sqrt(gamma_100_cover) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data)

# calculate AIC
AIC(m4_1, m4_2)

# checkout the best model
Anova(m4_1)

# Extract predictions gamma scale ---------------------------------------------
clima_pred_100m <- get_model_data(m3_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m <- get_model_data(m3_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_100m <- get_model_data(m3_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_100m <- get_model_data(m3_1, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_100m <- get_model_data(m3_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m <- get_model_data(m4_1, type = "pred", terms = "Prec_Varieb")


# 2) Plot the predictions  ---------------------------------------------------

# alpha-gamma ---------------------------------------------------------------

# SR -------------------------------------------------------------------------

# Climate ---------------------------------------------------------------------
Fig_Clima_10 <- ggplot(clima_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(pca1_clima, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Climate gradient") +
geom_line(linetype = 5, linewidth = 0.5, col = "#50A0C8")

Fig_Clima_100 <- ggplot(clima_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = clima_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(pca1_clima, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Climate gradient") +
  geom_line(data = clima_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Soil C ----------------------------------------------------------------------
Fig_Humus_10 <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(Corg_percent, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil C") +
  geom_line(linetype = 1, linewidth = 1, col = "#50A0C8")

Fig_Humus_100 <- ggplot(Humus_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = Humus_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(Corg_percent, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil C") +
  geom_line(data = Humus_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Litter % --------------------------------------------------------------------
Fig_Litter_10 <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(cover_litter, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Litter cover") +
  geom_line(linetype = 1, linewidth = 1, col = "#50A0C8")

Fig_Litter_100 <- ggplot(Litter_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = Litter_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(cover_litter, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Litter cover") +
  geom_line(data = Litter_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Soil pH ---------------------------------------------------------------------
Fig_pH_10 <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(pH, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil pH") +
  geom_line(linetype = 1, linewidth = 1, col = "#50A0C8")

Fig_pH_100 <- ggplot(pH_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = pH_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(pH, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil pH") +
  geom_line(data = pH_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Grazing ---------------------------------------------------------------------
Fig_grazing_10 <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(grazing_intencity, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19, position = position_jitter(width = 0.2)
  ) +
  labs(y = "Plant cover", x = "Grazing intencity") +
  geom_line(linetype = 5, linewidth = 0.5, col = "#50A0C8")

Fig_grazing_100 <- ggplot(grazing_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = grazing_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(grazing_intencity, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19, position = position_jitter(width = 0.2)
  ) +
  labs(y = "Plant cover", x = "Grazing intencity") +
  geom_line(data = grazing_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Precipitation CV ------------------------------------------------------------
Fig_precipCV_10 <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data,
    aes(Prec_Varieb, alpha_10_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Precipitation CV") +
  geom_line(linetype = 5, linewidth = 0.5, col = "#50A0C8")

Fig_precipCV_100 <- ggplot(precipCV_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = precipCV_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data,
    aes(Prec_Varieb, gamma_100_cover, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Precipitation CV") +
  geom_line(data = precipCV_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Combine all plots ----------------------------------------------------------

FigS8 <-
  Fig_Clima_10 + Fig_Clima_100 +
    Fig_precipCV_10 + Fig_precipCV_100 +
    Fig_Humus_10 + Fig_Humus_100 +
    Fig_pH_10 + Fig_pH_100 +
    Fig_Litter_10 + Fig_Litter_100 +
    Fig_grazing_10 + Fig_grazing_100 +
    plot_annotation(tag_levels = "a") +
    plot_layout(ncol = 2) & 
    scale_color_manual(values = habitat_colors) &
    ylab(NULL) &
    theme(
      plot.margin = margin(3, 1, 3, 20),
      plot.tag = element_text(size = 6, face = "bold"),
      plot.tag.position = c(0.15, 1.06)
    )

FigS8
