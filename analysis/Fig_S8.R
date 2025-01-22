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

# SR - species richness
# ENSPIE - evenness measure calculated as inverse Simpson using species cover
# cover - is cumulative plant cover

# "data/alpha_beta_gamma_community_variabl.csv" combines all diversity measures and plant cover
# alpha diversity measures (SR and ENSPIE) include doubled 10 m2 plots,
# thus "series" (i.e. 100m2 plots), nested in dataset (separate vegetation survey campaign)
# are fitted as a random effect
# gamma diversity measures (SR and ENSPIE)include 100m2 plots (i.e. the sample size is half of what we have for the 10m2 plots)
# beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

# Read climate data and compund climate variable from PCA analysis in "1_prepare_data/ PCA_environment.R"
climate_PCA <- read_csv("data/climate_PCA.csv")

# Read all environmental data
header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

# mean per series (per 100m2 plots)
header_mean <- header %>%
  select(c(
    series, zonality, habitat_broad,
    where(is.numeric)
  )) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# Prepare subset of data for alpha scale (10 m2 plots) -------------------------

alpha <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
    by = c("dataset", "plotID", "series", "subplot")
  ) %>%
  mutate(dataset = factor(dataset))

str(alpha)

# Remove NAs and select only needed variables
alpha_data <- alpha %>%
  dplyr::select(
    alpha_10_div, alpha_10_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
    dataset, series, habitat_broad,
    subplot
  ) %>%
  mutate(
    Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)
  ) %>%
  mutate(habitat = fct_relevel(habitat_broad, c(
    "saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"
  ))) %>%
  drop_na()

str(alpha_data)

# plot on a mean alpha per series to omit pseudoreplication of the plots:
alpha_mean <- alpha_data %>%
  dplyr::select(
    alpha_10_div, alpha_10_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
    dataset, series, habitat_broad
  ) %>%
  mutate(
    Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)
  ) %>%
  mutate(habitat = fct_relevel(habitat_broad, c(
    "saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"
  ))) %>%
  drop_na()

str(alpha_mean)

# Prepare subset of data for gamma scale (100 m2 plots) -------------------------
beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series")) %>%
  mutate(dataset = factor(dataset))

str(beta_gamma)

# selected variables, removed NAs
gamma_data <- beta_gamma %>%
  dplyr::select(
    gamma_100_div, gamma_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    BIO7, BIO15, BIO1, BIO12,
    pH, Corg_percent,
    dataset, series, habitat_broad, zonality
  ) %>%
  mutate(
    Tem_range = BIO7,
    Prec_Varieb = BIO15,
    Temprt = BIO1,
    Precipt = BIO12
  ) %>%
  mutate(habitat = fct_relevel(
    habitat_broad,
    c(
      "saline", "complex", "dry",
      "wet", "mesic", "fringe", "alpine"
    )
  )) %>%
  drop_na()

str(gamma_data)


# Plant cover data -----------------------------------------------------------

# todo: can this be simplified?
# alpha scale
tot_cover_10 <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(scale == 10 &
    metric == "cover") %>%
  rename(total_cover_10 = value) %>%
  select(dataset, series, subplot, total_cover_10) %>%
  mutate(dataset = factor(dataset))

alpha_data_cover <- alpha_data %>%
  left_join(tot_cover_10, by = c("dataset", "series", "subplot"))

# gamma scale
tot_cover_100 <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(scale == 100 &
    metric == "cover") %>%
  rename(total_cover_100 = value) %>%
  select(dataset, series, subplot, total_cover_100) %>%
  mutate(dataset = factor(dataset))

gamma_data_cover <- gamma_data %>%
  left_join(tot_cover_100, by = c("dataset", "series"))


# GLMM models ----------------------------------------------------------------

# alpha scale ----------------------------------------------------------------

m1_1 <- lmer(sqrt(total_cover_10) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

# check model
plot(m1_1) # heteroscedasticity
qqnorm(resid(m1_1))
qqline(resid(m1_1))

summary(m1_1)
Anova(m1_1)

m1_2 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

summary(m1_2)
Anova(m1_2)

plot_model(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)
plot_model(m1_2, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


# calculate AIC to compare models
AIC(m1_1, m1_2)

Anova(m1_1)

# Add precipiation variability -----------------------------------------------
# Select model
m2_1 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

m2_2 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

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
m3_1 <- lmer(sqrt(total_cover_100) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

# check model
plot(m3_1) # heteroscedasticity
qqnorm(resid(m3_1))
qqline(resid(m3_1))

summary(m3_1)
Anova(m3_1)

m3_2 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

summary(m3_2)
Anova(m3_2)

m3_3 <- lmer(sqrt(total_cover_100) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  cover_litter +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

Anova(m3_3)

# calculate AIC
AIC(m3_1, m3_2, m3_3)

# Checkout the best model
Anova(m3_1)
plot_model(m3_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


# Add precipiation variability -----------------------------------------------
# select model
m4_1 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

m4_2 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

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
    data = alpha_data_cover,
    aes(pca1_clima, total_cover_10, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Climate gradient")
geom_line(linetype = 5, linewidth = 0.5, col = "#50A0C8")

Fig_Clima_100 <- ggplot(clima_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = clima_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data_cover,
    aes(pca1_clima, total_cover_100, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Climate gradient") +
  geom_line(data = clima_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Soil C ----------------------------------------------------------------------
Fig_Humus_10 <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data_cover,
    aes(Corg_percent, total_cover_10, col = habitat),
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
    data = gamma_data_cover,
    aes(Corg_percent, total_cover_100, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil C") +
  geom_line(data = Humus_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Litter % --------------------------------------------------------------------
Fig_Litter_10 <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data_cover,
    aes(cover_litter, total_cover_10, col = habitat),
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
    data = gamma_data_cover,
    aes(cover_litter, total_cover_100, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Litter cover") +
  geom_line(data = Litter_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Soil pH ---------------------------------------------------------------------
Fig_pH_10 <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data_cover,
    aes(pH, total_cover_10, col = habitat),
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
    data = gamma_data_cover,
    aes(pH, total_cover_100, col = habitat),
    size = 1, alpha = 0.8, pch = 19
  ) +
  labs(y = "Plant cover", x = "Soil pH") +
  geom_line(data = pH_pred_100m, linetype = 1, linewidth = 1, col = "#D6604D")

# Grazing ---------------------------------------------------------------------
Fig_grazing_10 <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data_cover,
    aes(grazing_intencity, total_cover_10, col = habitat),
    size = 1, alpha = 0.8, pch = 19, position = position_jitter(w = 0.2)
  ) +
  labs(y = "Plant cover", x = "Grazing intencity") +
  geom_line(linetype = 5, linewidth = 0.5, col = "#50A0C8")

Fig_grazing_100 <- ggplot(grazing_pred_100m, aes(x, predicted)) +
  geom_ribbon(
    data = grazing_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D"
  ) +
  geom_point(
    data = gamma_data_cover,
    aes(grazing_intencity, total_cover_100, col = habitat),
    size = 1, alpha = 0.8, pch = 19, position = position_jitter(w = 0.2)
  ) +
  labs(y = "Plant cover", x = "Grazing intencity") +
  geom_line(data = grazing_pred_100m, linetype = 5, linewidth = 0.5, col = "#D6604D")

# Precipitation CV ------------------------------------------------------------
Fig_precipCV_10 <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(
    data = alpha_data_cover,
    aes(Prec_Varieb, total_cover_10, col = habitat),
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
    data = gamma_data_cover,
    aes(Prec_Varieb, total_cover_100, col = habitat),
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
