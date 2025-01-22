# Purpose: LMM analysis for spatial aggregation as response variable

# Load libraries -----------------------------------------------------------
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
climate_PCA <- read.csv("data/climate_PCA.csv")

# Read all environmental data
header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

# Read aggregation data
aggregation <- read_csv("data/aggregation.csv")

# mean per series (per 100m2 plots)
header_mean <- header %>%
  select(c(series, lat, lon, zonality, habitat_broad,
    where(is.numeric))) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# Prepare subset of data for beta scale --------------------------------------
beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series")) %>%
  mutate(dataset = factor(dataset))

beta_data <- beta_gamma %>%
  dplyr::select(dataset, series, habitat_broad, zonality,
    gamma_100_div, gamma_100_ENSPIE, gamma_100_cover,
    beta_100_div, beta_100_ENSPIE,
    lat, lon, pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
    cover_litter,
    Tem_range, Prec_Varieb,
    pH, Corg_percent,
  ) %>%
  mutate(mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na %>%
  left_join(aggregation, by = join_by(series)) %>%
  rename(aggregation = beta.BRAY.BAL)

str(beta_data)
summary(beta_data)


# Check correlation between variables -----------------------------------------
beta_data %>%
  dplyr::select(gamma_100_div, gamma_100_ENSPIE, gamma_100_cover, ,
    beta_100_div, beta_100_ENSPIE, aggregation) %>%
  cor() |>
  ggcorrplot::ggcorrplot(
    lab = TRUE, type = "lower",
    outline.color = "white",
    colors = c("red", "white", "blue")
  )

# GLLM model of aggregation --------------------------------------------------

# Exploration ----------------------------------------------------------------

m <- lmer(aggregation ~
  poly(pca1_clima, 2) +
  poly(pH, 2) +
  poly(Corg_percent, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# check model
plot(m) # heteroscedasticity
qqnorm(resid(m))
qqline(resid(m))

Anova(m)

# Model 1: all predictors (except precipitation CV) ---------------------------
# test quadratic effects of climate, soil C, pH, and litter

# poly(pca1_clima, 2)
# poly(pH, 2) is marginal

m1_1 <- lmer(aggregation ~ poly(pca1_clima, 2) +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_2 <- lmer(aggregation ~ pca1_clima +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_3 <- lmer(aggregation ~ poly(pca1_clima, 2) +
  pH +
  Corg_percent +
  cover_litter +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

m1_4 <- lmer(aggregation ~
  pca1_clima +
  pH +
  Corg_percent +
  cover_litter +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m1_1, m1_2, m1_3, m1_4) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

# Check the best model
Anova(m1_1)
# Anova(m1_4_ENSPIE)

# Model 2: Add Precipiation variability --------------------------------------

m2_1 <- lmer(aggregation ~
  poly(pca1_clima, 2) +
  poly(Prec_Varieb, 2) +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)


m2_2 <- lmer(aggregation ~
  poly(pca1_clima, 2) +
  Prec_Varieb +
  pH +
  Corg_percent +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m2_1, m2_2) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

# Check the best model
Anova(m2_1)


# Additional analysis for precipitation CV effects  ---------------------------
# Analysis whether precipitation variability adds explanatory power beyond the
# nonlinear effect of the climate gradient:
# Model m3_1 includes both linear and quadratic terms for the climate gradient.
# Model m3_2 includes the linear effects of both climate gradient and
# precipitation variability (i.e., the quadratic term for the climate gradient 
# was replaced by precipitation variability).
# If the latter model is better (e.g. AIC is smaller than two units), we have more
# support for claim that the precipitation variability effect is shown.

m3_1 <- lmer(aggregation ~
  poly(pca1_clima, 2) +
  (1 | dataset), data = beta_data)

m3_2 <- lmer(aggregation ~
  pca1_clima +
  Prec_Varieb +
  (1 | dataset), data = beta_data)

# calculate and compare AIC
AIC(m3_1, m3_2) %>%
  arrange(+AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

Anova(m3_1)
Anova(m3_2)

# Make model plots -----------------------------------------------------------

# final models used in the plots:
Anova(m1_1)
Anova(m2_1)


# Climate plot ---------------------------------------------------------------

clima_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.betaaggr_clima <- ggplot(clima_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, 
             aes(pca1_clima, aggregation, fill = habitat, col = habitat), 
             size = 1, alpha = 0.8, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Climate gradient') +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.betaaggr_clima


# Soil C plot ----------------------------------------------------------------
Humus_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")

Fig.betaaggr_soilC <- ggplot(Humus_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, 
             aes(Corg_percent, aggregation, fill = habitat, col = habitat), 
             size = 1, alpha = 0.8, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Soil C') +
  geom_line(linetype = 5, linewidth = 0.5, col = "black")

Fig.betaaggr_soilC


# Litter % plot --------------------------------------------------------------
Litter_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")

Fig.betaaggr_Litter <- ggplot(Litter_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, 
             aes(cover_litter, aggregation, fill = habitat, col = habitat), 
             size = 1, alpha = 0.8, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Litter cover') +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.betaaggr_Litter

# Soil pH plot ---------------------------------------------------------------
pH_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms = "pH[3.7:9, by=.001]")

Fig.betaaggr_soil.pH <- ggplot(pH_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, 
             aes(pH, aggregation, fill = habitat, col = habitat), 
             size = 1, alpha = 0.8, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Soil pH') +
  geom_line(linetype = 5, linewidth = 0.5, col = "black")

Fig.betaaggr_soil.pH


# Grazing plot --------------------------------------------------------------
grazing_pred_beta_aggr <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")

Fig.betaaggr_grazing <- ggplot(grazing_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, aes(grazing_intencity, aggregation,
    fill = habitat, col = habitat),
  size = 1, alpha = 0.8, pch = 21,
  position = position_jitter(w = 0.2)) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Grazing intencity') +
  geom_line(linetype = 5, size = 0.5, col = "black")

Fig.betaaggr_grazing

# Mowing plot ---------------------------------------------------------------
Fig.betaaggr_mowing <- ggplot(beta_data, aes(mowing, aggregation)) +
  geom_boxplot(color = "black") +
  labs(y = "Species aggregation", x = 'Mowing') +
  geom_point(aes(color = habitat, fill = habitat), pch = 21, 
             position = position_jitter(w = 0.1), size = 3, alpha = 0.8) +
  scale_color_manual(values = habitat_colors) + 
  scale_fill_manual(values = habitat_colors) +
  labs(color = 'Habitat type')

Fig.betaaggr_mowing

# Precipitation variability plot -------------------------------------------
precipCV_pred_beta_aggr <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")

Fig.betaaggr_precip.CV <- ggplot(precipCV_pred_beta_aggr, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data, 
             aes(Prec_Varieb, aggregation, fill = habitat, col = habitat), 
             size = 1, alpha = 0.8, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species aggregation", x = 'Precipitation CV') +
  geom_line(linetype = 5, size = 0.5, col = "black")

Fig.betaaggr_precip.CV

# Combine all plots ----------------------------------------------------------

Fig.betaaggr_clima +
  Fig.betaaggr_precip.CV +
  Fig.betaaggr_soilC +
  Fig.betaaggr_soil.pH +
  Fig.betaaggr_Litter +
  Fig.betaaggr_grazing +
  plot_annotation(tag_levels = 'a') +
  plot_layout(ncol = 3) & 
  theme( 
    plot.tag = element_text(size = 11, face = 'bold'),
    plot.tag.position = c(0.25, 1.06))
