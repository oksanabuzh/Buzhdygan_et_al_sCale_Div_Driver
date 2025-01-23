# Purpose: Produce FIg. 5 - relationships of beta SR and beta ENSPIE
#                           with the spatial aggregation proxy


# Load libraries -------------------------------------------------------------
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(sjPlot)
library(performance)
library(patchwork)

# Set theme for the plots
set_theme(base = theme_bw(), 
          axis.textsize.x = 0.8, 
          axis.textsize.y = 0.8, 
          axis.textcolor = "black",
          axis.title.color = "black",
          axis.title.size = 0.9, 
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
climate_PCA <- read_csv("data/climate_PCA.csv")

# Read all environmental data
header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read_csv("data/climate_PCA.csv"),
    by = "series"
  )

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

# Correlation among measures------
beta_data %>%
  select( # gamma_100_div, gamma_100_ENSPIE,
    beta_100_div, beta_100_ENSPIE, gamma_100_cover,
    aggregation) %>%
  rename("beta SR" = beta_100_div,
    "beta evenness" = beta_100_ENSPIE,
    "total plant cover" = gamma_100_cover,
    "species aggregation" = aggregation) %>%
  cor() |>
  ggcorrplot::ggcorrplot(
    lab = TRUE, type = "lower",
    outline.color = "white",
    colors = c("red", "white", "blue")
  )

# (1) beta SR model ----------------------------------------------------------

beta1 <- lmer(beta_100_div ~
  poly(gamma_100_cover, 2) +
  gamma_100_ENSPIE +
  aggregation +
  (1 | dataset), data = beta_data)


Anova(beta1)
summary(beta1)
check_collinearity(beta1)

# Plots -----------------------------------------------------------------------

# gamma_100_cover
Fig_SR_cover <- ggplot(ggeffects::ggpredict(beta1, terms = c("gamma_100_cover")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_cover, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness", x = 'Total cover')

Fig_SR_cover

# aggregation
Fig_SR_aggr <- ggplot(ggeffects::ggpredict(beta1, terms = c("aggregation[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(aggregation, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness", x = 'Intraspecific aggregation')

Fig_SR_aggr

# evenness
Fig_SR_evenness <- ggplot(ggeffects::ggpredict(beta1, terms = c("gamma_100_ENSPIE[0:25, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_ENSPIE, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(linewidth = 0.5, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness",
    x = expression(paste("Species evenness (at 100", m^{2}, ")")))

Fig_SR_evenness

# (2) beta ENSPIE model ------------------------------------------------------
beta_ENSPIE_1 <- lmer(beta_100_ENSPIE ~
  gamma_100_cover +
  aggregation +
  (1 | dataset), data = beta_data)

Anova(beta_ENSPIE_1)
summary(beta_ENSPIE_1)
check_collinearity(beta_ENSPIE_1)


# Plots -----------------------------------------------------------------------

# gamma_100_cover
ggplot(ggeffects::ggpredict(beta_ENSPIE_1, terms = c("gamma_100_cover")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_cover, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(size = 0.5, linetype = "longdash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Total cover')


# aggregation
Fig_ENSPIE_aggr <- ggplot(
  ggeffects::ggpredict(beta_ENSPIE_1, terms = c("aggregation[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(aggregation, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(size = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Intraspecific aggregation')

Fig_ENSPIE_aggr


# Combine all plots ----------------------------------------------------------

Fig_SR_cover +
  Fig_SR_evenness +
  Fig_SR_aggr + Fig_ENSPIE_aggr +
  plot_annotation(tag_levels = 'a') +
  plot_layout(ncol = 2) # & #ylab(NULL) &
theme( # plot.margin = margin(3, 1, 3, 20),
  plot.tag = element_text(size = 11, face = 'bold'),
  plot.tag.position = c(0.25, 1.06))
