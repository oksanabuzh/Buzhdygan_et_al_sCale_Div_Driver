# Purpose: Relationships of beta SR and beta ENSPIE with the spatial aggregation proxy


# libraries----
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(performance)
library(patchwork)


# prepare data----

# "data/climate_PCA.csv" contains scores for the compound climate variable,
# derived from the PCA analysis in "1_prepare_data/ PCA_environment.R"

# "data/Environm_variabl.csv" contains all environmental data

bray_curtis <- read_csv("data/aggregation.csv")
# for information on bray indices see
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12693

climate_PCA <- read.csv("data/climate_PCA.csv")

header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header)
names(header)


header_mean <- header %>%
  select(c(series, lat, lon, zonality, habitat_broad,
    where(is.numeric))) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()


beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series")) %>%
  mutate(dataset = factor(dataset))

str(beta_gamma)
names(beta_gamma)

# dataset is a separate vegetation survey campaign
beta_gamma$dataset

## data beta----
beta_data <- beta_gamma %>%
  dplyr::select(dataset, series, habitat_broad, zonality,
    gamma_100_div, gamma_100_ENSPIE, gamma_100_cover,
    beta_100_div, beta_100_ENSPIE,
    lat, lon, pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
  ) %>%
  mutate(Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na %>%
  left_join(bray_curtis, by = join_by(series)) %>%
  mutate(BRAY_tn = beta.BRAY.BAL, # balanced variation in abundance,
    # whereby the individuals of some species in one site
    # are substituted by the same number of individuals of different species in another site
    BRAY_nest = beta.BRAY.GRA, # some individuals are lost from one site to the other
    BRAY_dissiml = beta.BRAY)

str(beta_data)
summary(beta_data)


# Correlation among measures------
beta_data %>%
  select( # gamma_100_div, gamma_100_ENSPIE,
    beta_100_div, beta_100_ENSPIE, gamma_100_cover,
    #  mean_corner_diff, # mean value of between-subplot differences in cover among species
    #  sd_corner_diff,   # sd value of between-subplot differences in cover among species
    #  aggreg,
    BRAY_tn, # BRAY_nest, # BRAY_dissiml,
    # Jaccard_tn, Jaccard_nest #, Jaccard_dissiml
  ) %>%
  rename("beta SR" = beta_100_div,
    "beta evenness" = beta_100_ENSPIE,
    "total plant cover" = gamma_100_cover,
    #  "species turnover" = Jaccard_tn,
    #   "species nestedness" = Jaccard_nest,
    #  "Total occurance-based dissimilarity" = Jaccard_dissiml,
    "species aggregation" = BRAY_tn) %>%
  # rename(SR=sowndiv, "FG richness"=numfg) %>%
  cor() |>
  ggcorrplot::ggcorrplot(
    lab = TRUE, type = "lower",
    outline.color = "white",
    colors = c("red", "white", "blue")
  )

# (1) beta SR -----


beta1 <- lmer(beta_100_div ~
  poly(gamma_100_cover, 2) +
  gamma_100_ENSPIE +
  BRAY_tn + # BRAY_nest +
  # BRAY_dissiml +
  (1 | dataset), data = beta_data)


Anova(beta1)
summary(beta1)
check_collinearity(beta1)



## Plots ----

# gamma_100_cover
Fig_SR_cover <- ggplot(ggeffects::ggpredict(beta1, terms = c("gamma_100_cover")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_cover, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(size = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness", x = 'Total cover')

Fig_SR_cover

# BRAY_tn
Fig_SR_aggr <- ggplot(ggeffects::ggpredict(beta1, terms = c("BRAY_tn[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(BRAY_tn, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(size = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness", x = 'Intraspecific aggregation')

Fig_SR_aggr

# evenness
Fig_SR_evenness <- ggplot(ggeffects::ggpredict(beta1, terms = c("gamma_100_ENSPIE[0:25, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_ENSPIE, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(size = 0.5, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness",
    x = expression(paste("Species evenness (at 100", m^{2}, ")")))

Fig_SR_evenness


# BRAY_dissiml
ggplot(ggeffects::ggpredict(beta1, terms = c("BRAY_dissiml[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(BRAY_dissiml, beta_100_div),
    size = 1.5, alpha = 0.8, color = "#00AC7F", shape = 21, stroke = 0.8) +
  geom_line(size = 0.5, linetype = "longdash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta species richness", x = 'Cover-based dissimilarity')





# (2) beta ENSPIE -----
beta_ENSPIE_1 <- lmer(beta_100_ENSPIE ~
  gamma_100_cover +
  BRAY_tn + # BRAY_nest +
  # BRAY_dissiml +
  (1 | dataset), data = beta_data)

Anova(beta_ENSPIE_1)
summary(beta_ENSPIE_1)
check_collinearity(beta_ENSPIE_1)


# Plots ----

# gamma_100_cover
ggplot(ggeffects::ggpredict(beta_ENSPIE_1, terms = c("gamma_100_cover")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(gamma_100_cover, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(size = 0.5, linetype = "longdash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Total cover')


# BRAY_dissiml
ggplot(ggeffects::ggpredict(beta_ENSPIE_1, terms = c("BRAY_dissiml[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(BRAY_dissiml, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(size = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Cover-based dissimilarity')


# BRAY_tn
ggplot(ggeffects::ggpredict(beta_ENSPIE_1, terms = c("BRAY_tn[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(BRAY_tn, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(size = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Intraspecific aggregation')
