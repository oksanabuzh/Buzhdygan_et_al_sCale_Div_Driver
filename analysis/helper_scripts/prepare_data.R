# Purpose: produce Figure 2

# Load libraries --------------------------------------------------------------
library(tidyverse)

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

# mean per series (per 100m2 plots)
header_mean <- header %>%
  select(c(series, zonality, habitat_broad,
           where(is.numeric))) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# Prepare subset of data for alpha scale (10 m2 plots) -------------------------

alpha_data <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
            by = c("dataset", "plotID", "series", "subplot")
  )

# Remove NAs and select only needed variables
alpha_data <- alpha_data %>%
  dplyr::select(alpha_10_div, alpha_10_ENSPIE, alpha_10_cover,
                pca1_clima,
                grazing_intencity, mowing,
                cover_litter,
                Tem_range, Prec_Varieb,
                pH, Corg_percent,
                lat, lon,
                dataset, series, habitat_broad,
                subplot) %>%
  mutate(mowing = factor(mowing),
         dataset = factor(dataset),
  habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
                                                "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

# Prepare subset of data for gamma scale (100 m2 plots) -------------------------
gamma_data <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series"))

# selected variables, removed NAs
gamma_data <- gamma_data %>%
  dplyr::select(gamma_100_div, gamma_100_ENSPIE,
                gamma_100_cover,
                pca1_clima,
                grazing_intencity, mowing,
                cover_litter,
                lat, lon,
                Tem_range, Prec_Varieb, Temprt, Precipt,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>%
  mutate(
    mowing = factor(mowing),
    dataset = factor(dataset),
    habitat = fct_relevel(habitat_broad,
                               c("saline", "complex", "dry",
                                 "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

# Prepare subset of data for beta scale (100 m2 plots) -------------------------
beta_data <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series"))

# selected variables, removed NAs
beta_data <- beta_data %>%
  dplyr::select(beta_100_div, beta_100_ENSPIE,
                gamma_100_cover,
                pca1_clima,
                grazing_intencity, mowing,
                cover_litter,
                lat, lon,
                Tem_range, Prec_Varieb, Temprt, Precipt,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>%
  mutate(
    mowing = factor(mowing),
    dataset = factor(dataset),
    habitat = fct_relevel(habitat_broad,
                          c("saline", "complex", "dry",
                            "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()