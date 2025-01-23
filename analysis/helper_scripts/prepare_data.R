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

alpha <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
            by = c("dataset", "plotID", "series", "subplot")
  ) %>%
  mutate(dataset = factor(dataset))

# Remove NAs and select only needed variables
alpha_data <- alpha %>%
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,
                pca1_clima,
                grazing_intencity, mowing,
                cover_litter,
                Tem_range, Prec_Varieb,
                pH, Corg_percent,
                dataset, series, habitat_broad,
                subplot) %>%
  mutate(mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
                                                "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

str(alpha_data)

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
  dplyr::select(gamma_100_div, gamma_100_ENSPIE,
                pca1_clima,
                grazing_intencity, mowing,
                cover_litter,
                Tem_range, Prec_Varieb, Temprt, Precipt,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>%
  mutate(habitat = fct_relevel(habitat_broad,
                               c("saline", "complex", "dry",
                                 "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

str(gamma_data)


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
