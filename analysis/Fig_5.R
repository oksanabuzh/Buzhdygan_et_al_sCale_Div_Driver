# Purpose: Produce Fig. 5 - relationships of beta SR and beta ENSPIE
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
# This script prepared the following necessary data for all analyses
# - alpha_data: diversity, ENSPIE, cover and environmental variables for the 10m2 plots
# - gamma_data: diversity, ENSPIE, cover and environmental variables for the
#    100m2 plots
# - beta_data: diversity, ENSPIE, cover and environmental variables for the
#    beta scale
source("analysis/helper_scripts/prepare_data.R")

# Turn mowing variable into a factor
beta_data <- mutate(beta_data, mowing = factor(mowing))

# Check how the dataset looks like
beta_data

# Read aggregation data
aggregation <- read_csv("data/aggregation.csv")

# Join beta data with aggregation data
beta_data <- left_join(beta_data, aggregation, by = "series") %>% 
  rename(aggregation = beta.BRAY.BAL)

# add gamma enspie and diversity values to beta data
beta_data <- left_join(beta_data, gamma_data %>%  select(series, gamma_100_div, gamma_100_ENSPIE), by = "series")

# Correlation among measures --------------------------------------------------
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
  geom_line(linewidth = 0.5, linetype = "longdash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Total cover')


# aggregation
Fig_ENSPIE_aggr <- ggplot(
  ggeffects::ggpredict(beta_ENSPIE_1, terms = c("aggregation[0:1, by=.001]")),
  aes(x, predicted)) +
  geom_point(data = beta_data, aes(aggregation, beta_100_ENSPIE),
    size = 1.5, alpha = 0.8, color = "#00AC7F") +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  labs(y = "Beta evenness", x = 'Intraspecific aggregation')

Fig_ENSPIE_aggr


# Combine all plots ----------------------------------------------------------

Fig_SR_cover +
  Fig_SR_evenness +
  Fig_SR_aggr + Fig_ENSPIE_aggr +
  plot_annotation(tag_levels = 'a') +
  plot_layout(ncol = 2)  & 
  theme( 
  plot.tag = element_text(size = 11, face = 'bold'),
  plot.tag.position = c(0.25, 1.06))
