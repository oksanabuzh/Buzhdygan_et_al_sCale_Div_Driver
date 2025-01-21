# Purpose: Visualize R2 and produce Figures 3b
# Fig 3 b ----


library(tidyverse)

# Read and prepare data -------------------------------------------------------

# Read R2 of all models
R2_alpha_SR <- read_csv("results/R2_alpha_SR.csv") %>% select(Effect, Rsq) %>% rename(alpha_SR = Rsq)
R2_alpha_ENSPIE <- read_csv("results/R2_alpha_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(alpha_ENSPIE = Rsq)
R2_gamma_SR <- read_csv("results/R2_gamma_SR.csv") %>% select(Effect, Rsq) %>% rename(gamma_SR = Rsq)
R2_gamma_ENSPIE <- read_csv("results/R2_gamma_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(gamma_ENSPIE = Rsq)
R2_beta_SR <- read_csv("results/R2_beta_SR.csv") %>% select(Effect, Rsq) %>% rename(beta_SR = Rsq)
R2_beta_ENSPIE <- read_csv("results/R2_beta_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(beta_ENSPIE = Rsq)

# Join all R2 values by effect
R2_part_all <- purrr::reduce(list(R2_alpha_SR, R2_alpha_ENSPIE, R2_gamma_SR, 
                          R2_gamma_ENSPIE, R2_beta_SR, R2_beta_ENSPIE), 
              full_join, by = "Effect") 

# Filter and reorder factors for plotting and add categories to the drivers
R2_part_all <- R2_part_all %>%
  filter(!Effect == "Model") %>%
  mutate(Effect = fct_relevel(Effect, c(
    "poly(pca1_clima, 2)1",
    "poly(pca1_clima, 2)2",
    "poly(Prec_Varieb, 2)1",
    "poly(Prec_Varieb, 2)2",
    "poly(pH, 2)2",
    "poly(pH, 2)1",
    "pH",
    "poly(Corg_percent, 2)1",
    "poly(Corg_percent, 2)2",
    "Corg_percent",
    "poly(cover_litter, 2)1",
    "poly(cover_litter, 2)2",
    "grazing_intencity",
    "mowing1"
  ))) %>%
  arrange(Effect) %>%
  mutate(Driver = c(
    "Climate gradient",
    "Climate gradient",
    "Precipitation variability",
    "Precipitation variability",
    "Soil pH",
    "Soil pH",
    "Soil pH",
    "Soil C",
    "Soil C",
    "Soil C",
    "Litter cover",
    "Litter cover",
    "Grazing intencity",
    "Mowing")) %>%
  mutate(Driver_all = c(
    "Productivity / stress",
    "Productivity / stress",
    "Climate seasonality",
    "Climate seasonality",
    "Productivity / stress",
    "Productivity / stress",
    "Productivity / stress",
    "Productivity / stress",
    "Productivity / stress",
    "Productivity / stress",
    "Abiotic disturbance",
    "Abiotic disturbance",
    "Anthropogenic disturbance",
    "Anthropogenic disturbance"))


# Plot variance explained by driver -------------------------------------------
# Figure 3 b ------------------------------------------------------------------
fig_3b_data <- R2_part_all %>%
  group_by(Driver) %>%
  # sum up all numeric columns by driver
  summarise(
    across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  mutate(Driver = fct_relevel(Driver, c(
    "Mowing",
    "Grazing intencity",
    "Litter cover",
    "Soil pH",
    "Soil C",
    "Precipitation variability",
    "Climate gradient"
  ))) %>%
  pivot_longer(!Driver, names_to = "Diversity", values_to = "R2") %>%
  # separate the scale and the measure
  separate_wider_delim(cols = Diversity, delim = "_", 
                       names = c("scale", "type")) %>%
  mutate(scale = fct_relevel(scale, c("alpha", "gamma", "beta")),
         type = fct_relevel(type, c("SR", "ENSPIE")))

fig_3b_data

fig_3b <- ggplot(fig_3b_data, 
  aes(
    x = scale,
    y = Driver,
    fill = scale,
    size = R2)) +
  geom_point(pch = 21, colour = "#595959") +
  facet_wrap(~type) +
  scale_size_continuous(range = c(1, 22)) +
  scale_fill_manual(values = c(
    "alpha" = "#A6CEE3",
    "gamma" = "#D6604D",
    "beta" = "#00AC7F")) +
  labs(x = "Scale", y = "Biodiversity driver", size = R^2 ~ partial) +
  theme(legend.position = "right",
    legend.key = element_rect(colour = NA, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80"),
    axis.ticks = element_blank(),
    axis.text = element_text(colour = "black", size = 15)) +
  guides(color = "none", fill = "none") 

fig_3b


# Plot by overall driver ------------------------------------------------------
# Additional plot that is not part of the paper

R2_overall_drivers <- R2_part_all %>%
  group_by(Driver_all) %>%
  # sum up all numeric columns by driver
  summarise(
    across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  mutate(Driver_all = fct_relevel(Driver_all, c("Anthropogenic disturbance",
    "Abiotic disturbance",
    # "Environmental heterogeneity",
    "Productivity / stress"
  ))) %>%
  pivot_longer(!Driver_all, names_to = "Diversity", values_to = "R2") %>%
  separate_wider_delim(cols = Diversity, delim = "_", 
                       names = c("scale", "type")) %>%
  mutate(scale = fct_relevel(scale, c("alpha", "gamma", "beta")),
         type = fct_relevel(type, c("SR", "ENSPIE")))
  
ggplot(R2_overall_drivers, aes(
  x = scale,
  y = Driver_all,
  colour = Driver_all, 
  fill = Driver_all,
  size = R2)) +
  geom_point(pch = 21) +
  scale_size_continuous(range = c(1, 22)) +
  facet_wrap(~type) +
  scale_fill_manual(values = c(
    "Anthropogenic disturbance" = "#C49A00",
    "Abiotic disturbance" = "#F8CBAD",
    "Productivity / stress" = "#AFD191",
    "Environmental heterogeneity" = "#A6CEE3")) +
  scale_colour_manual(values = c("Anthropogenic disturbance" = "#C49A00",
    "Abiotic disturbance" = "#F8CBAD",
    "Productivity / stress" = "#AFD191",
    "Environmental heterogeneity" = "#A6CEE3")) +
  labs(x = "Scale", y = "Broad drivers", size = R^2 ~ parital) +
  theme(
    legend.key = element_rect(colour = NA, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80"),
    axis.ticks = element_blank(),
    axis.text = element_text(colour = "black", size = 15)) +
  guides(color = "none", fill = "none") # remove fill in a legend
