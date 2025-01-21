# Purpose: Make Fig 3 a, Fig. 4
# Barplots of the standardized effects of diversity drivers
# on alpha and gamma richnes and ENSPIE Fig 3 a
# Scatterplot of alpha vs gamma Fig. S4

library(tidyverse)
library(MetBrewer)
library(patchwork)

# Read and prepare data -------------------------------------------------------
alpha_st.eff <- read_csv("results/alpha_st.eff.csv") %>% 
  select(Response, Predictor, Std.Estimate)
gamma_st.eff <- read_csv("results/gamma_st.eff.csv") %>% 
  select(Response, Predictor, Std.Estimate)

alpha <- alpha_st.eff %>%
  rename(st.est = "Std.Estimate") %>%
  mutate(Measure = dplyr::recode(Response,
    'alpha_10_div' = 'SR',
    'log(alpha_10_ENSPIE)' = 'ENSPIE'))

gamma <- gamma_st.eff %>%
  rename(st.est = "Std.Estimate") %>%
  mutate(Measure = dplyr::recode(Response,
    'gamma_100_div' = 'SR',
    'log(gamma_100_ENSPIE)' = 'ENSPIE'))

# Bind datasets
Difference <- bind_rows(alpha = alpha, gamma = gamma, .id = "scale") %>%
  mutate(Driver = dplyr::recode(Predictor,
    clima_Comp = 'Climate gradient',
    PrecipCV_Comp = 'Precipitation variation',
    Corg_Comp = 'Soil C',
    Corg_percent = 'Soil C',
    pH = 'Soil pH',
    pH_Comp = 'Soil pH',
    Litter_Comp = 'Litter cover',
    grazing_intencity = 'Grazing intencity',
    mowing = 'Mowing')) %>%
  mutate(st.est_abs = abs(st.est))

Difference

# Turn character columns into factors in the right order for the plots 
Difference$Driver <- factor(Difference$Driver,
  levels = rev(c("Climate gradient",
    "Precipitation variation",
    "Soil C",
    "Soil pH",
    "Litter cover",
    "Grazing intencity",
    "Mowing"))
)

Difference$Measure <- factor(Difference$Measure, levels = c("SR", "ENSPIE"))
Difference$scale <- factor(Difference$scale, levels = c("gamma", "alpha"))

# Fig. 3 a --------------------------------------------------------------------

fig_3a <- ggplot(Difference, aes(y = Driver, x = st.est_abs, fill = scale)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  facet_wrap(~Measure, scales = "free_x", 
             labeller = labeller(Measure = c(
               SR = "Species richness", 
               ENSPIE = "ENSPIE"))) +
  scale_fill_manual(values = c("alpha" = "#A6CEE3",
    "gamma" = "#D6604D")) +
  labs(y = "Biodiversity driver", x = "Standardized effect size") +
  # recode the facet labels to be SR = Species richness and ENSPIE to be ENS and then PIE as subscript
    theme_bw() +
  theme(
    axis.title = element_text(face = "bold")
  )
fig_3a


# Figure 4 -------------------------------------------------------------------

# recalculate difference to get original oder of variables
Difference <- bind_rows(alpha = alpha, gamma = gamma, .id = "scale") %>%
  mutate(Driver = dplyr::recode(Predictor,
                                clima_Comp = 'Climate gradient',
                                PrecipCV_Comp = 'Precipitation variation',
                                Corg_Comp = 'Soil C',
                                Corg_percent = 'Soil C',
                                pH = 'Soil pH',
                                pH_Comp = 'Soil pH',
                                Litter_Comp = 'Litter cover',
                                grazing_intencity = 'Grazing intencity',
                                mowing = 'Mowing')) %>%
  mutate(st.est_abs = abs(st.est))

fig_4 <- Difference %>% 
  select(scale, Measure, Driver, st.est_abs) %>%
  pivot_wider(names_from = scale, values_from = st.est_abs) %>% 
  mutate(
    Measure = recode(
      Measure,
      SR = "Species richness",
      ENSPIE = "Evenness"
    )
  ) %>% 
  ggplot(aes(x = alpha, y = gamma,
  color = Driver, shape = Measure)) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 1) +
  geom_point(size = 4, stroke = 1.1) +
  scale_shape_manual(values = c(19, 21)) +
  theme_bw() +
  labs(
    x = expression(paste("Effect size at 10-", m^2, "plots")),
    y = expression(paste("Effect size at 100-", m^2, "plots")),
    shape = "Diversity measure",
    color = "Driver")
fig_4

