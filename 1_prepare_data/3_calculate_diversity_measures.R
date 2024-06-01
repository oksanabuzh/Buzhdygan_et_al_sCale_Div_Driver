# Purpose: Calculate alpha, beta, gamma diversity measures for the plots for
# analysis
library(tidyverse)
library(patchwork)
sp_data <- read_csv("data/database_analysis.csv")

# Sum up cover of the same species in different layers (shrub + herb layer)
cover_prep <- sp_data %>%
  filter(response == "cover") %>%
  group_by(dataset, group, species, series, subplot, scale) %>%
  summarize(value = sum(value))

# Calculate richness and ENSPIE 100 m2 ----------------------------------------

gamma_100 <- cover_prep %>%
  filter(scale == 100) %>%
  group_by(dataset, series, scale) %>%
  summarise(
    gamma_div = n_distinct(species),
    gamma_ENSPIE = vegan::diversity(value, index = "invsimpson")
  ) %>%
  ungroup()

# Calculate richness and ENSPIE 10 m2 ----------------------------------------

alpha <- cover_prep %>%
  filter(scale == 10) %>%
  group_by(dataset, series, subplot, scale) %>%
  summarise(
    alpha_div = n_distinct(species),
    alpha_ENSPIE = vegan::diversity(value, index = "invsimpson")
  ) %>%
  ungroup()

# Calculate richness and ENSPIE on 20 m2 -------------------------------------

gamma_20 <- cover_prep %>%
  filter(scale == 10) %>%
  mutate(scale = 20) %>%
  group_by(dataset, series, species, scale) %>%
  summarize(value = sum(value)) %>% # sum up the cover of the 2 10 m2 plots
  group_by(dataset, series, scale) %>%
  summarise(
    gamma_div = n_distinct(species),
    gamma_ENSPIE = vegan::diversity(value, index = "invsimpson")
  ) %>%
  ungroup()

# Calculate beta diversity for  20 m2 ---------------------------------------
beta_20 <- gamma_20 %>%
  select(-scale) %>%
  left_join(alpha, multiple = "all") %>%
  group_by(dataset, series, gamma_div, gamma_ENSPIE) %>% 
  summarize(alpha_div = mean(alpha_div),
         alpha_ENSPIE = mean(alpha_ENSPIE)) %>% 
  ungroup() %>% 
  mutate(
    beta_div = (gamma_div / alpha_div), # subplot beta (20m2)
    beta_ENSPIE = (gamma_ENSPIE / alpha_ENSPIE), # sub plot beta (20m2)
  ) %>%
  select(-c(gamma_div, alpha_div, gamma_ENSPIE, alpha_ENSPIE)) %>%
  mutate(scale = 20)

# Calculate beta diversity for  100 m2 ---------------------------------------

beta_100 <- gamma_100 %>%
  select(-scale) %>%
  left_join(alpha, multiple = "all") %>%
  group_by(dataset, series, gamma_div, gamma_ENSPIE) %>% 
  summarize(alpha_div = mean(alpha_div),
            alpha_ENSPIE = mean(alpha_ENSPIE)) %>% 
  ungroup() %>% 
  mutate(
    beta_div = (gamma_div / alpha_div), # plot beta (100m2)
    beta_ENSPIE = (gamma_ENSPIE / alpha_ENSPIE), # plot beta (100m2)
  ) %>%
  select(-c(gamma_div, alpha_div, gamma_ENSPIE, alpha_ENSPIE)) %>%
  mutate(scale = 100)

# combine beta diversity
beta_div <- beta_20 %>%
  bind_rows(beta_100) %>%
  pivot_longer(c(beta_div, beta_ENSPIE),
    names_to = "metric",
    values_to = "value"
  )

alpha <- alpha %>%
  pivot_longer(c(alpha_div, alpha_ENSPIE),
    names_to = "metric",
    values_to = "value"
  )

gamma_20 <- gamma_20 %>%
  pivot_longer(c(gamma_div, gamma_ENSPIE),
               names_to = "metric",
               values_to = "value"
  )

gamma_100 <- gamma_100 %>% 
  pivot_longer(c(gamma_div, gamma_ENSPIE),
               names_to = "metric",
               values_to = "value"
  )

div <- alpha %>%
  bind_rows(gamma_20) %>%
  bind_rows(gamma_100) %>%
  bind_rows(beta_div) %>% 
  separate(metric, into = c("type", "metric")) %>% 
  mutate(scale = factor(scale, levels = c(10, 20, 100))) %>%
  mutate(metric = fct_relevel(metric, "div", "ENSPIE"))

# Make figures ------------------------------------------------------------

rich_fig <- ggplot(
  data = div %>% filter(metric == "div" & type != "beta"),
  aes(x = scale, y = value)
) +
  geom_point(position = position_jitter(0.2), alpha = 0.5, color = "#C0C0C0") +
  geom_boxplot(aes(color = scale), alpha = 0.2) +
  labs(subtitle = "Species richness") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
# Careful, some high outlier values for 100 m2 scale are removed by ylim
pie_fig <- ggplot(
  data = div %>% filter(metric == "ENSPIE" & type != "beta"),
  aes(x = scale, y = value)
) +
  geom_point(position = position_jitter(0.2), alpha = 0.5, color = "#C0C0C0") +
  geom_boxplot(aes(color = scale), alpha = 0.2) +
  ylim(0, 30) +
  labs(subtitle = "ENSPIE") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Careful, some outliers are removed by ylim
beta_fig_20 <- ggplot(
  data = div %>% filter(type == "beta" & scale == 20),
  aes(x = factor(metric), y = value)
) +
  geom_point(position = position_jitter(0.2), alpha = 0.5, color = "#C0C0C0") +
  geom_boxplot(alpha = 0.2, color = "darkblue") +
  ylim(0, 2) +
  labs(subtitle = (expression(paste(italic(beta), "-Diversity - 20m2",
    sep = ""
  )))) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

# Careful, some data points are removed for high values by ylim
beta_fig_100 <- ggplot(
  data = div %>% filter(type == "beta" & scale == 100),
  aes(x = factor(metric), y = value)
) +
  geom_point(position = position_jitter(0.2), alpha = 0.5, color = "#C0C0C0") +
  geom_boxplot(color = "darkblue", alpha = 0.2) +
  ylim(0, 4) +
  labs(subtitle = (expression(paste(italic(beta), "-Diversity - 100m2",
    sep = ""
  )))) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

(rich_fig + pie_fig) / (beta_fig_20 + beta_fig_100)


# Write diversity measures ------------------------------------------------

write_csv(div, "data/alpha_beta_gamma.csv")
