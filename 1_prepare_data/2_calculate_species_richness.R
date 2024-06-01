# Puropse: Calculate species richness per plot, scale and species group from the
# cleaned raw data
library(tidyverse)
library(patchwork)

sp_data <- read_csv("data/database_analysis.csv")

# Optional: Turn off scientific notation for numbers
# options(scipen = "999")

grass_div <- sp_data %>%
  select(-c(layer, value, response)) %>%
  distinct() %>%
  mutate(scale = as.factor(as.character(scale))) %>%
  group_by(dataset, series, subplot, scale, group) %>%
  summarise(S = n_distinct(species)) %>%
  pivot_wider(names_from = group, values_from = S, values_fill = 0) %>%
  mutate(
    total_rich = sum(B + L + VP),
    crypto_rich = sum(B + L)
  ) %>%
  mutate(
    scale = as.numeric(as.character(scale))
  )

# plot species-area curves to check if they make sense
grass_div %>%
  group_by(dataset, series, scale) %>%
  summarize(
    B = mean(B),
    VP = mean(VP),
    L = mean(L),
    total_rich = mean(total_rich),
    crypto_rich = mean(crypto_rich)
  ) %>%
  ggplot(aes(x = scale, y = total_rich, group = series)) +
  geom_point() +
  geom_line() +
  facet_wrap(~dataset) +
  scale_x_log10()

# richness at different scales
grass_div %>%
  ggplot(aes(y = total_rich, x = factor(scale))) +
  geom_boxplot()

# write species richness table
write_csv(grass_div, "data/grass_rich.csv")
