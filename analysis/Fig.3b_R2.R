# Purpose: summertime and visualize R2
# Fig 3 b ----

# merge tables form the (G)LMMs outputs

library(tidyverse)



# model R2: marginal and conditional

Mod_R2_alpha_SR <- read.csv("results/Mod_R2_alpha_SR.csv")
Mod_R2_alpha_ENSPIE <- read.csv("results/Mod_R2_alpha_ENSPIE.csv")
Mod_R2_gamma_SR <- read.csv("results/Mod_R2_gamma_SR.csv")
Mod_R2_gamma_ENSPIE <- read.csv("results/Mod_R2_gamma_ENSPIE.csv")
Mod_R2_beta_SR <- read.csv("results/Mod_R2_beta_SR.csv")
Mod_R2_beta_ENSPIE <- read.csv("results/Mod_R2_beta_ENSPIE.csv")

Mod_R2 <- rbind(Mod_R2_alpha_SR[1, ], Mod_R2_alpha_ENSPIE[1, ],
  Mod_R2_gamma_SR[1, ], Mod_R2_gamma_ENSPIE[1, ],
  Mod_R2_beta_SR[1, ], Mod_R2_beta_ENSPIE[1, ]
) %>%
  select(-X)

row.names(Mod_R2) = c("alpha_SR", "alpha_ENSPIE",
  "gamma_SR", "gamma_ENSPIE",
  "beta_SR", "beta_ENSPIE")
Mod_R2

write.csv(Mod_R2, "results/Mod_R2_all.csv", row.names = TRUE)



# partial R2

R2_alpha_SR <- read_csv("results/R2_alpha_SR.csv") %>% select(Effect, Rsq) %>% rename(alpha_SR = Rsq)
R2_alpha_ENSPIE <- read_csv("results/R2_alpha_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(alpha_ENSPIE = Rsq)
R2_gamma_SR <- read_csv("results/R2_gamma_SR.csv") %>% select(Effect, Rsq) %>% rename(gamma_SR = Rsq)
R2_gamma_ENSPIE <- read_csv("results/R2_gamma_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(gamma_ENSPIE = Rsq)
R2_beta_SR <- read_csv("results/R2_beta_SR.csv") %>% select(Effect, Rsq) %>% rename(beta_SR = Rsq)
R2_beta_ENSPIE <- read_csv("results/R2_beta_ENSPIE.csv") %>% select(Effect, Rsq) %>% rename(beta_ENSPIE = Rsq)


R2_part_all <- R2_alpha_SR %>%
  full_join(R2_alpha_ENSPIE, by = "Effect") %>%
  full_join(R2_gamma_SR, by = "Effect") %>%
  full_join(R2_gamma_ENSPIE, by = "Effect") %>%
  full_join(R2_beta_SR, by = "Effect") %>%
  full_join(R2_beta_ENSPIE, by = "Effect") %>%
  mutate(Effect = fct_relevel(Effect, c("poly(pca1_clima, 2)1",
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
  mutate(Driver = c("Climate gradient",
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
  mutate(Driver_all = c("Productivity / stress",
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



R2_1 <- R2_part_all %>%
  group_by(Driver) %>%
  summarise(alpha_SR = sum(alpha_SR, na.rm = T), alpha_ENSPIE = sum(alpha_ENSPIE, na.rm = T),
    gamma_SR = sum(gamma_SR, na.rm = T), gamma_ENSPIE = sum(gamma_ENSPIE, na.rm = T),
    beta_SR = sum(beta_SR, na.rm = T), beta_ENSPIE = sum(beta_ENSPIE, na.rm = T)) %>%
  mutate(Driver = fct_relevel(Driver, c("Mowing",
    "Grazing intencity",
    "Litter cover",
    "Soil pH",
    "Soil C",
    "Precipitation variability",
    "Climate gradient"
  ))) %>%
  arrange(Driver) %>%
  pivot_longer(!Driver, names_to = "Diversity", values_to = "R2") %>%
  mutate(Diversity = fct_relevel(Diversity, c("alpha_SR",
    "gamma_SR",
    "beta_SR",
    "alpha_ENSPIE", "gamma_ENSPIE",
    "beta_ENSPIE")))

R2_1

ggplot(R2_1, # %>% filter(Diversity=="alpha_SR" | Diversity=="gamma_SR" | Diversity=="beta_SR")
  aes(x = Diversity,
    y = Driver,
    fill = Diversity,
    size = R2)) +
  geom_point(pch = 21, colour = "#595959") +
  #  geom_text(aes(label = round(R2,2)),  colour = "black",  size = 5) +
  # scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1, 22)) +
  scale_fill_manual(values = c("alpha_SR" = "#A6CEE3",
    "alpha_ENSPIE" = "#A6CEE3",
    "gamma_SR" = "#D6604D",
    "gamma_ENSPIE" = "#D6604D",
    "beta_SR" = "#00AC7F",
    "beta_ENSPIE" = "#00AC7F")) +

  #  scale_color_brewer(palette =  "Paired") +
  labs(x = NULL, y = NULL, size = R^2 ~ parital) +
  theme(legend.position = "right",
    legend.key = element_rect(colour = NA, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80"),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.text.x = element_text(colour = "black", size = 15)) +
  guides(color = "none", fill = "none") # remove fill in a legend




## by Driver_all

R2_2 <- R2_part_all %>%
  group_by(Driver_all) %>%
  summarise(alpha_SR = sum(alpha_SR, na.rm = T), alpha_ENSPIE = sum(alpha_ENSPIE, na.rm = T),
    gamma_SR = sum(gamma_SR, na.rm = T), gamma_ENSPIE = sum(gamma_ENSPIE, na.rm = T),
    beta_SR = sum(beta_SR, na.rm = T), beta_ENSPIE = sum(beta_ENSPIE, na.rm = T)) %>%
  mutate(Driver_all = fct_relevel(Driver_all, c("Anthropogenic disturbance",
    "Abiotic disturbance",
    "Environmental heterogeneity",
    "Productivity / stress"
  ))) %>%
  arrange(Driver_all) %>%
  pivot_longer(!Driver_all, names_to = "Diversity", values_to = "R2") %>%
  mutate(Diversity = fct_relevel(Diversity, c("alpha_SR", "alpha_ENSPIE",
    "beta_SR", "beta_ENSPIE",
    "gamma_SR", "gamma_ENSPIE")))

R2_2




ggplot(R2_2, aes(x = Diversity,
  y = Driver_all,
  colour = Driver_all, fill = Driver_all,
  size = R2)) +
  geom_point(pch = 21) +
  # geom_text(aes(label = round(R2,2)),  colour = "black",  size = 5) +
  # scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1, 22)) +
  scale_fill_manual(values = c("Anthropogenic disturbance" = "#C49A00",
    "Abiotic disturbance" = "#F8CBAD",
    "Productivity / stress" = "#AFD191",
    "Environmental heterogeneity" = "#A6CEE3")) +
  scale_colour_manual(values = c("Anthropogenic disturbance" = "#C49A00",
    "Abiotic disturbance" = "#F8CBAD",
    "Productivity / stress" = "#AFD191",
    "Environmental heterogeneity" = "#A6CEE3")) +
  #  scale_color_brewer(palette =  "Paired") +
  labs(x = NULL, y = NULL, size = R^2 ~ parital) +
  theme(legend.position = "right",
    legend.key = element_rect(colour = NA, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80"),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.text.x = element_text(colour = "black", size = 15)) +
  guides(color = "none", fill = "none") # remove fill in a legend








# Anova results ----
glmer_alpha_SR <- read_csv("results/glmer_alpha_SR.csv")
glmer_alpha_ENSPIE <- read_csv("results/glmer_alpha_ENSPIE.csv")
glmer_gamma_SR <- read_csv("results/glmer_gamma_SR.csv")
glmer_gamma_ENSPIE <- read_csv("results/glmer_gamma_ENSPIE.csv")
glmer_beta_SR <- read_csv("results/glmer_beta_SR.csv")
glmer_beta_ENSPIE <- read_csv("results/glmer_beta_ENSPIE.csv")


# summary results -----

summary_alpha_SR <- read_csv("results/summary_alpha_SR.csv")
summary_alpha_ENSPIE <- read_csv("results/summary_alpha_ENSPIE.csv")
summary_gamma_SR <- read_csv("results/summary_gamma_SR.csv")
summary_gamma_ENSPIE <- read_csv("results/summary_gamma_ENSPIE.csv")
summary_beta_SR <- read_csv("results/summary_beta_SR.csv")
summary_beta_ENSPIE <- read_csv("results/summary_beta_ENSPIE.csv")




###############


alpha_st.eff <- read.csv("results/alpha_st.eff.csv") %>% select(-X) %>% rename(P_sign = "X.1")
gamma_st.eff <- read.csv("results/gamma_st.eff.csv") %>% select(-X) %>% rename(P_sign = "X.1")



st.eff <- rbind(alpha_st.eff, gamma_st.eff)


plot1 <- ggplot(st.eff,
  aes(x = Response,
    y = Predictor,
    fill = Std.Estimate,
    size = Std.Estimate)) +
  geom_point(pch = 21, colour = "gray") +
  #  geom_text(aes(label = round(R2,2)),  colour = "black",  size = 5) +
  # scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1, 22)) +

  #  scale_color_brewer(palette =  "Paired") +
  labs(x = NULL, y = NULL, size = R^2 ~ parital) +
  theme(legend.position = "right",
    legend.key = element_rect(colour = NA, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey80"),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.text.x = element_text(colour = "black", size = 15)) +
  guides(color = "none", fill = "none") # remove fill in a legend

plot1 +
  scale_fill_manual(values = c("Mowing" = "#F8CBAD",
    "Grazing intencity" = "#F8CBAD",
    "Litter cover" = "#F8CBAD",
    "Soil C" = "#AFD191",
    "Soil pH" = "#AFD191",
    "Temperature range" = "#A6CEE3",
    "Precipitation variation" = "#A6CEE3",
    "Climate" = "#AFD191"))
