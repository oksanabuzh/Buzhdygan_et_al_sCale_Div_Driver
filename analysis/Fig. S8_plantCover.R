# Purpose: Figure S8

# libraries----
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(performance)
library(patchwork)


# data----
## data alpha----
alpha <- read_csv("Data/alpha_GLM.csv")
str(alpha)
names(alpha)
# dataset is a separate vegetation survey campaign
alpha$dataset <- factor(alpha$dataset)

# Remove NAs
alpha_data <- alpha %>%
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
    dataset, series, habitat_broad,
    subplot) %>%
  mutate(Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%

  drop_na

str(alpha_data)

# plot on a mean alpha per series to omit pseudoreplication of the plots:
alpha_mean <- read_csv("results/Div_NMDS_BRAY_Jaccard_Dataset.csv") %>%
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total, inclination,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
    dataset, series, habitat_broad) %>%
  mutate(Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na

str(alpha_mean)

## data gamma----
beta_gamma <- read_csv("Data/beta_gamma_GLM.csv")
str(beta_gamma)
names(beta_gamma)

# dataset is a separate vegetation survey campaign
beta_gamma$dataset <- factor(beta_gamma$dataset)

# Remove NAs
gamma_data <- beta_gamma %>%
  dplyr::select(gamma_100_div, gamma_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
    cover_litter,
    BIO7, BIO15,
    pH, Corg_percent,
    dataset, series, habitat_broad, zonality) %>%
  mutate(Tem_range = BIO7,
    Prec_Varieb = BIO15,
    mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na


# Plant cover -----

# alpha scale
tot_cover_10 <- read_csv("data/species_matrix_total.csv") %>%
  filter(scale == 10) %>%
  rowwise() %>%
  mutate(total_cover_10 = sum(c_across("Abietinella abietina":"Tribulus terrestris"), na.rm = T)) %>%
  ungroup() %>%
  select(dataset, series, subplot, total_cover_10) %>%
  mutate(dataset = factor(dataset))

str(alpha_data)
alpha_data_cover <- alpha_data %>%
  left_join(tot_cover_10, by = c("dataset", "series", "subplot"))

alpha_data_cover$total_cover_10
str(alpha_data_cover)

# gamma scale

tot_cover_100 <- read_csv("data/species_matrix_total.csv") %>%
  filter(scale == 100) %>%
  rowwise() %>%
  mutate(total_cover_100 = sum(c_across("Abietinella abietina":"Tribulus terrestris"), na.rm = T)) %>%
  ungroup() %>%
  select(dataset, series, total_cover_100) %>%
  mutate(dataset = factor(dataset))

gamma_data_cover <- gamma_data %>%
  left_join(tot_cover_100, by = c("dataset", "series"))

gamma_data_cover$total_cover_100


# GLMM ----

## alpha----


m1_1 <- lmer(sqrt(total_cover_10) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

# check model
plot(m1_1) # heteroscedasticity
qqnorm(resid(m1_1))
qqline(resid(m1_1))

summary(m1_1)
Anova(m1_1)

m1_2 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

summary(m1_2)
Anova(m1_2)

plot_model(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)
plot_model(m1_2, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


# calculate AIC
tlist <- AIC(m1_1, m1_2)
# tlist
arrange(tlist, +AIC)


Anova(m1_1)
# Anova(m1_1)

## Add Prec_Varieb ----
### Select model for the Prec_Varieb ----
m2_1 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)



m2_2 <- lmer(sqrt(total_cover_10) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset / series), data = alpha_data_cover)

# calculate AIC
tlist2 <- AIC(m2_1, m2_2)
# tlist
arrange(tlist2, +AIC)


Anova(m2_1)

### Predictions alpha----
clima_pred_10m <- get_model_data(m1_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_10m <- get_model_data(m1_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_10m <- get_model_data(m1_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_10m <- get_model_data(m1_1, type = "pred", terms = "pH[3.8:9, by=.001]")
grazing_pred_10m <- get_model_data(m1_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_10m <- get_model_data(m2_1, type = "pred", terms = "Prec_Varieb")


## gamma----


m3_1 <- lmer(sqrt(total_cover_100) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

# check model
plot(m3_1) # heteroscedasticity
qqnorm(resid(m3_1))
qqline(resid(m3_1))

summary(m3_1)
Anova(m3_1)

m3_2 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

summary(m3_2)
Anova(m3_2)



m3_3 <- lmer(sqrt(total_cover_100) ~
  poly(pca1_clima, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  cover_litter +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

Anova(m3_3)

# calculate AIC
tlist <- AIC(m3_1, m3_2, m3_3)
# tlist
arrange(tlist, +AIC)


Anova(m3_1)
plot_model(m3_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]", show.data = T)


## Add Prec_Varieb ----
### Select model for the Prec_Varieb ----
m4_1 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  poly(Prec_Varieb, 2) +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)



m4_2 <- lmer(sqrt(total_cover_100) ~
  pca1_clima +
  Prec_Varieb +
  poly(Corg_percent, 2) +
  poly(pH, 2) +
  poly(cover_litter, 2) +
  grazing_intencity + mowing +
  (1 | dataset), data = gamma_data_cover)

# calculate AIC
tlist2 <- AIC(m4_1, m4_2)
# tlist
arrange(tlist2, +AIC)


Anova(m4_1)


## Predictions gamma ----

clima_pred_100m <- get_model_data(m3_1, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Humus_pred_100m <- get_model_data(m3_1, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Litter_pred_100m <- get_model_data(m3_1, type = "pred", terms = "cover_litter[0:100, by=0.01]")
pH_pred_100m <- get_model_data(m3_1, type = "pred", terms = "pH[3.7:9, by=.001]")
grazing_pred_100m <- get_model_data(m3_1, type = "pred", terms = "grazing_intencity[-0.2:3.2, by=0.01]")
precipCV_pred_100m <- get_model_data(m4_1, type = "pred", terms = "Prec_Varieb")





# 2) Plot the predictions  ----

#         saline    complex       dry       wet       mesic        fringe       alpine
col = c("#4e3910", "#CC6600", "#e3c28b", "#CC99FF", "#0066FF", "#00B200", "#006600")


# alpha-gamma----

## SR----

### Clima ----

Fig_Clima_10 <- ggplot(clima_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(pca1_clima, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Climate gradient') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8")



Fig_Clima_100 <- ggplot(clima_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = clima_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(pca1_clima, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Climate gradient') +
  geom_line(data = clima_pred_100m, linetype = 5, size = 0.5, col = "#D6604D")

### Soil C ----
Fig_Humus_10 <- ggplot(Humus_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(Corg_percent, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Soil C') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8")


Fig_Humus_100 <- ggplot(Humus_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = Humus_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(Corg_percent, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Soil C') +
  geom_line(data = Humus_pred_100m, linetype = 1, size = 1, col = "#D6604D")


### Litter % ----

Fig_Litter_10 <- ggplot(Litter_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(cover_litter, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Litter cover') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8")

Fig_Litter_100 <- ggplot(Litter_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = Litter_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(cover_litter, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Litter cover') +
  geom_line(data = Litter_pred_100m, linetype = 5, size = 0.5, col = "#D6604D")

### Soil pH ----

Fig_pH_10 <- ggplot(pH_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(pH, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Soil pH') +
  geom_line(linetype = 1, size = 1, col = "#50A0C8")

Fig_pH_100 <- ggplot(pH_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = pH_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(pH, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Soil pH') +
  geom_line(data = pH_pred_100m, linetype = 1, size = 1, col = "#D6604D")


### Grazing ----

Fig_grazing_10 <- ggplot(grazing_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(grazing_intencity, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21, position = position_jitter(w = 0.2)) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Grazing intencity') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8")


Fig_grazing_100 <- ggplot(grazing_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = grazing_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(grazing_intencity, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21, position = position_jitter(w = 0.2)) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Grazing intencity') +
  geom_line(data = grazing_pred_100m, linetype = 5, size = 0.5, col = "#D6604D")


### Prec_Varieb ----

Fig_precipCV_10 <- ggplot(precipCV_pred_10m, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#A6CEE3") +
  geom_point(data = alpha_data_cover,
    aes(Prec_Varieb, total_cover_10, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Precipitation CV') +
  geom_line(linetype = 5, size = 0.5, col = "#50A0C8")

Fig_precipCV_100 <- ggplot(precipCV_pred_100m, aes(x, predicted)) +
  geom_ribbon(data = precipCV_pred_100m, aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.1, fill = "#D6604D") +
  geom_point(data = gamma_data_cover,
    aes(Prec_Varieb, total_cover_100, fill = habitat, col = habitat), size = 1, alpha = 0.8,
    pch = 21) +
  scale_fill_manual(values = col) + scale_color_manual(values = col) +
  labs(y = "Plant cover", x = 'Precipitation CV') +
  geom_line(data = precipCV_pred_100m, linetype = 5, size = 0.5, col = "#D6604D")



## Merge plots----

library(patchwork)

set_theme(base = theme_bw(), axis.textsize.x = 0.6, axis.textsize.y = 0.6, axis.textcolor = "black",
  axis.title.color = "black", axis.title.size = 0.65, legend.pos = "None",
  geom.linetype = 2)
FigS8 <-
  Fig_Clima_10 + Fig_Clima_100 +
    Fig_precipCV_10 + Fig_precipCV_100 +
    Fig_Humus_10 + Fig_Humus_100 +
    Fig_pH_10 + Fig_pH_100 +
    Fig_Litter_10 + Fig_Litter_100 +
    Fig_grazing_10 + Fig_grazing_100 +
    #
    #
    plot_annotation(tag_levels = 'a') +
    plot_layout(ncol = 2) & ylab(NULL) & theme(plot.margin = margin(3, 1, 3, 20),
    plot.tag = element_text(size = 6, face = 'bold'),
    plot.tag.position = c(0.15, 1.06))


FigS8
