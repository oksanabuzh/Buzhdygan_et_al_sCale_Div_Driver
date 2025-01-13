# Purpose: Descriptive and summary statistics. 
# In this script we will prepare data for the following plots:
# - Figure 1 b: Species richness for each scale and habitat type
# - Fig. S7: Species richness at 100-m² plots for each habitat type
# - Fig. S7b: Beta diversity at 100-m² plots for each habitat type
# - Fig. S3 a: Species richness for each scale
# - Fig. S3 b: ENSPIE for each scale
# - Fig. S3 c: R2 models for SR
# - Fig. S3 d: R2 models for ENSPIE
# - Fig. S1 a: Correlation among Temperature and precipitation
# - Fig. S1 b: Correlation among precipitation variability and climate PC
# - Fig. S2 a: Correlation among soil pH and climate PC
# - Fig. S2 b: Correlation among soil C and climate PC
# - Fig. S2 c: Correlation among litter cover and climate PC
# - Fig. S2 d: Correlation among soil C and litter cover
# - Fig. 43 e: Plant cover

library(tidyverse)
library(sjPlot)
library(lme4)

# Define habitat colors to be used in all plots to distinguish habitats
# define habitat colors
habitat_colors = c(
  saline = "#4e3910",
  complex = "#CC6600",
  dry = "#e3c28b",
  wet = "#CC99FF",
  mesic = "#0066FF",
  fringe = "#00B200",
  alpine = "#006600")

# Read and prepare data -------------------------------------------------------

# Read climate PCA data
climate_PCA <- read_csv("data/climate_PCA.csv")

# Read header data with environmental variables for each series and subplot
header <- read_csv("data/Environm_variabl.csv") %>%
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header)
names(header)

# Calculate mean values for 2 corners for each series
header_mean <- header %>%
  select(c(series, zonality, habitat_broad,
    where(is.numeric))) %>%
  group_by(series, zonality, habitat_broad) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# Species richness and ENSPIE data for alpha and gamma scale -------------------

# Prepare data for gamma scale ------------------------------------------------

# Beta and gamma diversity  
beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series"))

str(beta_gamma)

# dataset is a separate vegetation survey campaign
beta_gamma$dataset <- factor(beta_gamma$dataset)

# Remove NAs and select only needed columns for gamma scale (100 m²)
gamma_data <- beta_gamma %>%
  dplyr::select(gamma_100_div, gamma_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    BIO7, BIO15, BIO1, BIO12,
    pH, Corg_percent,
    dataset, series, habitat_broad, zonality) %>%
  mutate(Tem_range = BIO7,
    Prec_Varieb = BIO15,
    Temprt = BIO1,
    Precipt = BIO12,
    mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad, c("saline", "complex", "dry",
    "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

str(gamma_data)

# Prepare data for alpha scale ------------------------------------------------
alpha <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
    by = c("dataset", "plotID", "series", "subplot")
  )
str(alpha)
names(alpha)

# dataset is a separate vegetation survey campaign
alpha$dataset <- factor(alpha$dataset)

# Remove NAs
alpha_data <- alpha %>%
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
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
  drop_na()

str(alpha_data)

# Merge alpha & gamma ---------------------------------------------------------

a <- alpha_data %>%
  mutate(scale = "alpha", SR = alpha_10_div, ENSPIE = alpha_10_ENSPIE) %>%
  select(scale, SR, ENSPIE, habitat)

g <- gamma_data %>%
  mutate(scale = "gamma", SR = gamma_100_div, ENSPIE = gamma_100_ENSPIE) %>%
  select(scale, SR, ENSPIE, habitat)

alpha_gamma <- bind_rows(a, g)

str(alpha_gamma)


# Prepare data for beta scale ------------------------------------------------

# Remove NAs and select only needed columns for beta scale

beta_data <- beta_gamma %>%
  dplyr::select(beta_100_div, beta_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    # cover_shrub_total,     inclination,
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

str(beta_data)

# Calculate descriptive statistics -------------------------------------------

# Percentage of Zonal vs azonal habitats ------------------------------------

gamma_data %>%
  count(zonality) %>%
  mutate(pr = n / sum(n) * 100)

# Figure 1 b: Species richness for each scale and habitat type ---------------

# alpha and gamma scale
ggplot(alpha_gamma, aes(habitat, SR, group_by = scale, color = habitat)) +
  geom_point(aes(shape = scale, col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 0.9,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "Species richness", x = 'Grassland habitat type',
    color = "Habitat", shape = "Scale") +
  theme_bw()

# Fig. S7 ---------------------------------------------------------------------

# only gamma diversity
ggplot(gamma_data, aes(habitat, gamma_100_div, group_by = zonality, color = habitat)) +
  geom_point(aes(shape = zonality, col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 1.4,
      jitter.height = 0)) +
  # geom_boxplot(alpha=0, lwd=0.6, outlier.shape = NA)+
  #  stat_boxplot(geom ='errorbar', width = 0.5) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = c(21, 17)) +
  labs(y = expression(paste("Species richness at 100-", m^2, "plots")), 
       x = 'Grassland habitat type',
    color = "Habitat", shape = "Zonality ") +
  theme_bw()

# Fg. S7b
# only beta diversity
ggplot(beta_data, aes(habitat, beta_100_div, color = habitat)) +
  geom_point(aes(col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 1.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  stat_boxplot(geom = 'errorbar', width = 0.5) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "", x = 'Grassland habitat type',
    color = "Habitat") +
  theme_bw()

# Species richness plots for each scale (Fig. S3 a) --------------------------
# SR alpha and gamma scale
ggplot(alpha_gamma, aes(scale, SR, color = scale)) +
  geom_point(aes(fill = scale), col = "black",
    size = 2, alpha = 0.7, pch = 21,
    position = position_jitterdodge(jitter.width = 0.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  # stat_boxplot(geom ='errorbar', width = 0.5) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  scale_color_manual(values = c("#59A5CB", "#D6604D")) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "Species richness", x = 'Scale',
    color = "Scale", fill = "Scale") +
  theme_bw()

# ENSPIE plots for each scale (Fig. S3 b) ------------------------------------
ggplot(alpha_gamma, aes(scale, ENSPIE, color = scale)) +
  geom_point(aes(fill = scale), col = "black",
    size = 2, alpha = 0.7, pch = 21,
    position = position_jitterdodge(jitter.width = 0.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  # stat_boxplot(geom ='errorbar', width = 0.5) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  scale_color_manual(values = c("#59A5CB", "#D6604D")) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "ENSPIE", x = 'Scale',
    color = "Scale", fill = "Scale") +
  theme_bw()

# Plot R2 models  for each scale ----------------------------------------------

# Read all R2 results and combine them in one table --------------------------
R2_alpha_SR <- read_csv("results/partial_R2_M2_alpha_SR.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "alpha", measure = "SR") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_alpha_ENSPIE <- read_csv("results/partial_R2_M2_alpha_ENSPIE.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "alpha", measure = "ENSPIE") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_gamma_SR <- read_csv("results/partial_R2_M2_gamma_SR.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "gamma", measure = "SR") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_gamma_ENSPIE <- read_csv("results/partial_R2_M2_gamma_ENSPIE.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "gamma", measure = "ENSPIE") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

model_R2_all <- R2_alpha_SR %>%
  bind_rows(R2_alpha_ENSPIE) %>%
  bind_rows(R2_gamma_SR) %>%
  bind_rows(R2_gamma_ENSPIE)

model_R2_all

# Plot R2 for SR models (Fig. S3 c) ----------------------------------------

ggplot(model_R2_all %>% filter(measure == "SR"),
  aes(y = Rsq, x = scale, fill = scale), col = "black") +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .2,
    position = position_dodge(.9)) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  # xlab(bquote('Model' ~ m^-2~S^-1'))
  ylab(bquote('Model ' * R^2 * ' (marginal) for species richness')) +
  labs(x = 'Scale', fill = "Scale") +
  theme_bw()

# Plot R2 for ENSPIE models (Fig. S3 d) ------------------------------------

ggplot(model_R2_all %>% filter(measure == "ENSPIE"),
  aes(y = Rsq, x = scale, fill = scale), col = "black") +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .2,
    position = position_dodge(.9)) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  # xlab(bquote('Model' ~ m^-2~S^-1'))
  ylab(bquote('Model ' * R^2 * ' (marginal) for ' * ENS[PIE])) +
  labs(x = 'Scale', fill = "Scale") +
  theme_bw()


# Climate variables ------------------------------------------------------

# Correlation among Temperature and precipitation --------------------------
m1 <- lm(log(Precipt) ~ Temprt, data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m1)
par(mfrow = c(1, 1))

car::Anova(m1)
summary(m1)
gamma_data$Temprt

# Set a theme for all following plots
set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
  axis.title.color = "black", axis.title.size = 1.4, legend.pos = "None", geom.linetype = 2)

plot_model(m1, type = "pred", terms = "Temprt", # show.data=F,
  title = "", line.size = 0.5) + aes(linetype = "solid")

Temp <- get_model_data(m1, type = "pred", terms = "Temprt")
Temp

# todo: What is this and where does it come from
clima_pred_10m$x
clima_pred_10m$predicted
clima_pred_10m$conf.low
clima_pred_10m$conf.high

# Fig. S1 a --------------------------------------------------------------------

fig_temp <- ggplot(Temp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(Temprt, Precipt, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Precipitation, mm", x = 'Temperature') +
  xlim(0, 118) +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_temp

# Correlation among precipitation variability and climate PC -------------------

min(gamma_data$pca1_clima)

# Precipitation variability
m2 <- lm(Prec_Varieb ~ poly(pca1_clima, 2), data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))

car::Anova(m2)
summary(m2)

Prec_Var <- get_model_data(m2, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Prec_Var

# Fig. S1 b ------------------------------------------------------------------

fig_prec_var <- ggplot(Prec_Var, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, Prec_Varieb, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Precipitation variability", x = "Climate gradient (PC))") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_prec_var


# Corrlation among soil pH and climate PC -------------------------------------

m3 <- lm(log(pH) ~ pca1_clima, data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))

Anova(m3)
summary(m3)

ph <- get_model_data(m3, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
ph

# Fig. S2 a -------------------------------------------------------------------

fig_ph <- ggplot(ph, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, pH, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Soil pH", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_ph

# Correlation among soil C and climate PC -------------------------------------

m4 <- lm(log(Corg_percent) ~ poly(pca1_clima, 2), data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m4)
par(mfrow = c(1, 1))

Anova(m4)
summary(m4)

soilS <- get_model_data(m4, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
soilS

# Fig. S2 b -------------------------------------------------------------------

fig_soilC <- ggplot(soilS, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, Corg_percent, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Soil C", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_soilC

# Correlation among litter cover and climate PC ------------------------------

min(gamma_data$cover_litter)

m5 <- lm(log(cover_litter + 1) ~ poly(pca1_clima, 1), data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m5)
par(mfrow = c(1, 1))

car::Anova(m5)
summary(m5)

Litter <- get_model_data(m5, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")
Litter

# Fig. S2 c -------------------------------------------------------------------
fig_litter <- ggplot(Litter, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, cover_litter, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Litter cover", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_litter

# Correlation among soil C and litter cover ----------------------------------
m6 <- lm(log(cover_litter + 1) ~ Corg_percent, data = gamma_data)

# check model assumptions
par(mfrow = c(2, 2))
plot(m6)
par(mfrow = c(1, 1))

car::Anova(m6)
summary(m6)

Corg <- get_model_data(m6, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")
Corg

# Fig. S2 d -------------------------------------------------------------------

fig_Corg <- ggplot(Corg, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(Corg_percent, cover_litter, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Litter cover", x = "Soil C") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

fig_Corg

# Plant cover ---------------------------------------------------------------

# todo: This data does not exist!
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


m7 <- glmer(alpha_10_div ~
  poly(total_cover_10, 2) +
  (1 | dataset / series), family = "poisson", data = alpha_data_cover)
# check model
plot(m7)
qqnorm(resid(m7))
qqline(resid(m7))

sum(residuals(m7, type = "pearson")^2) / df.residual(m7)
check_overdispersion(m7)


Anova(m7)
summary(m7)


m8 <- glmer(gamma_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), family = "poisson", data = gamma_data_cover)

# check model
plot(m8)
qqnorm(resid(m8))
qqline(resid(m8))

check_overdispersion(m8)
# use negative binomial

m8b <- glmer.nb(gamma_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), data = gamma_data_cover)

check_overdispersion(m8b)

car::Anova(m8b)
summary(m8b)


beta_data_cover <- beta_data %>%
  left_join(tot_cover_100, by = c("dataset", "series"))

m9 <- lmer(beta_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), data = beta_data_cover)


# check model
plot(m9)
qqnorm(resid(m9))
qqline(resid(m9))

car::Anova(m9)
summary(m9)



# plots
min(alpha_data_cover$total_cover_10)
max(alpha_data_cover$total_cover_10)

alphaSR_PlantCover <- get_model_data(m7, type = "pred", terms = "total_cover_10[7.122:221, by=.001]")
alphaSR_PlantCover

# Fig. S4 a -----
Fig.alphaSR_PlantCover <- ggplot(alphaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = alpha_data_cover, aes(total_cover_10, alpha_10_div, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#64ABCE")

Fig.alphaSR_PlantCover



min(gamma_data_cover$total_cover_100)
max(gamma_data_cover$total_cover_100)

gammaSR_PlantCover <- get_model_data(m8b, type = "pred", terms = "total_cover_100[7.122:221, by=.001]")
gammaSR_PlantCover

# Fig. S4 b -----

Fig.gammaSR_PlantCover <- ggplot(gammaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data_cover, aes(total_cover_100, gamma_100_div,
    fill = habitat, col = habitat),
  size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#D6604D")

Fig.gammaSR_PlantCover



m9 <- lmer(beta_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), data = beta_data_cover)


betaSR_PlantCover <- get_model_data(m9, type = "pred", terms = "total_cover_100[7.122:221, by=.001]")
betaSR_PlantCover

# Fig. S4 c -----

Fig.betaSR_PlantCover <- ggplot(betaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data_cover, aes(total_cover_100, beta_100_div,
    fill = habitat, col = habitat),
  size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaSR_PlantCover



############################# End
