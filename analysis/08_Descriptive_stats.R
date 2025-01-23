# Purpose: Descriptive and summary statistics,
#          generate Fig. 1 b,c;  Fig. S1 a-b; Fig. S2 a-d;
#                   Fig. S3 a,b,c,d; Fig. S4 a,b; Fig. S7.

# Load libraries --------------------------------------------------------------
library(tidyverse)
library(sjPlot)
library(lme4)
library(performance)

# Define habitat colors to be used in all plots to distinguish habitats
habitat_colors = c(
  saline = "#4e3910",
  complex = "#CC6600",
  dry = "#e3c28b",
  wet = "#CC99FF",
  mesic = "#0066FF",
  fringe = "#00B200",
  alpine = "#006600")

# Set theme for the model plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1,
          axis.textsize.y = 1,
          axis.textcolor = "black",
          axis.title.color = "black",
          axis.title.size = 1.4,
          legend.pos = "None",
          geom.linetype = 2)


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

# Prepare subset of data for gamma scale (100 m2 plots) -------------------------
beta_gamma <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "gamma" | type == "beta") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header_mean, by = c("dataset", "series")) %>% 
  mutate(dataset = factor(dataset))

# selected variables and remove NAs
gamma_data <- beta_gamma %>%
  dplyr::select(gamma_100_div, gamma_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    Tem_range, Prec_Varieb, Temprt, Precipt,
    pH, Corg_percent,
    dataset, series, habitat_broad, zonality) %>%
  mutate(mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad,
    c("saline", "complex", "dry",
      "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

# Prepare subset of data for alpha scale (10 m2 plots) -------------------------
alpha <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type == "alpha") %>%
  unite("metric", c(type, scale, metric), sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  full_join(header,
    by = c("dataset", "plotID", "series", "subplot")
  ) %>% 
  mutate(dataset = factor(dataset))

# selected variables and remove NAs
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
  mutate(habitat = fct_relevel(habitat_broad,
    c("saline", "complex", "dry",
      "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

# Combine alpha and gamma data ------------------------------------------------

a <- alpha_data %>%
  # pivot_longer(c("alpha_10_div", "alpha_10_ENSPIE"), names_to = "Diver") %>%
  mutate(scale = "alpha", SR = alpha_10_div, ENSPIE = alpha_10_ENSPIE) %>%
  select(scale, SR, ENSPIE, habitat)

str(a)

g <- gamma_data %>%
  # pivot_longer(c("alpha_10_div", "alpha_10_ENSPIE"), names_to = "Diver") %>%
  mutate(scale = "gamma", SR = gamma_100_div, ENSPIE = gamma_100_ENSPIE) %>%
  select(scale, SR, ENSPIE, habitat)

alpha_gamma <- a %>%
  bind_rows(g)

str(alpha_gamma)

# Prepare subset of data for beta scale --------------------------------------

beta_data <- beta_gamma %>%
  dplyr::select(beta_100_div, beta_100_ENSPIE,
    pca1_clima,
    grazing_intencity, mowing,
    cover_litter,
    Tem_range, Prec_Varieb,
    pH, Corg_percent,
    dataset, series, habitat_broad) %>%
  mutate(mowing = factor(mowing)) %>%
  mutate(habitat = fct_relevel(habitat_broad,
    c("saline", "complex", "dry",
      "wet", "mesic", "fringe", "alpine"))) %>%
  drop_na()

# Calculate descriptive statistics ---------------------------------------------

# % of zonal vs azonal habitats -----------------------------------------------
gamma_data %>%
  count(zonality) %>%
  mutate(pr = n / sum(n) * 100)

# Plot SR for each habitat ----------------------------------------------------
# Fig. 1 b -------------------------------------------------------------------
ggplot(alpha_gamma, aes(habitat, SR, group_by = scale, color = habitat)) +
  geom_point(aes(shape = scale, col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 0.9,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  #  stat_boxplot(geom ='errorbar', width = 0.5) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "Species richness", x = 'Grassland habitat type',
    color = "Habitat", shape = "Scale") +
  theme_bw()

# Fig. 1 c -------------------------------------------------------------------
ggplot(beta_data, aes(habitat, beta_100_div, color = habitat)) +
  geom_point(aes(col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 1.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  stat_boxplot(geom = 'errorbar', width = 0.5) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "", x = 'Grassland habitat type',
    color = "Habitat") +
  theme_bw()

# Fig. S7 --------------------------------------------------------------------

# Zonal/Azonal
ggplot(gamma_data, aes(habitat, gamma_100_div, group_by = zonality, color = habitat)) +
  geom_point(aes(shape = zonality, col = habitat),
    size = 2, alpha = 0.9,
    position = position_jitterdodge(jitter.width = 1.4,
      jitter.height = 0)) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = c(21, 17)) +
  labs(y = expression(paste("Species richness at 100-", m^2, "plots")), 
       x = 'Grassland habitat type',
    color = "Habitat", shape = "Zonality ") +
  theme_bw()


# Plot SR for each scale ------------------------------------------------------
# Fig. S3 a -------------------------------------------------------------------
ggplot(alpha_gamma, aes(scale, SR, color = scale)) +
  geom_point(aes(fill = scale), col = "black",
    size = 2, alpha = 0.7, pch = 21,
    position = position_jitterdodge(jitter.width = 0.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  scale_color_manual(values = c("#59A5CB", "#D6604D")) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "Species richness", x = 'Scale',
    color = "Scale", fill = "Scale") +
  theme_bw()

# Plot ENSPIE for each scale --------------------------------------------------
# Fig. S3 b -------------------------------------------------------------------
ggplot(alpha_gamma, aes(scale, ENSPIE, color = scale)) +
  geom_point(aes(fill = scale), col = "black",
    size = 2, alpha = 0.7, pch = 21,
    position = position_jitterdodge(jitter.width = 0.4,
      jitter.height = 0)) +
  geom_boxplot(alpha = 0, lwd = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  scale_color_manual(values = c("#59A5CB", "#D6604D")) +
  scale_shape_manual(values = c(21, 19)) +
  labs(y = "ENSPIE", x = 'Scale',
    color = "Scale", fill = "Scale") +
  theme_bw()

# plot model R2 for each scale-------------------------------------------------

# Fig. S3 c SR ----------------------------------------------------------------
R2_alpha_SR <- read_csv("results/R2_alpha_SR.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "alpha", measure = "SR") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_alpha_ENSPIE <- read_csv("results/R2_alpha_ENSPIE.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "alpha", measure = "ENSPIE") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_gamma_SR <- read_csv("results/R2_gamma_SR.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "gamma", measure = "SR") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

R2_gamma_ENSPIE <- read_csv("results/R2_gamma_ENSPIE.csv") %>%
  filter(Effect == "Model") %>%
  mutate(scale = "gamma", measure = "ENSPIE") %>%
  select(scale, measure, Rsq, upper.CL, lower.CL)

# Combine all R2 models
model_R2_all <- R2_alpha_SR %>%
  bind_rows(R2_alpha_ENSPIE) %>%
  bind_rows(R2_gamma_SR) %>%
  bind_rows(R2_gamma_ENSPIE)

model_R2_all

ggplot(model_R2_all %>% filter(measure == "SR"),
  aes(y = Rsq, x = scale, fill = scale), col = "black") +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .2,
    position = position_dodge(.9)) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  ylab(bquote('Model ' * R^2 * ' (marginal) for species richness')) +
  labs(x = 'Scale', fill = "Scale") +
  theme_bw()

# Fig. S3 d ENSPIE ------------------------------------------------------------

ggplot(model_R2_all %>% filter(measure == "ENSPIE"),
  aes(y = Rsq, x = scale, fill = scale), col = "black") +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .2,
    position = position_dodge(.9)) +
  scale_fill_manual(values = c("#59A5CB", "#D6604D")) +
  ylab(bquote('Model ' * R^2 * ' (marginal) for ' * ENS[PIE])) +
  labs(x = 'Scale', fill = "Scale") +
  theme_bw()


# Climate variables -----------------------------------------------------------
# Correlation among Temperature and precipitation -----------------------------

m1 <- lm(log(Precipt) ~ Temprt, data = gamma_data)

# check model
par(mfrow = c(2, 2))
plot(m1)
par(mfrow = c(1, 1))

car::Anova(m1)
summary(m1)

# Fig. S1 a Temperature ------------------------------------------------------
Temp <- get_model_data(m1, type = "pred", terms = "Temprt")

Fig.Temp <- ggplot(Temp, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(Temprt, Precipt, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Precipitation, mm", x = 'Temperature') +
  xlim(0, 118) +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.Temp

# Correlation among Climate PC and precipitation variability -----------------

m2 <- lm(Prec_Varieb ~ poly(pca1_clima, 2), data = gamma_data)

# check model
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))

car::Anova(m2)
summary(m2)

# Fig. S1 b Precipitation variability ----------------------------------------
Prec_Var <- get_model_data(m2, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.Prec_Var <- ggplot(Prec_Var, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, Prec_Varieb, fill = habitat, 
                                    col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Precipitation variability", x = "Climate gradient (PC))") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.Prec_Var

# soil pH ---------------------------------------------------------------------

m3 <- lm(log(pH) ~ pca1_clima, data = gamma_data)

# check model

par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))

car::Anova(m3)
summary(m3)

# Fig. S2 a Soil pH ----------------------------------------------------------
ph <- get_model_data(m3, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.ph <- ggplot(ph, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, pH, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Soil pH", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.ph

# soil C ----------------------------------------------------------------------

m4 <- lm(log(Corg_percent) ~ poly(pca1_clima, 2), data = gamma_data)

# check model
par(mfrow = c(2, 2))
plot(m4)
par(mfrow = c(1, 1))

car::Anova(m4)
summary(m4)

# Fig. S2 b soil PCarbon -----------------------------------------------------
soilC <- get_model_data(m4, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.soilC <- ggplot(soilC, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, Corg_percent, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + scale_color_manual(values = habitat_colors) +
  labs(y = "Soil C", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.soilC


# Litter % --------------------------------------------------------------------

m5 <- lm(log(cover_litter + 1) ~ poly(pca1_clima, 1), data = gamma_data)

# check model
par(mfrow = c(2, 2))
plot(m5)
par(mfrow = c(1, 1))

car::Anova(m5)
summary(m5)

# Fig. S2 c Litter cover -----------------------------------------------------
Litter <- get_model_data(m5, type = "pred", terms = "pca1_clima[-1.2:4.8, by=.001]")

Fig.Litter <- ggplot(Litter, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(pca1_clima, cover_litter,
                                    fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Litter cover", x = "Climate gradient (PC)") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.Litter

# Soil C % and litter cover ---------------------------------------------------
m6 <- lm(log(cover_litter + 1) ~ Corg_percent, data = gamma_data)

# check model

par(mfrow = c(2, 2))
plot(m6)
par(mfrow = c(1, 1))

car::Anova(m6)
summary(m6)

# Fig. S2 d Soil carbon and litter -------------------------------------------
Corg <- get_model_data(m6, type = "pred", terms = "Corg_percent[0:9.5, by=.001]")

Fig.Corg <- ggplot(Corg, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data, aes(Corg_percent, cover_litter, 
                                    fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Litter cover", x = "Soil C") +
  geom_line(linetype = 1, linewidth = 1, col = "black")

Fig.Corg


# Plant cover -----------------------------------------------------------------

# alpha scale
tot_cover_10 <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(scale == 10 &  metric == "cover") %>%
  rename(total_cover_10 = value) %>%
  select(dataset, series, subplot, total_cover_10) %>%
  mutate(dataset = factor(dataset))

alpha_data_cover <- alpha_data %>%
  left_join(tot_cover_10, by = c("dataset", "series", "subplot"))

# gamma scale
tot_cover_100 <- read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(scale == 100 &  metric == "cover") %>%
  rename(total_cover_100 = value) %>%
  select(dataset, series, subplot, total_cover_100) %>%
  mutate(dataset = factor(dataset))

gamma_data_cover <- gamma_data %>%
  left_join(tot_cover_100, by = c("dataset", "series"))

# alpha model
m7 <- glmer(alpha_10_div ~
  poly(total_cover_10, 2) +
  (1 | dataset / series), family = "poisson", data = alpha_data_cover)

# check model
plot(m7)
qqnorm(resid(m7))
qqline(resid(m7))

# Check overdispersion
sum(residuals(m7, type = "pearson")^2) / df.residual(m7)
check_overdispersion(m7)

car::Anova(m7)
summary(m7)

# gamma model
m8 <- glmer(gamma_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), family = "poisson", data = gamma_data_cover)

# check model
plot(m8)
qqnorm(resid(m8))
qqline(resid(m8))

check_overdispersion(m8)

# use negative binomial because of overdispersion
m8b <- glmer.nb(gamma_100_div ~
  poly(total_cover_100, 2) +
  (1 | dataset), data = gamma_data_cover)

check_overdispersion(m8b)

car::Anova(m8b)
summary(m8b)

# beta model
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


# Cover plots -----------------------------------------------------------------
alphaSR_PlantCover <- get_model_data(m7, type = "pred", terms = "total_cover_10[7.122:221, by=.001]")

# Fig. S4 a Cover alpha - richness -------------------------------------------
Fig.alphaSR_PlantCover <- ggplot(alphaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = alpha_data_cover, aes(total_cover_10, alpha_10_div, fill = habitat, col = habitat),
    size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#64ABCE")

Fig.alphaSR_PlantCover


# Fig. S4 b plant cover gamma ------------------------------------------------
gammaSR_PlantCover <- get_model_data(m8b, type = "pred", terms = "total_cover_100[7.122:221, by=.001]")

Fig.gammaSR_PlantCover <- ggplot(gammaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = gamma_data_cover, aes(total_cover_100, gamma_100_div,
    fill = habitat, col = habitat),
  size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#D6604D")

Fig.gammaSR_PlantCover


# Fig. S4 c plant cover beta -------------------------------------------------
betaSR_PlantCover <- get_model_data(m9, type = "pred", terms = "total_cover_100[7.122:221, by=.001]")

Fig.betaSR_PlantCover <- ggplot(betaSR_PlantCover, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(data = beta_data_cover, aes(total_cover_100, beta_100_div,
    fill = habitat, col = habitat),
  size = 3, alpha = 0.7, pch = 21) +
  scale_fill_manual(values = habitat_colors) + 
  scale_color_manual(values = habitat_colors) +
  labs(y = "Species richness", x = 'Plant cover') +
  geom_line(linetype = 1, linewidth = 1, col = "#00AC7F")

Fig.betaSR_PlantCover
