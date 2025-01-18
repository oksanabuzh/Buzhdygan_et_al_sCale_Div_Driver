# Purpose: Tests of Residual Spatial Correlation

# helpful information:
# https://theoreticalecology.wordpress.com/2012/05/12/spatial-autocorrelation-in-statistical-models-friend-or-foe/
# https://rdrr.io/cran/DHARMa/man/testSpatialAutocorrelation.html

# libraries----
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)


# prepare data----

# "data/climate_PCA.csv" contains scores for the compound climate variable, 
# derived from the PCA analysis in "1_prepare_data/ PCA_environment.R"

# "data/Environm_variabl.csv" contains all environmental data


climate_PCA <- read.csv("data/climate_PCA.csv")

header <- read_csv("data/Environm_variabl.csv") %>% 
  full_join(
    read.csv("data/climate_PCA.csv"),
    by = "series"
  )

str(header) 
names (header)

# prepare subset of data for alpha scale (10-m2 plots)

alpha <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>%
  filter(type=="alpha")%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header, 
            by=c("dataset", "plotID", "series", "subplot")
  ) %>% 
  mutate(dataset=factor(dataset))

str(alpha) 
names (alpha)

# dataset is a separate vegetation survey campaign
alpha$dataset 


# Select variables, remove NAs  
alpha_data <- alpha %>% 
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,  
                lat, lon, pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad,
                subplot) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  
  drop_na

str(alpha_data)

# plot on a mean alpha per series to omit pseudoreplication of the plots:
alpha_mean <- alpha_data %>% 
  dplyr::select(alpha_10_div, alpha_10_ENSPIE,  
                lat, lon, pca1_clima, 
                grazing_intencity, mowing, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

str(alpha_mean)

## data gamma----

header_mean <- header %>% 
  select(c(series, lat, lon, zonality, habitat_broad, 
           where(is.numeric))) %>% 
  group_by(series,  zonality, habitat_broad) %>% 
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  %>% 
  ungroup()

beta_gamma <-read_csv("data/alpha_beta_gamma_community_variabl.csv") %>% 
  filter(type=="gamma" | type=="beta" )%>% 
  unite("metric", c(type, scale, metric), sep="_") %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  full_join(header_mean, by=c("dataset", "series") )%>% 
  mutate(dataset=factor(dataset))

str(beta_gamma) 
names (beta_gamma)

# dataset is a separate vegetation survey campaign
beta_gamma$dataset 

# Select variables, remove NAs 

gamma_data <- beta_gamma %>% 
  dplyr::select(gamma_100_div, gamma_100_ENSPIE, 
                lat, lon, pca1_clima, 
                grazing_intencity, mowing, 
                # cover_shrub_total,     inclination, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

## data beta----
beta_data <- beta_gamma %>% 
  dplyr::select(beta_100_div, beta_100_ENSPIE, 
                lat, lon, pca1_clima, 
                grazing_intencity, mowing, 
                # cover_shrub_total,     inclination, 
                cover_litter,
                BIO7, BIO15,
                pH, Corg_percent,
                dataset, series, habitat_broad, zonality) %>% 
  mutate(Tem_range = BIO7,
         Prec_Varieb = BIO15,
         mowing=factor(mowing)) %>% 
  mutate(habitat =fct_relevel(habitat_broad, c("saline", "complex", "dry", 
                                               "wet" , "mesic", "fringe", "alpine"))) %>% 
  drop_na

str(beta_data)




# Morans's I tests ----  
#------------------------------------------------------------------------------#
# alpha ----
#------------------------------------------------------------------------------#
# alpha SR ----
#------------------------------------------------------------------------------#
m1_3_a <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)


m2_1_a <- glmer (alpha_10_div ~ 
                 poly(pca1_clima, 2) +
                 poly(Prec_Varieb, 2) +
                 poly(Corg_percent,2)+
                 pH +
                 poly(cover_litter,2) +
                 grazing_intencity + mowing +
                 (1|dataset/series), family = "poisson", data = alpha_data)


## Get residuals ----

#  For lme4, re.form = NULL simulats residuals conditional on fitted rendom effcets
## re.form specify which random effects to condition on when predicting. 
# If NULL, include all random effects; if NA or ~0, include no random effects.
# https://rdrr.io/cran/DHARMa/man/testSpatialAutocorrelation.html

# (1) get randomized residuals.
res.sim_alpha_m1 <- DHARMa::simulateResiduals(m1_3_a, re.form = NULL)
res.sim_alpha_m2 <- DHARMa::simulateResiduals(m2_1_a, re.form = NULL)

# randomized residuals from DHARMa is better than deviance residuals, because  deviance residuals are not homogeneous)
# https://stats.stackexchange.com/questions/507934/testing-the-spatiale-autocorrelation-on-the-residuals-of-the-mixed-effect-logist?newreg=e8a0041e387743139c3e9885b71d62eb
# https://rdrr.io/cran/DHARMa/man/testSpatialAutocorrelation.html
residuals(m1_3_a, type="deviance")

# (2)  generate a matrix of inverse distance weights.  
# In the matrix, entries for pairs of points that are close together 
# are higher than for pairs of points that are far apart.  
# We use latitude and longitude for each plot, generate a distance matrix, 
# then take inverse of the matrix values and replace the diagonal entries with zero:
dM_a =  1/as.matrix(dist(cbind(alpha_data$lon, alpha_data$lat)))
diag(dM_a) <- 0
str(dM_a)
# We have created a matrix where each off-diagonal entry [i, j] in the matrix is equal to 1/(distance between point i and point j). 


# (3) calculate Moran’s I (DHARMa works using ape package)
DHARMa::testSpatialAutocorrelation(res.sim_alpha_m1, distMat = dM_a)
DHARMa::testSpatialAutocorrelation(res.sim_alpha_m2, distMat = dM_a)

#------------------------------------------------------------------------------#
# alpha ENSPIE----
#------------------------------------------------------------------------------#
m1_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)

m2_1_ENSPIE <- lmer(log(alpha_10_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb, 2) +
                      pH +
                      Corg_percent +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset/series),  data = alpha_data)


## Get residuals ----

# (1) get randomized residuals.
res.sim_alpha_m1_ENSPIE <- DHARMa::simulateResiduals(m1_1_ENSPIE, re.form = NULL)
res.sim_alpha_m2_ENSPIE <- DHARMa::simulateResiduals(m2_1_ENSPIE, re.form = NULL)



# (2)  generate a matrix of inverse distance weights.  
# In the matrix, entries for pairs of points that are close together 
# are higher than for pairs of points that are far apart.  
# We use latitude and longitude for each plot, generate a distance matrix, 
# then take inverse of the matrix values and replace the diagonal entries with zero:
dM_a =  1/as.matrix(dist(cbind(alpha_data$lon, alpha_data$lat)))
diag(dM_a) <- 0
str(dM_a)
# We have created a matrix where each off-diagonal entry [i, j] in the matrix is equal to 1/(distance between point i and point j). 

# (3) calculate Moran’s I (DHARMa works using ape package)
DHARMa::testSpatialAutocorrelation(res.sim_alpha_m1_ENSPIE, distMat = dM_a)
DHARMa::testSpatialAutocorrelation(res.sim_alpha_m2_ENSPIE, distMat = dM_a)

#------------------------------------------------------------------------------#
###  gamma -----
#------------------------------------------------------------------------------#
# gamma SR ----
#------------------------------------------------------------------------------#
m1_1_g <- glmer.nb(gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

m2_1_g <- glmer.nb(gamma_100_div ~ 
                    poly(pca1_clima, 2) +
                    poly(Prec_Varieb, 2) +
                    poly(Corg_percent,2)+
                    poly(pH, 2) +
                    poly(cover_litter,2) +
                    grazing_intencity + mowing + 
                    (1|dataset), data = gamma_data)

# (1) get randomized residuals.
res.sim_gamma_m1 <- DHARMa::simulateResiduals(m1_1_g, re.form = NULL)
res.sim_gamma_m2 <- DHARMa::simulateResiduals(m2_1_g, re.form = NULL)

# (2)  generate a matrix of inverse distance weights. 
dM_g =  1/as.matrix(dist(cbind(gamma_data$lon, gamma_data$lat)))
diag(dM_g) <- 0
str(dM_g)

# (3) Moran’s I
DHARMa::testSpatialAutocorrelation(res.sim_gamma_m1, distMat = dM_g)
DHARMa::testSpatialAutocorrelation(res.sim_gamma_m2, distMat = dM_g)


#------------------------------------------------------------------------------#
# gamma ENSPIE ----
#------------------------------------------------------------------------------#
m1_3_ENSPIE_g <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)

m2_1_ENSPIE_g <- lmer(log(gamma_100_ENSPIE) ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = gamma_data)


# (1) get randomized residuals.
res.sim_gamma_m1_ENSPIE <- DHARMa::simulateResiduals(m1_3_ENSPIE_g, re.form = NULL)
res.sim_gamma_m2_ENSPIE <- DHARMa::simulateResiduals(m2_1_ENSPIE_g, re.form = NULL)

# (2)  generate a matrix of inverse distance weights. 
dM_g =  1/as.matrix(dist(cbind(gamma_data$lon, gamma_data$lat)))
diag(dM_g) <- 0
str(dM_g)

# (3) Moran’s I
DHARMa::testSpatialAutocorrelation(res.sim_gamma_m1_ENSPIE, distMat = dM_g)
DHARMa::testSpatialAutocorrelation(res.sim_gamma_m2_ENSPIE, distMat = dM_g)


#------------------------------------------------------------------------------#
# beta ----
#------------------------------------------------------------------------------#
# beta SR ----
#------------------------------------------------------------------------------#

m1_1_b <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)


m2_1_b <- lmer (beta_100_div ~ 
                poly(pca1_clima, 2) +
                poly(Prec_Varieb, 2) +
                poly(Corg_percent,2)+
                poly(pH, 2) +
                poly(cover_litter,2) +
                grazing_intencity + mowing +
                (1|dataset),  data = beta_data)
## Get residuals ----

# (1) get randomized residuals.
res.sim_beta_m1 <- DHARMa::simulateResiduals(m1_1_b, re.form = NULL)
res.sim_beta_m2 <- DHARMa::simulateResiduals(m2_1_b, re.form = NULL)

# (2)  generate a matrix of inverse distance weights. 
dM_b =  1/as.matrix(dist(cbind(beta_data$lon, beta_data$lat)))
diag(dM_b) <- 0
str(dM_b)

# (3) Moran’s I
DHARMa::testSpatialAutocorrelation(res.sim_beta_m1, distMat = dM_b)
DHARMa::testSpatialAutocorrelation(res.sim_beta_m2, distMat = dM_b)


#------------------------------------------------------------------------------#
# beta ENSPIE ----
#------------------------------------------------------------------------------#

m1_1_ENSPIE_b <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(pH, 2) +
                      poly(Corg_percent,2) +
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)
m2_1_ENSPIE_b <- lmer(beta_100_ENSPIE ~ 
                      poly(pca1_clima, 2) +
                      poly(Prec_Varieb,2) +
                      poly(pH,2) +
                      poly(Corg_percent,2)+
                      poly(cover_litter,2) +
                      grazing_intencity + mowing +
                      (1|dataset),  data = beta_data)

## Get residuals ----
# (1) get randomized residuals.
res.sim_beta_m1_ENSPIE <- DHARMa::simulateResiduals(m1_1_ENSPIE_b, re.form = NULL)
res.sim_beta_m2_ENSPIE <- DHARMa::simulateResiduals(m2_1_ENSPIE_b, re.form = NULL)

# (2)  generate a matrix of inverse distance weights. 
dM_b =  1/as.matrix(dist(cbind(beta_data$lon, beta_data$lat)))
diag(dM_b) <- 0
str(dM_b)

# (3) Moran’s I
DHARMa::testSpatialAutocorrelation(res.sim_beta_m1_ENSPIE, distMat = dM_b)
DHARMa::testSpatialAutocorrelation(res.sim_beta_m2_ENSPIE, distMat = dM_b)
