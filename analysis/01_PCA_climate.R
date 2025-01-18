# Purpose: PCA analysis for extracting compaud variable of climate gradient

# packages -----
library(tidyverse)
library(factoextra)

# data----
header_data <- read_csv("data/Environm_variabl.csv")

# selected climate variables
climate <- header_data %>%
  group_by(series) %>% 
  summarise(BIO1=mean(BIO1), BIO12=mean(BIO12)) %>% 
  mutate(
    Temprt = as.vector(scale(BIO1, scale = TRUE)),
    Precipt = as.vector(scale(BIO12, scale = TRUE))
  ) %>%
  select(Temprt, Precipt, series) %>%
  drop_na()

nrow(climate)

# perform PCA
climate_pca <- prcomp(climate %>% select(-series))

climate_pca
# round(climate_pca$rotation[, 1], 2)

summary(climate_pca)

# Scree plot
fviz_eig(climate_pca, 
         addlabels = TRUE
         )

# 2-dimension plot of PCA coordinates
fviz_pca_var(climate_pca,
 # axes = c(1, 2),
  col.var = "contrib", # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE # Avoid text overlapping
)



fviz_pca_var(climate_pca,
  axes = c(1, 1),
  col.var = "contrib", # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE # Avoid text overlapping
)

climate_pca_var <- get_pca_var(climate_pca)
climate_pca_ind <- get_pca_ind(climate_pca)
climate_pca_ind$coord

pca1_clima <- climate_pca_ind$coord[, 1] # select 1st PC

# join with raw data 
climate <- bind_cols(climate, pca1_clima = pca1_clima)


head(climate_pca_var$contrib)
head(climate_pca_var$cor)

# Plot correlations of variables with the PC1 and PC2
library("corrplot")
corrplot(climate_pca_var$cor, is.corr=F) 

## Plot contribution of variables to the PC1 and PC2
library(ggcorrplot)
ggcorrplot(climate_pca_var$contrib,
           hc.order = TRUE, # type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue")
)

climate_pca_var$cor



write_csv(climate, "data/climate_PCA.csv")
