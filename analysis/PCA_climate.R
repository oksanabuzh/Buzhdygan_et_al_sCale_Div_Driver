# Purpose: PCA for climate gradient

library(tidyverse)
library(factoextra)

header_data <- read_csv("data/headers.csv")

# Climate -----------------------------------------------------------------

climate <- header_data %>%
  group_by(series) %>%
  summarise(BIO1 = mean(BIO1), BIO12 = mean(BIO12)) %>%
  mutate(
    Temprt = as.vector(scale(BIO1, scale = TRUE)),
    Precipt = as.vector(scale(BIO12, scale = TRUE))
  ) %>%
  select(Temprt, Precipt, series) %>%
  drop_na()

nrow(climate)

climate_pca <- prcomp(climate %>% select(-series))

climate_pca
# round(climate_pca$rotation[, 1], 2)

summary(climate_pca)

fviz_eig(climate_pca
  , addlabels = TRUE
)

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

pca1_clima <- climate_pca_ind$coord[, 1]

climate <- bind_cols(climate, pca1_clima = pca1_clima)


head(climate_pca_var$contrib)
head(climate_pca_var$cor)


library("corrplot")
corrplot(climate_pca_var$cor, is.corr = F)

library(ggcorrplot)
ggcorrplot(climate_pca_var$contrib,
  hc.order = TRUE, type = "lower",
  lab = TRUE, lab_size = 3,
  colors = c("red", "white", "blue")
)

climate_pca_var$cor



# Soil texture ------------------------------------------------------------

# Create a subset of the texture data and calculate the mean values
# for each plot
soil_text <- header_data %>%
  mutate(
    sand = as.vector(scale(sand, scale = T)),
    silt = as.vector(scale(silt, scale = T)),
    clay = as.vector(scale(clay, scale = T)),
    #
    pH = as.vector(scale(pH, scale = T)),
    C = as.vector(scale(Corg_percent, scale = T))
  ) %>%
  dplyr::select(
    sand, silt, clay,
    pH, C,
    plotID
  ) %>%
  drop_na()

# Run the PCA and look at the results
soil_text_pca <- prcomp(soil_text %>% select(-plotID))

soil_text_pca
# round(soil_text_pca$rotation[, 1], 2)

summary(soil_text_pca)

# Plots
fviz_eig(soil_text_pca)

fviz_pca_var(soil_text_pca,
  axes = c(1, 2),
  # col.var = "cor", # Color by contributions to the PC
  # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE # Avoid text overlapping
)

fviz_pca_var(soil_text_pca,
  axes = c(2, 1),
  col.var = "contrib", # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "yellow", "red", "blue"),
  repel = TRUE # Avoid text overlapping
)

soil_text_var <- get_pca_var(soil_text_pca)
soil_text_var$contrib
soil_text_var$cor
soil_text_ind <- get_pca_ind(soil_text_pca)
soil_text_ind
soil_text_ind$coord

round(soil_text_var$cor, 2)
head(soil_text_var$contrib)
head(soil_text_var$cor)


library("corrplot")
corrplot(soil_text_var$cor, is.corr = F)
corrplot(soil_text_var$contrib, is.corr = F)

col <- colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))

corrplot(soil_text_var$contrib, is.corr = F, method = "color",
  col = col(6),
  # type="lower", # order="hclust",
  addCoef.col = "black", number.cex = 0.7, # Add coefficient of correlation
  tl.col = "black", tl.srt = 45 # , #Text label color and rotation
  # Combine with significance
  # p.mat = p.mat, sig.level = 0.01, insig = "blank",
  # hide correlation coefficient on the principal diagonal
  #  diag=FALSE
)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(soil_text_var$cor, is.corr = F, method = "color",
  col = col(200),
  # type="lower", # order="hclust",
  addCoef.col = "black", number.cex = 0.7, # Add coefficient of correlation
  tl.col = "black", tl.srt = 45 # , #Text label color and rotation
  # Combine with significance
  # p.mat = p.mat, sig.level = 0.01, insig = "blank",
  # hide correlation coefficient on the principal diagonal
  #  diag=FALSE
)


soil_text_var$cor

library(ggcorrplot)
ggcorrplot(soil_text_var$contrib,
  # hc.order = TRUE,
  type = "lower",
  lab = TRUE, lab_size = 3)

# Extract the PC axis coordinates and bind them to the dataset
pca1_soil_text <- soil_text_ind$coord[, 1]
pca2_soil_text <- soil_text_ind$coord[, 2]

soil_text <- bind_cols(soil_text,
  pca1_soil_text = pca1_soil_text,
  pca2_soil_text = pca2_soil_text
)




# Combine all PCAs --------------------------------------------------------

headers_PCA <- full_join(header_data, climate, by = "series")
headers_PCA

PCA_all <- headers_PCA %>%
  select(pca1_soil_text, pca2_soil_text, pca1_clima)

corl1 <- round(cor(PCA_all,
  method = c("pearson"), use = "pairwise.complete.obs"
), 2)

corl1

library(ggcorrplot)
ggcorrplot(corl1,
  hc.order = TRUE, type = "lower",
  lab = TRUE, lab_size = 3,
  colors = c("red", "white", "blue")
)


write_csv(climate, "data/climate_PCA.csv")
