# Folder "data"

## About

Contains data that are used for analysis 

#### Structure

| Document                               | What's inside                            |
| -------------------------------------- |----------------------------------------- |
| Environm_variab.csv                    | Environmental data                       |
| alpha_beta_gamma_community_variabl.csv | Plant diversity measures and plant cover for each scale |
| climate_PCA                            | scores from the [PCA analysis](analysis/PCA_climate.R)  |       


#### Metadata

##### Environm_variab.csv 

dataset

plotID

series

subplot

habitat_broad

zonality

BIO1

BIO12

BIO7

BIO15

grazing_intencity

mowing

cover_litter

pH

Corg_percent

#### alpha_beta_gamma_community_variabl.csv
Combines all diversity measures and plant cover for each scale:
- alpha diversity and cover measures (SR, ENSPIE, and cover) include doubled 10 m2 plots
- gamma diversity and cover measures (SR, ENSPIE, and cover) include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
- beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

SR - species richness

ENSPIE - evenness measure calculated as inverse Simpson using species cover

cover - cumulative cover of plant community




##### climate_PCA