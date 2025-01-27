# Folder "data"

## About

Contains data that are used for analysis 

### Structure

| Document                               | Description                                                  |
| -------------------------------------- | ------------------------------------------------------------ |
| Environm_variab.csv                    | Environmental data                                           |
| alpha_beta_gamma_community_variabl.csv | Plant diversity measures and plant cover for each scale      |
| climate_PCA.csv                        | scores from the [PCA analysis](../analysis/01_PCA_climate.R) |
| aggregation.csv                        | intraspecific aggregation proxy for each plot


## Metadata

### `Environm_variab.csv`

This dataset contains all measured environmental data, that were used in the paper

#### Variables

| Short name        | type      | Long name                 | Description                                                                                                                                                                                                                        |
| ----------------- | --------- | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| dataset           | numeric   | Dataset ID                | Dataset ID, sampled by different teams and years                                                                                                                                                                                   |
| series            | character | Series ID                 | 100 m2 plot that includes two 10 m2 plots that are nested within it                                                                                                                                                                |
| subplot           | character | Subplot                   | one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner                                                                                               |
| plotID            | character | PlotID                    | Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots                                                                                                                      |
| habitat_broad     | character | Grassland habitat type    | Broad types of grassland habitats, distinguished based on EUNIS habitat type system                                                                                                                                                |
| zonality          | character | Zonality of vegetation    | Zonality:  zonal (natural) vegetation types, azonal (seminatural) grassland types                                                                                                                                                  |
| Temprt            | numeric   | Mean annual temperature   | Mean annual temperature, extracted from the CHELSA climate database using plot coordinates                                                                                                                                         |
| Precipt           | numeric   | Mean annual precipitation | Mean annual precipitation, extracted from the CHELSA climate database using plot coordinates                                                                                                                                       |
| Tem_range         | numeric   | Temperature annual range  | Temperature annual range, calculated as the difference between the maximum temperature of the warmest month and minimum temperature of the coldest month. It was extracted from the CHELSA climate database using plot coordinates |
| Prec_Varieb       | numeric   | Precipitation seasonality | Precipitation seasonality, calculated as the coefficient of variation in precipitation across year, extracted from the CHELSA climate database using plot coordinates                                                              |
| grazing_intencity | numeric   | Grazing intensity         | Visual estimation of grazing intensity, ranging from 0 - no grazing to 3 - high grazing intensity                                                                                                                                  |
| mowing            | character | Presence of mowing        | Indication if there were recent evidences of mowing                                                                                                                                                                                |
| cover_litter      | numeric   | Cover litter (%)          | Percentage cover of litter after virtually removing all vegetation                                                                                                                                                                 |
| pH                | numeric   | soil pH                   | Soil pH level                                                                                                                                                                                                                      |
| Corg_percent      | numeric   | C total (%)               | Soil content of organic carbon, %                                                                                                                                                                                                  |

### `alpha_beta_gamma_community_variabl.csv`

Abbreviations used:

SR - species richness
ENSPIE - evenness measure calculated as inverse Simpson using species cover
cover - cumulative cover of plant community

Combines all diversity measures and plant cover for each scale:

- alpha diversity and cover measures (SR, ENSPIE, and cover) include doubled 10 m2 plots
- gamma diversity and cover measures (SR, ENSPIE, and cover) include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
- beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

#### Variables

| Short name | type      | Long name                                           | Description                                                                                                                                                                                                                                                                      |
| ---------- | --------- | --------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| dataset    | numeric   | Dataset ID                                          | Dataset ID, sampled by different teams and years                                                                                                                                                                                                                                 |
| series     | character | Series ID                                           | 100 m2 plot that includes two 10 m2 plots that are nested within it                                                                                                                                                                                                              |
| subplot    | character | Subplot                                             | one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner                                                                                                                                             |
| plotID     | character | PlotID                                              | Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots                                                                                                                                                                    |
| scale      | integer   | Spatial scale                                       | Grain size of the sampled plots: 10 is 10 m2 plots; 100 is 100 m2 plots, NA - nonaplicable (for beta diversity)                                                                                                                                                                  |
| type       | character | Scale type as used in the paper                     | Name of spatial scale: alpha - diversity and cover measures include doubled 10 m2 plots; gamma - diversity and cover measures include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity); beta diversity measures are calculated as gamma/alpha |
| metric     | character | Measure of plant community                          | Measure of biodiversity and cover of plant community: SR - species richness; ENSPIE - evenness measure calculated as inverse Simpson using species cover; cover - cumulative cover of plant community                                                                            |
| value      | numeric   | Value for the respective measure of plant community | Value for the respective measure of plant community |

### `climate_PCA.csv`

Contains scores for the compound climate variable, derived from the [PCA analysis](../analysis/01_PCA_climate.R),
see `../analysis/Readme.md` for details.

### `aggregation.csv`

Contains beta.BRAY.BAL - the proxy of intraspecific aggregation  for each plot. Spatial intraspecific aggregation was estimated by comparing dissimilarity in species covers between the two corners (i.e., two 10-m2 plots) within each 100-m2 plot. For this, we calculated  beta.BRAY.BAL - the balanced variation component of Bray–Curtis dissimilarity in species cover using ‘betapart’ package in R (Baselga & Orme, 2012). This measure is independent of total community abundance (total plant cover in our study) and measures the balanced variation in species abundance between two quadrats, i.e. when  cover increases for some species and decreases for others, maintaining similar total cover across quadrats, including also species turnover, where abundance of one species is replaced by other species (Baselga, 2017). Higher dissimilarity in covers of taxa between the two 10-m2 corners within the same 100-m2 plot implies higher intraspecific aggregation.

Baselga, A. (2017). Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution, 8(7), 799–808. https://doi.org/10.1111/2041-210X.12693

Baselga, A., & Orme, C. D. L. (2012). Betapart: An R package for the study of beta diversity. Methods in Ecology and Evolution, 3(5), 808–812. https://doi.org/10.1111/j.2041-210X.2012.00224.x

#### Variables

| Short name    | type      | Long name   | Description                                                         |
| ------------- | --------- | ----------- | ------------------------------------------------------------------- |
| series        | character | Series ID   | 100 m2 plot that includes two 10 m2 plots that are nested within it |
| beta.BRAY.BAL | numeric   | aggregation | Proxy of intraspecific aggregation of plant community               |
