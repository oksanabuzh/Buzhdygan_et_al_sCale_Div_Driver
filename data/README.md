# Folder "data"

## About

Contains data that are used for analysis 

#### Structure

| Document                               | What's inside                            |
| -------------------------------------- |----------------------------------------- |
| Environm_variab.csv                    | Environmental data                       |
| alpha_beta_gamma_community_variabl.csv | Plant diversity measures and plant cover for each scale |
| climate_PCA                            | scores from the [PCA analysis](analysis/PCA_climate.R)  |       


## Metadata

### -> Environm_variab.csv 
This dataset contains all measured environmental data, that were used in the paper

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| dataset	  |numeric | Dataset ID	| Dataset ID, sampled by different teams and years |
|series	| character	| Series ID	| 100 m2 plot that includes two 10 m2 plots that are nested within it |
|subplot	| character	| Subplot 	| one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner |
| plotID	| character	| PlotID	| Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots|
|habitat_broad	| 	character	| Grassland habitat type	| Broad types of grassland habitats, distinguished based on EUNIS habitat type system |
|zonality	| 	character		| Zonality of vegetation		| Zonality:  zonal (natural) vegetation types, azonal (seminatural) grassland types |
|BIO1		| numeric	| 	Mean annual temperature		| Mean annual temperature, extracted from the CHELSA climate database using plot coordinates |
|BIO12	| 	numeric		| Mean annual precipitation		| Mean annual precipitation, extracted from the CHELSA climate database using plot coordinates |
|BIO7	| 	numeric		| Temperature annual range		| Temperature annual range, calculated as the difference between the maximum temperature of the warmest month and minimum temperature of the coldest month. It was extracted from the CHELSA climate database using plot coordinates |
|BIO15	| 	numeric		| Precipitation seasonality	| 	Precipitation seasonality, calculated as the coefficient of variation in precipitation across year, extracted from the CHELSA climate database using plot coordinates |
|grazing_intencity	| 	numeric		| Grazing intensity 	| 	Visual estimation of grazing intensity, ranging from 0 - no grazing to 3 - high grazing intensity |
|mowing		| character		| Presence of mowing		| Indication if there were recent evidences of mowing |
|cover_litter	| 	numeric		| Cover litter (%) 	| Percentage cover of litter after virtually removing all vegetation |
|pH	| 	numeric		|  soil pH		|  Soil pH level |
|Corg_percent		| numeric		|  C total (%)		|  Soil content of organic carbon, % |



### -> alpha_beta_gamma_community_variabl.csv
Combines all diversity measures and plant cover for each scale:
- alpha diversity and cover measures (SR, ENSPIE, and cover) include doubled 10 m2 plots
- gamma diversity and cover measures (SR, ENSPIE, and cover) include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
- beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

SR - species richness

ENSPIE - evenness measure calculated as inverse Simpson using species cover

cover - cumulative cover of plant community

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| dataset	  |numeric | Dataset ID	| Dataset ID, sampled by different teams and years |
|series	| character	| Series ID	| 100 m2 plot that includes two 10 m2 plots that are nested within it |
|subplot	| character	| Subplot 	| one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner |
|plotID	| character	| PlotID	| Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots|
|scale	| integer |	Spatial scale | Grain size of the sampled plots: 10 is 10 m2 plots; 100 is 100 m2 plots, NA - nonaplicable (for beta diversity)|
|type	| character	| Scale type as used in the paper	| Name of spatial scale: alpha - diversity and cover measures include doubled 10 m2 plots; gamma - diversity and cover measures include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity); beta diversity measures are calculated as gamma/alpha |
|metric	| character	| Measure of plant community |	Measure of biodiversity and cover of plant community: SR - species richness; ENSPIE - evenness measure calculated as inverse Simpson using species cover; cover - cumulative cover of plant community |
|value	| numeric	| Value for the respective measure of plant community	| Value for the respective measure of plant community|


### -> climate_PCA
Contains scores for the compound climate variable, derived from the [PCA analysis](analysis/PCA_climate.R)