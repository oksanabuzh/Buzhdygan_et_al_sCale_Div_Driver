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

### Environm_variab.csv 

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

### alpha_beta_gamma_community_variabl.csv
Combines all diversity measures and plant cover for each scale:
- alpha diversity and cover measures (SR, ENSPIE, and cover) include doubled 10 m2 plots
- gamma diversity and cover measures (SR, ENSPIE, and cover) include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
- beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

SR - species richness

ENSPIE - evenness measure calculated as inverse Simpson using species cover

cover - cumulative cover of plant community




Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| dataset	  |numeric | Dataset ID	| Dataset ID, sampled by different teams and years |
|series	| character	| Series ID	| 100 m2 plot that includes two 10 m2 plots that are nested within it |
|subplot	| character	| Subplot 	| one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner |
|plotID	| character	| PlotID	v Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots|
|scale	| integer |	Spatial scale | Grain size of the sampled plots: 10 is 10 m2 plots; 100 is 100 m2 plots, NA - nonaplicable (for beta diversity)|
|type	| character	| Scale type as used in the paper	| Name of spatial scale: alpha - diversity and cover measures include doubled 10 m2 plots; gamma - diversity and cover measures include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity); beta diversity measures are calculated as gamma/alpha |
|metric	| character	| Measure of plant community	Measure of biodiversity and cover of plant community: SR - species richness; ENSPIE - evenness measure calculated as inverse Simpson using species cover; cover - cumulative cover of plant community |
|value	| numeric	| Value for the respective measure of plant community	Value for the respective measure of plant community|


### climate_PCA