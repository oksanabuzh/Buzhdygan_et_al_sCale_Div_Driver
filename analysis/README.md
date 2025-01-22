# Folder "analysis"

Contains scripts for PCA analysis, mixed models, and figures.

# Files:

`01_PCA_climate.R` - performs PCA analysis for creating the composite variable of climate gradient, consisting of increasing precipitation and decreasing temperature.

`02_GLMM_alpha_10_div.R` - performs GLMM analysis for the 10-m2 plots for diversity measures (SR - species richness, ENSPIE - evenness).

`03_GLMM_gamma_100_div.R` - performs GLMM analysis for the 100-m2 plots for diversity measures.

`04_GLMM_beta_div.R` - performs GLMM analysis for beta diversity.

`05_a_Stand_Effect_size_alpha.R` - obtains the standardized coefficients of the predictor effects on diversity measures at the 10-m2 plots.

`05_b_Stand_Effect_size_gamma.R` - obtains the standardized coefficients of the predictor effects on diversity measures at the 100-m2 plots.

`06_ResidSpatCorr.R` - tests residual spatial autocorrelation.

`07_aggregation.R` - performs LMM analysis for spatial aggregation as response variable (100-m2 plots).

`08_Descriptive_stats.R` - performs descriptive and summary statistics, creates plots for the figures Fig. 1 b,c; Fig. S1 a-b; Fig. S2 a-d; Fig. S3 a,b,c,d; Fig. S4 a,b; Fig. S7.

`Fig_2.R` - creates plots in Fig.2 for the results of mixed models testing the effects of environmental drivers on species richness and evenness at the 10-m2 and 100-m2 plots.

`Fig_3a_4.R` - creates plots on relative strengths (standardized effect size) of the effects of each environmental driver on local diversity measures (species richness and ENSPIE) for the 10-m2 and 100-m2 plots (Fig.3a); and for the scale-dependent effects of environmental drivers on species richness and evenness at small (10 m2) and larger (100 m2) fine-grain plots (Fig.4).

`Fig_3b_R2.R` - creates plots in Fig.3b on relative importance (partial R2) of the environmental drivers in governing β-diversity measures and local diversity in 10-m2 and 100-m2.

`Fig_5.R` - creates plots in Fig. 5 for the results of mixed models testing the relationships between β-diversity and proximate factors: total cover, evenness and intraspecific aggregation.

`Fig_S5_S6.R` - creates plots in Fig.S5 and Fig.S6 for the results of mixed models testing the effects of environmental drivers on species richness and evenness at the 10-m2 and 100-m2 plots.

`Fig_S8.R` - creates plots in Fig.S8 for the results of mixed models testing the effects of environmental drivers on plant total cover at the 10-m2 and 100-m2 plots.