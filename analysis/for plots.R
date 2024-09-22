# Purpose: prepare data for plots
# Run first 
## 02_GLMM_alpha_10_div, 
## 03_GLMM_gamma_100_div,
## 04_GLMM_beta_div

# Re-scale predictions 

clima_pred_10m_rsc <- clima_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_10m_rsc <- Humus_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_10m_rsc <- Litter_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_10m_rsc <- pH_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_10m_rsc <- grazing_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_10m_rsc <- precipCV_pred_10m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

clima_pred_10m_Ensp_rsc <-clima_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_10m_Ensp_rsc <-Humus_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_10m_Ensp_rsc <- Litter_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_10m_Ensp_rsc <- pH_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_10m_Ensp_rsc <- grazing_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_10m_Ensp_rsc <- precipCV_pred_10m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))


clima_pred_100m_rsc <- clima_pred_100m  %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_100m_rsc <- Humus_pred_100m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_100m_rsc <- Litter_pred_100m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_100m_rsc <- pH_pred_100m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_100m_rsc <- grazing_pred_100m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_100m_rsc <- precipCV_pred_100m %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))


clima_pred_100m_Ensp_rsc <- clima_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_100m_Ensp_rsc <- Humus_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_100m_Ensp_rsc <- Litter_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_100m_Ensp_rsc <- pH_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_100m_Ensp_rsc <- grazing_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_100m_Ensp_rsc <- precipCV_pred_100m_Ensp %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))


clima_pred_beta_rsc <- clima_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_beta_rsc <- Humus_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_beta_rsc <- Litter_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_beta_rsc <- pH_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_beta_rsc <- grazing_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_beta_rsc <- precipCV_pred_beta%>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

clima_pred_beta_ENSPIE_rsc <- clima_pred_beta_ENSPIE %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Humus_pred_beta_ENSPIE_rsc <- Humus_pred_beta_ENSPIE %>%   
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

Litter_pred_beta_ENSPIE_rsc <- Litter_pred_beta_ENSPIE %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

pH_pred_beta_ENSPIE_rsc <- pH_pred_beta_ENSPIE %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

grazing_pred_beta_ENSPIE_rsc <- grazing_pred_beta_ENSPIE %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))

precipCV_pred_beta_ENSPIE_rsc <- precipCV_pred_beta_ENSPIE %>% 
  mutate(x1=x/max(x), predicted1=predicted/max(predicted))


# Plot standardized effects for SR -----


Fig_SR_clima  <- ggplot(clima_pred_10m_rsc, aes(x1, predicted1)) +
  labs(y="Species richness", x='Climate gradient')+
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=clima_pred_100m_rsc, linetype=1, size=1, col="#D6604D") +
  geom_line(data=clima_pred_beta_rsc, linetype=5, size=0.5, col="#00966F")


Fig_SR_Humus  <- ggplot(Humus_pred_10m_rsc, aes(x1, predicted1)) +
  labs(y="Species richness", x='Soil C')+
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=Humus_pred_100m_rsc, linetype=1, size=1, col="#D6604D") +
  geom_line(data=Humus_pred_beta_rsc, linetype=5, size=0.5, col="#00966F")

Fig_SR_Litter  <- ggplot(Litter_pred_10m_rsc, aes(x1, predicted1)) +
  labs(y="Species richness", x='Litter cover')+
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=Litter_pred_100m_rsc, linetype=1, size=1, col="#D6604D") +
  geom_line(data=Litter_pred_beta_rsc, linetype=1, size=1, col="#00966F")


Fig_SR_pH <- ggplot(pH_pred_10m_rsc, aes(x1, predicted1)) +
  labs(y="Species richness", x='Soil pH')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") +
  geom_line(data=pH_pred_100m_rsc, linetype=1, size=0.5, col="#D6604D") +
  geom_line(data=pH_pred_beta_rsc, linetype=1, size=1, col="#00966F")


Fig_SR_grazing <- ggplot(grazing_pred_10m_rsc, aes(x1, predicted1)) +
  scale_fill_manual(values = col)+  scale_color_manual(values = col) +
  labs(y="Species richness", x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") +
  geom_line(data=grazing_pred_100m_rsc, linetype=5, size=0.5, col="#D6604D") +
  geom_line(data=grazing_pred_beta_rsc, linetype=5, size=0.5, col="#00966F")


Fig_SR_precipCV <- ggplot(precipCV_pred_10m_rsc, aes(x1, predicted1)) +
  labs(y="Species richness", x='Precipitation CV')+
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=precipCV_pred_100m_rsc, linetype=1, size=1, col="#D6604D")+
  geom_line(data=precipCV_pred_beta_rsc, linetype=5, size=0.5, col="#00966F") 


# Plot standardized effects for ENSPIE -----

Fig_ENSPIE_Clima  <- ggplot(clima_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Climate gradient')+ 
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=clima_pred_100m_Ensp_rsc, linetype=1, size=1, col="#D6604D")+
  geom_line(data=clima_pred_beta_ENSPIE_rsc, linetype=1, size=0.5, col="#00966F") 


Fig_ENSPIE_Humus <- ggplot(Humus_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil C')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") +
  geom_line(data=Humus_pred_100m_Ensp_rsc, linetype=5, size=0.5, col="#D6604D")+
  geom_line(data=Humus_pred_beta_ENSPIE_rsc, linetype=1, size=1, col="#00966F") 


Fig_ENSPIE_Litter  <- ggplot(Litter_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Litter cover')+
  geom_line(linetype=1, size=0.5, col="#50A0C8") +
  geom_line(data=Litter_pred_100m_Ensp_rsc, linetype=1, size=1, col="#D6604D")+
  geom_line(data=Litter_pred_beta_ENSPIE_rsc, linetype=1, size=1, col="#00966F") 


Fig_ENSPIE_pH <- ggplot(pH_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Soil pH')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") +
  geom_line(data=pH_pred_100m_Ensp_rsc, linetype=5, size=0.5, col="#D6604D")+
  geom_line(data=pH_pred_beta_ENSPIE_rsc, linetype=5, size=0.5, col="#00966F") 



Fig_ENSPIE_grazing  <- ggplot(grazing_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Grazing intencity')+
  geom_line(linetype=5, size=0.5, col="#50A0C8") +
  geom_line(data=grazing_pred_100m_Ensp_rsc, linetype=1, size=0.5, col="#D6604D")+
  geom_line(data=grazing_pred_beta_ENSPIE_rsc, linetype=5, size=0.5, col="#00966F")  


Fig_ENSPIE_precipCV <- ggplot(precipCV_pred_10m_Ensp_rsc, aes(x1, predicted1)) +
  labs(y=expression(paste("ENS"[PIE])), x='Precipitation CV')+
  geom_line(linetype=1, size=1, col="#50A0C8") +
  geom_line(data=precipCV_pred_100m_Ensp_rsc, linetype=1, size=1, col="#D6604D") +
  geom_line(data=precipCV_pred_beta_ENSPIE_rsc, linetype=5, size=0.5, col="#00966F") 




# Merge plots -----
set_theme(base = theme_bw(), axis.textsize.x = 0.6, axis.textsize.y = 0.6, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 0.65,legend.pos = "None", 
          geom.linetype = 2)


Fig3 <- 
  Fig_SR_clima +
  Fig_SR_precipCV + 
  Fig_SR_Humus +
  Fig_SR_pH +
  Fig_SR_Litter + 
  Fig_SR_grazing +



  #
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol=1) & ylab(NULL) & theme(plot.margin = margin(3, 1, 3, 20), 
                                           plot.tag = element_text(size = 6, face = 'bold'), 
                                           plot.tag.position = c(0.15, 1.06))

Fig3


