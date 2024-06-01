# Purpose Fig 3 a, Fig. 4

## Barplots of the standardized effects of diversity drivers  
### on alpha and gamma richnes and ENSPIE Fig 3 a
## Scatterplot of alpha vs gamma Fig. S4

library(tidyverse)
library(MetBrewer)
library(patchwork)

# data----
alpha_st.eff <-read.csv ("results/alpha_st.eff.csv")%>% rename(P_sign="X.1") %>% select(-X,  -P.Value, -P_sign) 
gamma_st.eff <-read.csv ("results/gamma_st.eff.csv") %>% rename(P_sign="X.1")%>% select(-X,  -P.Value, -P_sign) 


alpha <- alpha_st.eff %>% 
  rename(alpha.st.est="Std.Estimate") %>% 
  mutate(Measure=dplyr::recode(Response, 
                                   'alpha_10_div'='SR',
                                    'log(alpha_10_ENSPIE)'='ENSPIE')) %>% 
  rename(st.est=alpha.st.est)

gamma <- gamma_st.eff %>% 
              rename(gamma.st.est="Std.Estimate")  %>% 
              mutate(Measure=dplyr::recode(Response, 
                                   'gamma_100_div'='SR',
                                    'log(gamma_100_ENSPIE)'='ENSPIE')) %>% 
  rename(st.est=gamma.st.est)

alpha    
gamma

# Bind datasets
Difference <-   bind_rows(alpha=alpha, gamma=gamma, .id="scale") %>% 
  mutate(Driver= dplyr::recode(Predictor,  
                               clima_Comp="Climate gradient",
                               clima_Comp='Climate gradient',
                                   PrecipCV_Comp='Precipitation variation',
                                   Corg_Comp='Soil C',
                               Corg_percent ='Soil C',
                                   pH='Soil pH',
                                   pH_Comp='Soil pH',
                                   Litter_Comp='Litter cover',
                                   grazing_intencity='Grazing intencity',
                                   mowing='Mowing' )) %>% 
  mutate(St.est=abs(st.est))

Difference

Difference$Driver <- factor(Difference$Driver, 
                            levels = rev(c("Climate gradient",
                                           "Precipitation variation", 
                                           "Soil C", 
                                           "Soil pH", 
                                           "Litter cover", 
                                           "Grazing intencity",
                                           "Mowing"))
                            ) 

Difference$Measure <- factor(Difference$Measure, levels=c("SR", "ENSPIE"))
Difference$scale <- factor(Difference$scale, levels=c("gamma", "alpha"))

# Fig. 3 a -----

plot <- ggplot(Difference, aes(y = Driver, x = St.est, fill = scale)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  facet_wrap(~Measure, scales = "free_x") +
  # MetBrewer::scale_fill_met_d("Kandinsky") +
  scale_fill_manual(values = c("alpha"= "#A6CEE3",
                               "gamma"= "#D6604D"))+
  theme_bw()
plot
ggsave("plot.png", plot, width = 10, height = 5, units = "in", dpi = 300)



plot_SR <- ggplot(Difference %>% filter(Measure=="SR"), 
                  aes(y = Driver, x = St.est, fill = scale)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  labs(y=" ", x=" ")+
  xlim(0,0.42)+
  scale_fill_manual(values = c("alpha"= "#A6CEE3",
                               "gamma"= "#D6604D"))+
  theme_bw()+
  theme(legend.position="None",
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 13))
plot_SR



plot_ENSPIE <- ggplot(Difference %>% filter(Measure=="ENSPIE"), 
                  aes(y = Driver, x = St.est, fill = scale)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  labs(y=" ", x=" ")+
  xlim(0,0.6)+
  scale_fill_manual(values = c("alpha"= "#A6CEE3",
                               "gamma"= "#D6604D"))+
  theme_bw()+
  theme(legend.position="None",
        axis.text.x = element_text(size = 11), 
        axis.text.y =element_blank())
plot_ENSPIE

# Combine plots

## Plots combined ----
# for alpha SR (10 m^2)

plot_SR + plot_ENSPIE +
  plot_layout(ncol=2)



Difference <-   bind_rows(alpha=alpha, gamma=gamma, .id="scale") %>% 
  mutate(Driver= dplyr::recode(Predictor,  
                               clima_Comp="Climate gradient",
                               clima_Comp='Climate gradient',
                               PrecipCV_Comp='Precipitation variation',
                               Corg_Comp='Soil C',
                               Corg_percent ='Soil C',
                               pH='Soil pH',
                               pH_Comp='Soil pH',
                               Litter_Comp='Litter cover',
                               grazing_intencity='Grazing intencity',
                               mowing='Mowing' )) %>% 
  mutate(St.est=abs(st.est))

Mismatch <- Difference %>% 
  select(scale, Measure, Driver, St.est) %>% 
  pivot_wider(names_from = scale, values_from = St.est)


SR <- ggplot(Mismatch %>% filter(Measure=="SR"), 
       aes(x = alpha, y = gamma, 
           color = Driver)) + 
             # fill=dummy_fill, shape=AG_BG_stock)) +
  geom_abline(intercept = 0, slope = 1, col="gray", linewidth=1) +
    geom_point(size=3, stroke=1.1 ) +
#  scale_shape_manual(values=MyShape2) +
#  scale_color_manual(values=myPalette)+
  # scale_fill_manual(values=myPalette_with_white)+
  # facet_wrap( ~  in_out_flux, scales = "free") +
  theme_bw() +  # guides(fill = "none") +
  labs(x="Effect size at alpha scale", y= "Effect size at gamma scale", title="Species richness",
      # shape="AG/BG subnetwork", 
      col="Driver")


ENSPIE <- ggplot(Mismatch %>% filter(Measure=="ENSPIE"), 
             aes(x = alpha, y = gamma, 
                 color = Driver)) + 
  # fill=dummy_fill, shape=AG_BG_stock)) +
  geom_abline(intercept = 0, slope = 1, col="gray" , linewidth=1) +
  geom_point(size=3, stroke=1.1 ) +
  #  scale_shape_manual(values=MyShape2) +
  #  scale_color_manual(values=myPalette)+
  # scale_fill_manual(values=myPalette_with_white)+
  # facet_wrap( ~  in_out_flux, scales = "free") +
  theme_bw() +  # guides(fill = "none") +
  labs(x="Effect size at alpha scale", y= "Effect size at gamma scale", title=expression(ENS["PIE"]),
       # shape="AG/BG subnetwork", 
       col="Driver")

SR + ENSPIE +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides = "collect", ncol=2) & # theme(legend.position = 'right') +
  # plot_layout(ncol=2) & # ylab(NULL) & 
  theme(legend.position = 'right',
        plot.margin = margin(10, 30, 6, 30), 
        plot.tag = element_text(size = 12, face = 'bold'), 
        plot.tag.position = c(0.08, 1.06))



# Fig. 4 -----
ggplot(Mismatch %>% 
         mutate(Measure= dplyr::recode(Measure,  
                                      SR="Species richness",
                                      ENSPIE="ENSPIE")), 
                aes(x = alpha, y = gamma, 
                     color = Driver, shape=Measure)) + 
  # fill=dummy_fill, shape=AG_BG_stock)) +
  geom_abline(intercept = 0, slope = 1, col="gray" , linewidth=1) +
  geom_point(size=4, stroke=1.1 ) +
  scale_shape_manual(values=c(19,21)) +
  #  scale_color_manual(values=myPalette)+
  # scale_fill_manual(values=myPalette_with_white)+
  # facet_wrap( ~  in_out_flux, scales = "free") +
  theme_bw() +  # guides(fill = "none") +
  labs(x="Effect size at alpha scale", y= "Effect size at gamma scale",
       shape="Diversity measure", 
       col="Driver")

