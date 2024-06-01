library(tidyverse)


k.dat<-read.csv ("data/headers.csv", header=T)

str(k.dat) 
names (k.dat)


PCA <- read_csv("data/PCA2.csv")

k.dat <- k.dat %>% 
  full_join(PCA, by = "plotID") 


corel <- k.dat %>% 
  # select_if(is.numeric)
   select(
     built_up_2km, build_up_5km,
     heat_index,
     cover_litter,
     burning, grazing, mowing, ploughing, abandonment, encroachment, grazing_intencity, 
          economic_use, 
          cover_shrub_total)




summary(corel)

str(corel)
names(corel)


corl1 <- round(cor(corel, use="pairwise.complete.obs", method = c("pearson")),2)
corl1

write.csv(corl1, "data/correlation.csv",  row.names = T)

library(corrplot)
x11(height=9,width=8.5)
corrplot(corl1, type = "upper",  
         tl.col = "black", tl.srt = 50)


library(ggcorrplot)

x11(height=9,width=8.5)
ggcorrplot(corl1, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, 
           colors = c("red", "white", "blue"))



# Correlations among predictors (own data)

varb <- k.dat%>%
  select(altitude, 
         # heat_index, 
         microrelief,
        # depth_mean,
         cover_shrubs, cover_litter, cover_stones,
         height_herb_mean,
         grazing, mowing, ploughing, abandonment,
        # sand, silt, clay,
        pH,
       # CaCO3,
       C_percent, # C total%
       Corg,  # C organic
       cover_total_sum
       #, eunis  # habitat type
       )


str(varb) 
varb$heat_index <-  as.numeric(varb$heat_index)

corl2 <- round(cor(na.omit(varb), method = c("pearson")),1)
corl2
write.csv(round(corl2,1), "correl2.csv", row.names = T)


library(ggcorrplot)

x11(height=9,width=8.5)
ggcorrplot(corl2, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, 
           colors = c("red", "white", "blue"))

