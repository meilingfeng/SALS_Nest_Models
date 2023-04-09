
library(tidyverse)
library(sf)
library(gstat)
library(lme4)

## Set file path to data
# -------------------------------------------
dat_path<-"G:/My Drive/Research/SHARP/Data/"


## Load in Data
# -------------------------------------------

#Load nest fates with environmental covariates
fates<-read.csv(paste0(dat_path,"nest_fate_model_dataset.csv"))%>%
  #remove unimportant variables
  dplyr::select(-c("TERRBRD","MISSING"))
#convert missing data to 0, these are land cover classes that were not found within nest boundaries.
fates[is.na(fates)]<-0

#Load nest point locations
nests<-st_read(paste0(dat_path,"Nest Locations/nest_locations_11_10_22.shp"))%>%
  filter(id%in%fates$id)
coords<-nests%>%
  st_transform("EPSG:4269")%>%
  sf::st_coordinates()%>%
  cbind(nests)

#Add lat long to fates and covariates
fates<-left_join(coords[,c("id","X","Y")],fates,by="id")


#variogram of residuals-no apparent spatial auto corr
mod<-glm(fate~uvvr_mean+HIMARSH+LOMARSH+Y,data = fates,family = "binomial")
summary(mod)

fates$res<-mod$residuals  
dat<-dplyr::select(fates,Y,X,res,site)
gram<-variogram(res~1,data=dat[dat$site=="ER",])         
plot(gram)
fit.gram<-fit.variogram(gram, model = vgm(psill=2, model="Sph", range = 30, nuggest=0.5))
plot(gram,fit.gram)



#Plot predictors against response, look at shape of relationships
fates_long<-pivot_longer(fates,c("Y","HIMARSH","LOMARSH","PHRG","uvvr_mean","UPLND","STRM","POOL"),"var","value")
ggplot(fates_long, aes(x = value, y = fate)) +
  geom_point()+
  geom_smooth(method="gam", method.args = list( family = "binomial"))+
  facet_wrap(~var,scales = "free")
## which suggest a strong non-linearlity, although an odd one                    
ggplot(fates, aes(x = sqrt(uvvr_mean), y = fate)) +
  geom_smooth(method="lm", method.args = list( family = "binomial")) +
  theme_bw()
## I could be persuaded to use a sqrt term, but what does it mean?
#Maybe use a quadratic term for HIMARSH



#check for correlation between predictors
# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Plotting the correlation matrix
pairs(st_drop_geometry(fates[,c("HIMARSH","LOMARSH","uvvr_mean","Y")]),
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines

install.packages("PerformanceAnalytics")

library(PerformanceAnalytics)

chart.Correlation(st_drop_geometry(fates[,c("HIMARSH","LOMARSH","uvvr_mean","Y")]), histogram = TRUE, method = "pearson")

#check for independence assuptions
par(mfrow=c(2,2))
plot(mod)
par(mfrow=c(1,1))

