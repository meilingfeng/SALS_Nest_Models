

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)

## Plotting spatial data 
library(spData)
library(tmap)
library(viridis)

### Set up
# -------------------------------------------
source("05d_Model_Selection_Evaluation.R")

# spatial layers

# Raster predictor layers for a particular zone
predictors<-rast(unlist(file_list_all_zones[[8]]))
#name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","HIMARSH","LOMARSH")

# East coast 
data(us_states)
states<-us_states%>%filter(REGION=="Norteast")%>%st_transform(crs(predictors[[1]]))

# Nest locations
nests<-st_as_sf(pres_dat, coords=c("longitude","latitude"),crs=crs(predictors[[1]]))%>%
  mutate(fate = as.factor(ifelse(is.na(fate), 0, fate)))


tmaptools::palette_explorer()

## Map of observations
#--------------------------------------------------
base<-tm_shape(states) + tm_borders()  
  base+
    tm_shape(filter(nests,bp=="b")) + tm_dots(col="black")+
  tm_shape(filter(nests,bp=="p")) + tm_dots(col="fate",palette=c("darkgreen","brown"))

# you can overlay raster maps
tm_shape(india) + tm_borders()  +
  tm_shape(map_india) + tm_raster()


## combine mulitple panes

tmap_arrange(map_world + tm_shape(north_america) + tm_fill(), 
             tm_shape(north_america) + tm_fill())

## say we want different colours for different countries
tm_shape(north_america) + tm_fill(col = "name_long")

# to explore the inbuilt colours, look at.
##tmaptools::palette_explorer()

## It is easy to add a compass and scale bar

na_map <- tm_shape(north_america) + 
  tm_fill(col = "name_long", title = "Country name") ## note legend title
na_map + tm_compass(position = c("right", "bottom"), type = "4star") +
  tm_scale_bar()
# We might also want to change the layout (equivalent to theme in ggplot)

na_map + tm_compass(position = c("right", "bottom"), type = "4star") +
  tm_scale_bar() + tm_layout(bg.color = "aliceblue")
## there are also built in themes
na_map + tm_style("classic")

#Map of predictors
#-------------------------------------------------




#variable contributions
#--------------------------------------------