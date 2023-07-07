library(tidyverse)
library(terra)
library(sf)


########
#sample nest precitions at random veg sites
#############

zones<-c(2:6)
out_list<-list()

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. Load nest prediction rasters
pres_list<-unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Placement"),~list.files(.,pattern = "pres_BRTpreds_30mb.tif$",full.names=T)))
surv_list<-unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Success"),~list.files(.,pattern = "surv_BRTpreds_30mb.tif$",full.names=T)))

#read as raster layers
pres<-map(pres_list,rast)
surv<-map(surv_list,rast)

## 3. Load vegetation points
veg<-st_read(paste0(path_out, "Final_outputs/Veg_locations/veg_locations_12_29_22.shp"))%>%
  rename(id=veg_id)


### Sample probabilities at veg points
#--------------------------------------------

## 1. Probability of nest presence
for(i in 1:length(pres)) {
 
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(pres[[i]], vect(veg), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(lyr1))%>%
    dplyr::select(id,presence=lyr1)
}

#empty region dataframes indicate no SALS nests in that region
#combine regional dataframes
pres2<-do.call("rbind",out_list)%>%
  distinct(id,.keep_all = T)


## 2. Probability of nest success
for(i in 1:length(surv)) {
    #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(surv[[i]], vect(veg), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(lyr1))%>%
    dplyr::select(id,success=lyr1)
}

#Combine regional dataframes
surv2<-do.call("rbind",out_list)%>%
  distinct(id,.keep_all = T)



### nest presence, nest survival, and vegetation data into one dataframe
#---------------------------------------------------------------------------------------
final_dat<-left_join(veg,pres2, by='id')%>%
  left_join(surv2,by='id')%>%
  #dplyr::select(id,ndvi,hom_txt,ent_txt,cor_txt)%>%
  distinct(id,.keep_all = T)%>%
  filter(crd_typ!=1)%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         Year=year(mdy(date)))

if(!file.exists(paste0(path_out,"Final_outputs/Nest_Predictions_at_Veg_Points.csv"))){
  write.csv(st_drop_geometry(final_dat),paste0(path_out,"Final_outputs/Nest_Predictions_at_Veg_Points.csv"),row.names = F)
}


