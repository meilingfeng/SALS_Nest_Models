library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(exactextractr)


#############################################################
## Summarize environmental predictors at nest points
#############################################################



## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


# Load all nest observations and environmental predictor data
if(!exists("veg_class")){
  source("04a_Format_Load_Predictors.R")
}

#load all the species we want to focus on: Saltmarsh sparrow, Seaside sparrow, Clapper Rail, Willet, Nelson's sparrow and potential saltmarsh sparrow hybrids
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
for (j in 1:length(speciesnames)){
  nests<-st_read(paste0(path_out,"Intermediate_outputs/Nests/",speciesnames[j],"_nests_nonbuffed.shp"),stringsAsFactors = FALSE)

## a. Marsh Vegetation Classes
#-------------------------------

#Create an empty list to hold a summary of vegetation in nests for each regional layer (will hold a list of summary tables, one for each region)
out_list<-list()
#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(vg_cls[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(vg_cls[[1]]))), bind = T)%>%
    mutate(region=zones[i])%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    left_join(veg_class,by="value")%>%
    filter(!is.na(value))%>%
    dplyr::select(-value,-Year,-site,-fate)
}

#empty region dataframes indicate no SALS nests in that region

#combine vegetation values for nests across all regions into 1 dataframe
veg_prop<-do.call("rbind",out_list)%>%
  distinct(id,latitude,longitude,.keep_all = T)


## b. mean UVVR ratio at each nest
#-------------------------------
uvvr_mean<-terra::extract(uvvr, vect(st_transform(nests,crs(uvvr))), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_mean)


## c. Change in UVVR at each nest
#-------------------------------
uvvr_diff2<-terra::extract(uvvr_diff, vect(st_transform(nests,crs(uvvr_diff))), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_diff)



## d. NDVI 
#-------------------------------
for(i in 1:length(ndvi)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(ndvi[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(ndvi[[1]]))), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,ndvi=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
ndvi2<-do.call("rbind",out_list)



## e. PCA
#-------------------------------
for(i in 1:length(pca)) {
  #rename raster values
  layer<-dplyr::rename_with(pca[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(pca[[1]]))), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,pca=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
pca2<-do.call("rbind",out_list)




## g. Entropy (of NDVI)
#-------------------------------
for(i in 1:length(txt_entro)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_entro[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(txt_entro[[1]]))), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,ent_txt=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
txt_entro2<-do.call("rbind",out_list)




## h. Correlation (of NDVI)
#-------------------------------
for(i in 1:length(txt_corr)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_corr[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(txt_corr[[1]]))), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,cor_txt=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
txt_corr2<-do.call("rbind",out_list)



## i. Precipitation
#-------------------------------
precip2<-terra::extract(precip, vect(st_transform(nests,crs(precip))), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,precip=focal_mean)



## j. Tidal Restrictions
#--------------------------------
tideres2<-terra::extract(tideres, vect(st_transform(nests,crs(tideres))), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,tideres)


## k. Elevation
#-------------------------------
for(i in 1:length(dem)) {
  #rename raster values
  layer<-dplyr::rename_with(dem[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(st_transform(nests,crs(dem[[1]]))), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,elevation=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
dem2<-do.call("rbind",out_list)

### join NDVI, PCA, Texture, UVVR, and vegetation classes at each nest into 1 table
#---------------------------------------------------------------------------------------
final_dat_local<-left_join(nests,ndvi2, by='id')%>%
  left_join(pca2,by='id')%>%
  left_join(uvvr_diff2,by='id')%>%
  left_join(uvvr_mean,by='id')%>%
  left_join(txt_entro2,by='id')%>%
  left_join(txt_corr2,by='id')%>%
  left_join(dplyr::select(veg_prop,-latitude,-longitude,-bp),by='id')%>%
  left_join(precip2,by='id')%>%
  left_join(tideres2,by='id')%>%
  left_join(dem2,by='id')%>%
  #dplyr::select(id,ndvi,hom_txt,ent_txt,cor_txt)%>%
  distinct(id,.keep_all = T)
final_dat_local$cor_txt[is.na(final_dat_local$cor_txt)]<-0


write.csv(st_drop_geometry(final_dat_local),paste0(path_out,"Final_outputs/",speciesnames[j],"_nest_vars_local.csv"),row.names = F)
}
