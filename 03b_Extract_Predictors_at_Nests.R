library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(tmap)
library(exactextractr)

dat_path<-"D:/Data/"
path_out<-"D:/Outputs/"

if(!exists("nests")){
  source("03a_Format_Data_Nests_and_Predictors.R")
}


## Summarize at nest points
# ---------------------------------------------

## a. Marsh Vegetation Classes
#-------------------------------

#Create an empty list to hold a summary of vegetation in nests for each regional layer (will hold a list of summary tables, one for each region)
out_list<-list()
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(vg_cls[[1]]))
#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(vg_cls[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
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
veg_prop<-do.call("rbind",out_list)



## b. mean UVVR ratio at each nest
#-------------------------------
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(uvvr))
uvvr_mean<-terra::extract(uvvr, vect(nests), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_mean)


## c. Change in UVVR at each nest
#-------------------------------
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(uvvr_diff))
uvvr_diff2<-terra::extract(uvvr_diff, vect(nests), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_diff)#=uvvr_mean)#For test



## d. NDVI 
#-------------------------------
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(ndvi[[1]]))
for(i in 1:length(ndvi)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(ndvi[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
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
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(pca[[1]][[1]]))
for(i in 1:length(pca)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(pca[[i]][[1]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,pca=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
pca2<-do.call("rbind",out_list)



## f. Homogeneity (of NDVI)
#-------------------------------
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(txt_homo[[1]]))
for(i in 1:length(txt_homo)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_homo[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,hom_txt=value)
}

#empty region dataframes indicate no SALS nests in that region

#combine NDVI values for nests across all regions into 1 dataframe
txt_homo2<-do.call("rbind",out_list)


## g. Entropy (of NDVI)
#-------------------------------
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(txt_entro[[1]]))
for(i in 1:length(txt_entro)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_entro[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
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
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(txt_corr[[1]]))
for(i in 1:length(txt_corr)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_corr[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
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
precip<-rast(paste0(dat_path,"Precip/PRISM_ppt_30yr_normal_800mM4_annual_bil.bil"))
#set coord system of nest points to raster layer
nests<-st_transform(nests,crs(precip))

precip2<-terra::extract(precip, vect(nests), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,precip=PRISM_ppt_30yr_normal_800mM4_annual_bil)

### join NDVI, PCA, Texture, UVVR, and Proportion of vegetation classes at each nest into 1 table *Fix homogeneity and add precipitation
#---------------------------------------------------------------------------------------
final_dat_local<-left_join(nests,ndvi2, by='id')%>%
  left_join(pca2,by='id')%>%
  left_join(uvvr_diff2,by='id')%>%
  left_join(uvvr_mean,by='id')%>%
  left_join(txt_entro2,by='id')%>%
  left_join(txt_corr2,by='id')%>%
  left_join(dplyr::select(veg_prop,-latitude,-longitude),by='id')%>%
  left_join(precip2,by='id')%>%
  #dplyr::select(id,ndvi,hom_txt,ent_txt,cor_txt)%>%
  distinct(.keep_all = T)
final_dat_local$cor_txt[is.na(final_dat_local$cor_txt)]<-0

write.csv(final_dat_local,paste0(path_out,"SALS_nest_vars_local.csv"),row.names = F)


