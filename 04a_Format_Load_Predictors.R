library(sf)
library(tidyverse)
library(tidyterra)
#library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package


#################################################################################################
# Combine nest data and background data into one observation dataset
# Create a buffered version of the observation dataset
# load the environmental predictors at their original resolutions to sample at the observation points
#################################################################################################



## 3. Load buffered environmental data at their original resolutions (generated in script 02a)
# --------------------------------------------------------------------------


## Environmental Predictor 1: Marsh Vegetation Classes

#Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
zones<-c(1:8)

#list raster files
#file_list<-unlist(map(paste0(path_out,"Intermediate_outputs/HIMARSH"),~list.files(.,pattern = "vegcls_buff.tif$",full.names=T))) #maybe only need buffer for layers not derived from NAIP
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))


#read as raster layers
vg_cls<-map(file_list,rast)

#vegetation class definitions
veg_class<-data.frame(value=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))




## Environmental Predictor 2: Unvegetated-Vegetated ratio
uvvr<-  rast(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier_buff.tif"))%>%
  dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything())  
uvvr_diff<-  rast(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier_buff.tif"))%>%
  dplyr::rename_with(function(x){x<-'uvvr_diff'},.cols = everything()) 



## Environmental Predictor 3: NDVI
file_list1<-unlist(map(paste0(path_out,"Intermediate_outputs/NDVI/"),~list.files(.,pattern = "zeroed_NDVI_buff.tif$",full.names=T)))
ndvi<-map(file_list1,rast)



## Environmental Predictor 4: Principal Component of raw NAIP reflection bands (RBG)
#file_list2<-unlist(map(paste0(path_out,"Intermediate_outputs/PCA/"),~list.files(.,pattern = "PC1_buff.tif$",full.names=T)))
file_list2<-unlist(map(paste0(path_out,"Intermediate_outputs/PCA/"),~list.files(.,pattern = "PC1.tif$",full.names=T)))
pca<-map(file_list2,rast)



## Environmental Predictor 5: Precipitation
precip<-rast(paste0(path_out,"Intermediate_outputs/Precip/precip_buff.tif"))



## Environmental Predictor 6: tidal restrictions
tideres<- rast(paste0(path_out,"Intermediate_outputs/Tidal_restriction/tideres_buff.tif"))%>%
  dplyr::rename_with(function(x){x<-'tideres'},.cols = everything()) 



## Environmental Predictor 7: Texture
txt_corr<-map(unlist(map(paste0(path_out,"Intermediate_outputs/Texture/"),~list.files(.,pattern = "cor.*[0-9]_buff.tif$",full.names=T))),rast)
txt_entro<-map(unlist(map(paste0(path_out,"Intermediate_outputs/Texture/"),~list.files(.,pattern = "ent.*[0-9]_buff.tif$",full.names=T))),rast)


## Environmental Predictor 8: Elevation
file_list3<-unlist(map(paste0(path_out,"Intermediate_outputs/Elevation/"),~list.files(.,pattern = "DEM_buff.tif$",full.names=T)))
dem<-map(file_list3,rast)

