library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package


#################################################################################################
# Combine nest data and background data into one observation dataset
# Create a buffered version of the observation dataset
# load the environmental predictors at their original resolutions to sample at the observation points
#################################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"



## 1. Format Nest Observations
# -------------------------------------
# load cleaned focal species nest observation shapefile

nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_dist_err_removed.shp"))%>%
  mutate(bp="p")%>%
  dplyr::select(-lat_bin)

# output a csv file for nest coordinates and their fate, if recorded
if(!file.exists(paste0(path_out,"Final_outputs/Nest_Coords_fates_SALS_06_21_23.csv"))){
write.csv(st_drop_geometry(nests),paste0(path_out,"Final_outputs/Nest_Coords_fates_SALS_06_21_23.csv"),row.names=F)
}




## 2. Add random veg and background points to nest data
# --------------------------------------------------------------------------
#load veg and background points
source("03_background_selection.R")

# crop to nest extent
#bg_points2<-st_crop(bg_points,nests)%>%mutate(Region=NA)
bg_points2<-bg_points# dont need to for bg since we selected only from zones with nests
veg2<-st_crop(veg,nests)

#make sure columns align
names(nests)
names(veg2)
names(bg_points2)

#add points to nest data
nests<-rbind(nests,bg_points2,veg2)%>%
  distinct(.keep_all = TRUE)

#nest and background buffers sf object
nests_buff<- nests%>%
  #buffer nest and background points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)




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
  




