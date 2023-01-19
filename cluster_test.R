library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(tmap)
library(exactextractr)
library(GLCMTextures)
#run test, get full zone 1 NDVI into cluster as well as UVVR data, finish loading pca

# test data
#uvvr<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=60)%>%
#  dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything())
#values(uvvr)<-runif(ncell(uvvr),0,2)


ndvi<-list()
ndvi[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
values(ndvi[[1]])<-runif(ncell(ndvi[[1]]),0,1)
ndvi[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
values(ndvi[[2]])<-runif(ncell(ndvi[[2]]),0,1)

#pca<-list()
#pca[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
#values(pca[[1]])<-runif(ncell(pca[[1]]),0,1)
#pca[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
#values(pca[[2]])<-runif(ncell(pca[[2]]),0,1)

#t_veg<-list()
#t_veg[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
#values(t_veg[[1]])<-round(runif(ncell(t_veg[[1]]),1,9))
#t_veg[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
#values(t_veg[[2]])<-round(runif(ncell(t_veg[[2]]),1,9))
#vg_cls<-t_veg

## Set file path to data
# -------------------------------------------
dat_path<-"D:/Research/SHARP/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Research/SHARP/Outputs/"




## 1. Format Nest Observations
# -------------------------------------
# load nest observation shapefile
if(file.exists(paste0(path_out,"nest.rda"))){
  load(paste0(path_out,"nest.rda"))
}

if(!file.exists(paste0(path_out,"nest.rda"))){
nests<-st_read(paste0(dat_path,"Nest_Locations/nest_locations_01_3_23.shp"))%>%
  # filter records to just SALS or species of interest
  filter(Species=="SALS"&
           # also filter records missing coordinate information or that have coordinate errors
           crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  # convert all coordinates to degrees (NAD83) and create Long and Lat columns
  st_transform("EPSG:4269")%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         #create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                        fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(name=id,latitude=Lat,longitude=Long,fate,Year,site=site_cd)

# output table file
#write.csv(st_drop_geometry(nests),paste0(dat_out,"Nest_Coords_fates_SALS_01_5_23.csv"),row.names=F)
nests<-rename(nests,id=name)

nests_buff<- nests%>%
  #buffer nest points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)

}
## 2. Sample environmental data at nests using their original resolutions
# --------------------------------------------------------------------------

## 2. a. Load in data
#---------------------

if(!file.exists(paste0(path_out,"enviro.rda"))){

## Environmental Predictor 1: Marsh Vegetation Classes

#Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
zones<-c(1:8)

#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

#read as raster layers
vg_cls<-map(file_list,rast)

#set coordinate systems
# Set all coordinate systems to ESPG:4269 (NAD 1983 geographic coord system)
#for(i in 1:length(vg_cls)){
#vg_cls[[i]] <- terra::project(vg_cls[[i]],"EPSG:4269")
#}

#get cell area
area<-prod(res(vg_cls[[1]]))

#vegetation class definitions
veg_class<-data.frame(value=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))




## Enviromental Predictor 2: Unvegetated-Vegetated ratio
# transformed uvvr coordinate system from WGS to NAD83 in ArcPro
# just want band 4 (UVVR). Band 1 is unveg, 2 is veg, and 3 is proportion of water.

# A. UVVR mean across 2014-2018
uvvr<-rast(paste0(dat_path,"UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))%>%
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything()) 
# values above 2 are not accurate, set to NA
uvvr[uvvr>2]<-NA

# B. UVVR change (2014-2018)
uvvr14<-rast(paste0(dat_path,"UVVR/uvvr_14_utm18_mean_ext.tif"))%>%
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_14'},.cols = everything())
# values above 2 are not accurate, set to NA
uvvr14[uvvr14>2]<-NA
uvvr18<-rast(paste0(dat_path,"UVVR/uvvr_18_utm18_mean_ext.tif"))%>%
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_18'},.cols = everything())
uvvr18[uvvr18>2]<-NA
# subtract 2018 from 2014 to get change
uvvr_diff<-uvvr18-uvvr14
uvvr_diff<-uvvr_diff%>%
  dplyr::rename_with(function(x){x<-'uvvr_diff'},.cols = everything())



## Environmental Predictor 3: NAIP

#list raster files
file_list1<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))
file_list2<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "PCA.tif$",full.names=T)))

#read as raster layers
ndvi<-map(file_list1,rast)
pca<-map(file_list2,rast)

#set values below 0 to 0 - barren land or water
for(i in 1:length(ndvi)){
  dat<-ndvi[[i]]
  dat[dat<0]<-0
  ndvi[[i]]<-dat
}
}

## 3. Calculate texture surfaces
# --------------------------------------------
if(file.exists(paste0(path_out,"enviro.rda")) &
   !file.exitst(paste0(path_out,"txt.rda"))){
  
  load(paste0(path_out,"enviro.rda"))

  #save(nests,nests_buff,file = paste0(path_out,"nest.rda"))
  #save(dat_path,zones,area,path_out,file = paste0(path_out,"param.rda"))))

## 3. a. Raster Quanitization 
#----------------------------
# use equal range method when comparing across several datasets (can set gloabl min and max range using max_val and min_val) - splits data into equal ranges
# alternative option is "equal prob" which splits by quantiles, used in original paper (Haralick and Shanmugam 1973)
# NDVI of live plants ranges 0:1, set range and breaks at 0.1 intervals
txt_homo<-txt_entro<-txt_corr<-list()

for (i in 1:length(ndvi)){
  ndvi_rq<- quantize_raster(r = ndvi[[i]], n_levels = 10, method = "equal range")
  
  ## 3. b. calc texture metric surfaces
  #----------------------------
  # use window size of 3x3m 
  # calculate across multiple shifts (invariant)
  txt<- glcm_textures(ndvi_rq, w = c(3,3), n_levels =10, 
                      metrics = c("glcm_homogeneity", 
                                  "glcm_entropy", "glcm_correlation"),
                      quantization = "none", shift = list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)))
  
  txt_homo[[i]]<-txt[[1]]
  txt_entro[[i]]<-txt[[2]]
  txt_corr[[i]]<-txt[[3]]
  
  }

}


if(!file.exists(paste0(path_out,"1_hom_txt.tif"))){
  for(i in 1:length(txt_homo)){
writeRaster(txt_homo[[i]],file=paste0(path_out,i,"_hom_txt.tif"),overwrite=T)
writeRaster(txt_entro[[i]],file=paste0(path_out,i,"_ent_txt.tif"),overwrite=T)
writeRaster(txt_corr[[i]],file=paste0(path_out,i,"_cor_txt.tif"),overwrite=T)
  }
}
  


if(!file.exists(paste0(path_out,"uvvr_noOutlier.tif"))){
  writeRaster(uvvr,file = paste0(path_out,"uvvr_mean_noOutlier.tif"),overwrite=T)
  writeRaster(uvvr_diff,file = paste0(path_out,"uvvr_diff_noOutlier.tif"),overwrite=T)
}

if(!file.exists(paste0(path_out,"nest.rda"))){
save(nests,nests_buff,file = paste0(path_out,"nest.rda"))
}
