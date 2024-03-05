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

## 1. setting for loops for all the species
# -------------------------------------
#load all the species we want to focus on: Saltmarsh sparrow, Seaside sparrow, Clapper Rail, Willet, Nelson's sparrow and potential saltmarsh sparrow hybrids
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
for (j in 1:length(speciesnames)){

## 1. Format Nest Observations
# -------------------------------------
# load cleaned focal species nest observation shapefile
nests<-st_read(paste0(path_out,"Intermediate_outputs/Nests/",speciesnames[j],"_nests_background_selected.shp"))%>%
    mutate(bp="p")%>%
    dplyr::select(-Region)
n
## 2. Add random veg and background points to nest data
# --------------------------------------------------------------------------
#load veg and background points

#load background point coordinates
bg_all<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m.csv"))

#make background points an sf shapefile
bg_points<-st_as_sf(bg_all,coords=c("x","y"),crs="EPSG:26918")%>%
  st_transform("EPSG:4269")%>%
  mutate(longitude=sf::st_coordinates(.)[,1],
         latitude=sf::st_coordinates(.)[,2])%>%
  dplyr::select(-region)


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
st_write(nests,paste0(path_out,"Intermediate_outputs/Nests/",speciesnames[j],"_nests_nonbuffed.shp"),delete_layer = T)

#nest and background buffers sf object
nests_buff<- nests%>%
  #buffer nest and background points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)
st_write(nests_buff,paste0(path_out,"Intermediate_outputs/Nests/",speciesnames[j],"_nests_buffed.shp"),delete_layer = T)


}

