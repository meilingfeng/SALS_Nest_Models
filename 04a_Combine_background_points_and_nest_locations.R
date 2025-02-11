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
path_out<-"D:/Nest_Models/Outputs/"

## setting for loops for all the species
# -------------------------------------
#load all the species we want to focus on: Saltmarsh sparrow, Seaside sparrow, Clapper Rail, Willet, Nelson's sparrow
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")

for (j in 1:length(speciesnames)){

## 1. Format Nest Observations
# -------------------------------------
# load cleaned focal species nest observation shapefile
nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))%>%
    mutate(bp="p")%>%
  select(id,site=site_cd, bp, Year, Region, fate, Long,Lat)

    
## 2. Add random veg and background points to nest data
# --------------------------------------------------------------------------
#load background point coordinates
bg_points<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"))

#make background points an sf vector
bg_points<-st_as_sf(bg_points,coords=c("x","y"),crs="EPSG:26918")%>%
  mutate(Long=sf::st_coordinates(.)[,1],
         Lat=sf::st_coordinates(.)[,2])


#make sure columns align
names(nests)
names(bg_points)


#add bg points to nest data
nests<-rbind(nests,bg_points)%>%
  distinct(.keep_all = TRUE)


#and write the combined point data to file
st_write(nests,paste0(path_out,"Intermediate_outputs/Nest_locations/",speciesnames[j],"_nest_pres_bg_points.shp"),delete_layer = T)



# Also write as a buffered sf object
nests_buff<- nests%>%
  #buffer nest and background points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)
st_write(nests_buff,paste0(path_out,"Intermediate_outputs/Nest_locations/",speciesnames[j],"_nest_pres_bg_buff.shp"),delete_layer = T)
}

