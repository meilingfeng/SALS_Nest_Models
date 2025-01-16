library(sf)
library(tidyverse)
library(terra)#updated version of raster package
library(dismo)


##############################################################################################
#Start Species Specific Analysis Here. Summarize data availability after spatial data cleaning.
##############################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## 1. Format Nest Observations
# -------------------------------------

# load nest observation shapefile
nests_all<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%#already in UTM18 NAD83 "EPSG:26918"
  st_buffer(50)
st_write(nests_all,paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23_50buf.shp"))

nests_all2<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  st_transform(crs=4269)%>% #NAD83 for lat long
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])%>%
  st_drop_geometry()%>%
  select(id,Species,Year,Site, State, long,lat)
write.csv(nests_all2,paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.csv"))

