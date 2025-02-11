
library(tidyverse)
library(terra)
library(sf)

# get the latitudinal extent of each focal species breeding range within the Maine to Virginia (Delmarva Peninsula) study region

# use these extents to limit background selection and nesting predictions.


#file path names
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#Focal species list
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")

range_limits<-data.frame(species=speciesnames,ymin=rep(NA,length(speciesnames)),ymax=rep(NA,length(speciesnames)))

# SALS has a northern limit in the study region 
# get y max from nest records, exceeds the current range limit 
# source (Ruskin et al 2023) doi/10.1002/ece3.10532
nests_SALS<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_valid_nest_locations_2010_2024.shp"))%>%
  st_transform(crs="EPSG:26918")

range_limits[range_limits$species=="SALS","ymax"]<-ext(nests_SALS)[4]


# NESP has a southern breeding limit in the study range (they have the highest north breeding range)
# get y min from most southerly record of a breeding male in sandy neck beach park in Barnstable MA
# source (Ruskin et al 2023) doi/10.1002/ece3.10532
range_NESP<-st_read(paste0(dat_path,"Breeding_Ranges/NESP_southern_breeding_range_limit.shp"))%>%
  st_transform(crs="EPSG:26918")
range_limits[range_limits$species=="NESP","ymin"]<-ext(range_NESP)[3]


# SESP has a northern breeding limit within the study region 
# get latitudinal max from range map 
# source: U.S. Geological Survey (USGS) - Gap Analysis Project (GAP), 2018, Seaside Sparrow (Ammodramus maritimus) bSESPx_CONUS_2001v1 Range Map: U.S. Geological Survey data release, https://doi.org/10.5066/F7988614.
range_SESP<-st_read(paste0(dat_path,"Breeding_Ranges/bSESPx_CONUS_Range_2001v1.shp"))%>%
  st_transform(crs="EPSG:26918")
range_limits[range_limits$species=="SESP","ymax"]<-ext(range_SESP)[4]


# WILL completely overlaps the study region, don't need to limit


# CLRA has northern breeding limit in the study region
# base on breeding range 
# source: U.S. Geological Survey (USGS) - Gap Analysis Project (GAP), 2018, Clapper Rail (Rallus longirostris) bCLRAx_CONUS_2001v1 Habitat Map: U.S. Geological Survey data release, https://doi.org/10.5066/F7959FV4.
range_CLRA<-st_read(paste0(dat_path,"Breeding_Ranges/CLRA_northern_breeding_range_limit.shp"))%>%
  st_transform(crs="EPSG:26918")
range_limits[range_limits$species=="CLRA","ymax"]<-ext(range_CLRA)[4]


write.csv(range_limits,paste0(path_out,"Intermediate_outputs/ME_VA_range_limits_focal_spp.csv"),row.names=F)
