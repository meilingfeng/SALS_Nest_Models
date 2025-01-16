
library(sf)
library(tidyverse)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(dismo)
source("C:/Users/10788/Desktop/SaltMarsh/SHARP/Functions/gridSample_sf.R")


##############################################################################################
#Start Species Specific Analysis Here. Summarize data availability after spatial data cleaning.
##############################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## 1. Format Nest Observations
# -------------------------------------
#Record nesting data summary for each species of interest
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
species_sum<-as.data.frame(speciesnames)
# load nest observation shapefile
nests_all<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  # remove records that have coordinate errors (records missing location data already removed)
  filter(crd_typ!=1&Year>=2010)



for (j in 1:length(speciesnames)){
  nests<-nests_all%>%
    # filter records based on species
    filter(Species==speciesnames[j])%>%
    # convert all coordinates to decimal degrees (NAD83) and create Long and Lat columns
    st_transform("EPSG:4269")%>%
    mutate(Long = sf::st_coordinates(.)[,1],
           Lat = sf::st_coordinates(.)[,2],
           # create binary nest success variable for remaining records
           fate=case_when(fate=="FLEDGED" ~ 1,
                          fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
    dplyr:: select(id,latitude=Lat,longitude=Long,fate,Year,site=site_cd)%>%
    distinct(id,.keep_all=T)
  if(speciesnames[j]=="NESP"){
    nests=nests[-273,] #remove one record at site: 38101_HRC081:] site code and id aren't match, causing issues in distance calculation. 
  }

## 2. Summary of data availability
#---------------------------------------------------------
  species_sum$nest_n[j]<-nrow(nests) 
  species_sum$fate_n[j]<-nrow(nests[!(is.na(nests$fate)),])
  species_sum$site_n[j]<-length(unique(nests_all$site))
  species_sum$year_n[j]<-list(sort(unique(nests$Year)))
  species_sum
  
#  speciesnames nest_n fate_n site_n                                                           year_n
#  1         SALS   2476   2301     42 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  2         SESP   1078   1036     42             2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019
#  3         CLRA    394    325     42 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021
#  4         WILL    357    284     42       2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  5         NESP    273    271     42       2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  6         HYBR    100     99     42                                     2011, 2012, 2013, 2014, 2015
  
  ## 3. write a new file with only nests past 2010 to assume consistent sampling protocol
  #--------------------------------------------------------------------------------
  st_write(nests,paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_nests_2010_2020.shp"),append = FALSE)
}
