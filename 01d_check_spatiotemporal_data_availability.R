
library(sf)
library(tidyverse)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(dismo)
source("C:/Users/mefen/OneDrive/Documents/Github/UCONN/SHARP/Functions/gridSample_sf.R")


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
nests_all<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  # remove records missing coordinate information or that have coordinate errors
  filter(crd_typ!=1&mssng_l_!=1&mssng_c!=1&Year>=2010)
  
nests<-nests_all%>%
  # filter records to just SALS or another species of interest
  filter(Species=="SALS")%>%
  # convert all coordinates to decimal degrees (NAD83) and create Long and Lat columns
  st_transform("EPSG:4269")%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
  # create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                        fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(id,latitude=Lat,longitude=Long,fate,Year,site=site_cd)%>%
  distinct(id,.keep_all=T)




## 2. Summary of data availability
#---------------------------------------------------------
nrow(nests) 
#There are 3039 nest observations
nrow(nests[!(is.na(nests$fate)),])
#There are 2807 nest fate observations

length(unique(nests_all$site))
#53 sites

sort(unique(nests$Year))
#2006-2020, we may want to just use data from 2010 and later to assume standardized sampling protocol

nests<-nests%>%filter(Year>2009)
nrow(nests) 
#There are 2476 nest observations since 2010
nrow(nests[!(is.na(nests$fate)),])
#There are 2301 nest fate observations since 2010

length(unique(nests$site))
#31 sites
t<-summarise(group_by(nests,site),n=n())




## 3. write a new file with only nests past 2010 to assume consistent sampling protocol
#--------------------------------------------------------------------------------
st_write(nests,paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020.shp"))

