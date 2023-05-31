library(tidyverse)
library(dismo)
library(raster)
library(sf)
library(lubridate)

reso<-30
#file path names
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"
bg<-list()

#functions
#source("C:/Users/mefen/OneDrive/Documents/Github/UCONN/SHARP/Functions/randomPoints_terra_sf.R")




## 1. Use random veg plots in each zone
#----------------------------------------------------------------------------
veg<-st_read(paste0(path_out,"Final_outputs/Veg_locations/veg_locations_12_29_22.shp"))%>%
  filter(PontTyp=="Random",
         # also filter records missing coordinate information or that have coordinate errors
         crd_typ!=1)%>%
  st_transform("EPSG:4269")%>% #"EPSG:26918"
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         bp="v",
         Year=year(mdy(date)),
         fate=NA)%>% #mark as a random veg location
  #st_transform("EPSG:26918")%>%
  dplyr:: select(id=veg_id,latitude=Lat,longitude=Long,site=Site,bp,Year,fate)

if(!file.exists(paste0(path_out,"Intermediate_outputs/random_veg_points.csv"))){
write.csv(st_drop_geometry(veg),paste0(path_out,"Intermediate_outputs/random_veg_points.csv"),row.names=F)
}





## 2. Select random background points in each zone
#----------------------------------------------------------------------------
# if bg points dont exist, compute. Otherwise load the file.
if(!file.exists(paste0(path_out,"Intermediate_outputs/background_points_",reso,"m.csv"))){
  
#abring in mask of marsh area for each zone
load(paste0(path_out,"/predictor_files_all_zones_",reso,"m.rds"))
  masks<-list()
  for(i in 1:length(file_list_all_zones)){
  masks[[i]]<-raster(file_list_all_zones[[i]][[1]])
  }
  
  # a) bring in mask of marsh area for each zone 
  #list raster files
  file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))
  #read as raster layers
  masks<-map(file_list,rast)
  
  # b) set area outside marsh to NA
  for(i in 1:length(masks)){
    mask<-masks[[i]]
    mask[mask==0|mask==9|mask==7]<-NA
    masks[[i]]<-mask}
  # c) get proportion of unoccupied marsh area for each zone to assign number of random points
  cell_prop<-c()
  # c) set cells with nests to NA
  for(i in 1:length(masks)){
    mask<-masks[[i]]
    mask[which(!is.na(extract(mask, nests,cells=T)))]<-NA
    cell_prop[i]<-ncell(mask[!is.na(mask)])
    masks[[i]]<-mask
  }
  
  # d) Get a count of total SHARP survey points in each cell to assign sampling effort weights
  # load nest observation shapefile
  sampling<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
    # filter records missing coordinate information or that have coordinate errors
    filter(crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
    # convert all coordinates to degrees (NAD83) and create Long and Lat columns
    st_transform("EPSG:4269")
  
  weights<-list()
  
  for(i in 1:length(masks)){
  weights<-extract(mask, sampling,cells=T)
    
    
# select n random points out of 10,000 relative to the proportional area of each zone
  cell_prop<-cell_prop/sum(cell_prop)
n<-round(cell_prop*10000)
sum(n)

# set seed to assure that the examples will always
# have the same random sample.
for(i in 1:length(masks)){
set.seed(1963)
bg[[i]] <- as.data.frame(as.matrix(randomPoints(mask=masks[[i]], n[i])))%>% #tryf=50
  mutate(region=i,
         bp="b") #mark as a background location
}

bg_all<-do.call("rbind",bg)

#And inspect the results by plotting

# set up the plotting area for two maps
par(mfrow=c(1,2))
plot(!is.na(masks[[6]]), legend=FALSE)
points(bg[[6]], cex=0.5)

#write coordinates to csv
write.csv(bg_all,paste0(path_out,"Intermediate_outputs/background_points_30m.csv"),row.names=F)
}

#load background point coordinates
bg_all<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_30m.csv"))

#make background points sf shapefile
bg_points<-st_as_sf(bg_all,coords=c("x","y"),crs="EPSG:26918")%>%
  st_transform("EPSG:4269")%>%
  mutate(longitude=sf::st_coordinates(.)[,1],
         latitude=sf::st_coordinates(.)[,2],
         Year=NA,
         site=NA,
         fate=NA,
         id=paste0("b",c(1:nrow(.)),"r",region))%>%
  dplyr::select(-region)



