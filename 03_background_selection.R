library(tidyverse)
library(dismo)
library(raster)
library(terra)
library(sf)
library(lubridate)

reso<-30
#file path names
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"
bg<-list()

#functions
#source("C:/Users/mefen/OneDrive/Documents/Github/UCONN/SHARP/Functions/randomPoints_terra_sf.R")

### NOTE #####
# run this code through the source in 04a_Format_Load_Data_Nests_and Predictors.R


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
  
# a) bring in mask of marsh area for each zone
load(paste0(path_out,"/predictor_files_all_zones_",reso,"m.rds"))
  masks<-list()
  for(i in 1:length(file_list_all_zones)){
  masks[[i]]<-raster(file_list_all_zones[[i]][[1]])
  }


  
  
# b) set area outside marsh and cells with nests to NA
  for(i in 1:length(masks)){
    mask<-masks[[i]]
    mask[mask==0|mask==9|mask==7|mask==8]<-NA
    mask[unique(raster::extract(mask, nests,cell=T)[,1])]<-NA
    masks[[i]]<-mask
    }

# c) Give each cell a probability of being selected for background
  # load nest observation shapefile
  sampling<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
    # filter records missing coordinate information or that have coordinate errors
    filter(crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
    # convert all coordinates to degrees (NAD83) and create Long and Lat columns
    st_transform("EPSG:4269")
  
  # give each non NA cell a probability of random selection
  # cells that have more SHARP nest points having greater weight (probability of selection)
  for(i in 1:length(masks)){
    mask<-masks[[i]]
  # Get a count of total SHARP survey points in each cell to assign sampling effort weights
  weights<-raster::extract(mask, sampling,cell=T)[,1]
  weights<-data.frame(cells=weights[!is.na(weights)])%>%
    group_by(cells)%>%
    #counts the number of SALS nests in each cell
    summarise(n=n())%>%
    ungroup()%>%
    # divide nest count by max nest count to standardize weights to 1 (probabilites)
    mutate(weight=n/max(n))
    # set all cells with no records to default of 0.01 probability of selection
  mask[!(is.na(mask))]<-0.01
    # weight all others by their standardized record count
  mask[weights$cells]<-weights$weight
  masks[[i]]<-mask
  }
    
# d) select n random points that equal the total number of nest records in each zone
    #first assign each nest to a zone (which points overlap each correll raster zone layer)
  nests$Region<-NA
  nests$site<-substr(nests$id,1,2)
  for(i in 1:length(masks)){
    nests[which(!is.na(extract(masks[[i]], nests))),"Region"]<-i
  }
  nests<-ungroup(mutate(group_by(nests,site),Region=max(Region,na.rm=T)))
  nests[nests$site=="PA",]$Region<-3
  
  # total number of nests in each zone is total background points for each zone
  n<-st_drop_geometry(summarise(group_by(nests,Region),n=n()))
  # add a few extra points in case some are on the raster edge and dont sample finer res data
  n$n<-round(n$n+(n$n*0.1))

  
  for(i in n$Region){
    mask<-masks[[i]]
  # adjustment that accounts for the total non-NA cells to sample from
    tf<-round(ncell(mask[!is.na(mask)])/n[n$Region==i,]$n)
  # set seed to assure that the examples will always have the same random sample.
  set.seed(1963)
  bg[[i]] <- as.data.frame(as.matrix(randomPoints(mask=masks[[i]], n[n$Region==i,]$n, prob = T,tryf=tf)))%>% 
    mutate(region=i,
          bp="b") #mark as a background location
  }

  bg_all<-do.call("rbind",bg)

#And inspect the results by plotting

# set up the plotting area for two maps
par(mfrow=c(1,2))
plot(!is.na(masks[[4]]), legend=FALSE)
points(bg[[4]], cex=0.5)

#write coordinates to csv
write.csv(bg_all%>%mutate(Year=NA,site=NA,fate=NA,id=paste0("b",c(1:nrow(.)),"r",region)),
          paste0(path_out,"Intermediate_outputs/background_points_30m.csv"),row.names=F)
}

#load background point coordinates
bg_all<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_30m.csv"))

#make background points sf shapefile
bg_points<-st_as_sf(bg_all,coords=c("x","y"),crs="EPSG:26918")%>%
  st_transform("EPSG:4269")%>%
  mutate(longitude=sf::st_coordinates(.)[,1],
         latitude=sf::st_coordinates(.)[,2])%>%
  dplyr::select(-region)



