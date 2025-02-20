library(tidyverse)
library(dismo)
library(raster)
library(terra)
library(sf)
library(lubridate)


#file path names
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#objects
#Focal species list
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")

#Focal species breeding range limits within the study region
range_limits<-read.csv(paste0(path_out,"Intermediate_outputs/ME_VA_range_limits_focal_spp.csv"))

#list with 30m res tidal marsh layers
load(paste0(path_out,"/predictor_files_all_zones_30m.rds"))

# load all species nest observations
sampling<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_29_25.shp"))%>%
  # filter records missing coordinate information or that have coordinate errors
  filter(ot_bnds!=1&wrng_st!=1&iso_rec!=1&mssng_l_!=1&mssng_c!=1)%>%
  # create Long and Lat columns
  st_transform("EPSG:26918")



## Select random background points in each geographic zone for each focal species
#----------------------------------------------------------------------------
# a) get locations where nests are present for each species
for (j in 1:length(speciesnames)){
  # If background points don't already exist for the species, compute the following.
  #if(!file.exists(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"))){
  
    #reset object to hold bg points from each zone 
    bg<-list()
    
  # load cleaned focal species nest observation shapefile
  nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))
  
  
  
    
# b) bring in a mask of the tidal marsh area in 8 geographic zones along the Atlantic coast
  masks<-list()
  for(i in 1:length(file_list_all_zones)){
  masks[[i]]<-rast(file_list_all_zones[[i]][[1]])
  }
  
  
  
# c) set cells outside the marsh boundary to NA (do not select these as background locations)
  for(i in 1:length(masks)){
    mask<-masks[[i]]
    
    #raster values for upland and open water. Set these as outside the marsh boundary. 
    mask[mask==0|mask==9|mask==7]<-NA 
  
    #store marsh boundary mask
    masks[[i]]<-mask
    
  }  
  
#check NA count
global(masks[[3]],fun="isNA")
  

masks2<-masks
# d) Give each cell a probability of being selected for background points based on SHARP sampling intensity
  # give all non-NA cells a probability of being randomly selected as a background location
  # cells that have more nest observations (assumed to be greater sampling intensity) having greater weight (probability of selection)
  for(i in 1:length(masks2)){
    mask<-masks2[[i]]
    # set all cells with no records to default of 0.01 probability of selection (lowest probability)
    mask[cells(mask)]<-0.01
    # Count the total nest observations in each cell
    weights<-terra::extract(mask, sampling, cell=T)[,3]
    #cells that have nest observations
    weights<-data.frame(cells=weights[!is.na(weights)])%>%
      group_by(cells)%>%
      #counts the number of nests in each cell
      summarise(n=n())%>%
      ungroup()%>%
      # divide nest count by max nest count to standardize weights to 0-1 scale (probabilites)
      mutate(weight=n/max(n))%>%
      # filter for cells that are not NA (cells(raster) returns cell numbers of non-NA cells)
      filter(cells%in%cells(mask))%>%
      select(-n)
    # give all cells with records their standardized nest count
    mask[weights$cells]<-weights$weight
    # give the same nest count to all adjacent cells
    adj_cells<-as.data.frame(adjacent(mask,weights$cells,directions="queen", pairs=T, include=T, symmetrical=T))%>%
      left_join(weights,by=c("V1"="cells"))%>%
      left_join(weights,by=c("V2"="cells"))%>%
      filter(if_any(c(weight.x,weight.y),is.na))
    adj_cells1<-select(adj_cells,cells=V1,weight=weight.y)%>%
      filter(!is.na(weight))
    adj_cells2<-select(adj_cells,cells=V2,weight=weight.x)%>%
      filter(!is.na(weight))
    adj_cells3<-rbind(adj_cells1,adj_cells2)
    mask[adj_cells3$cells]<-adj_cells3$weight
    
    #store probability mask
    masks2[[i]]<-mask
  }

#check NA count
global(masks[[3]],fun="isNA")
  


# d) set cells outside the species breeding range and cells with species' nest observations to NA (exclude from background point selection)
  masks3<-masks2
  
  for(i in 1:length(masks3)){
    mask<-masks3[[i]]
   
    #cells that contain nests
    mask[unique(terra::extract(mask, nests, cell=T)[,3])]<-NA
    
    # set cells outside the species' breeding range to NA
if(any(!is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))){
  # if the species has a southern latitude limit
        # and the limit extends into the geographic zone, mask where it stops in that zone
  if(which(!is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))==1 &
          # highest latitude in the zone is NOT lower than lowest latitude of range
     !(ext(mask)[4]<range_limits[range_limits$species==speciesnames[j],]$ymin)){
    mask<-terra::mask(mask,ext(ext(mask)[1:2],range_limits[range_limits$species==speciesnames[j],]$ymin,ext(mask)[4]), touches=T, updatevalue=NA)
  }
        # if the limit does not reach the zone, mask the entire zone
  if(which(!is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))==1 &
          # highest latitude in the zone IS lower than lowest latitude of range
     ext(mask)[4]<range_limits[range_limits$species==speciesnames[j],]$ymin){
    values(mask)<-NA
  }
  # if the species has a northern latitude limit
        # and the limit extends into the geographic zone, mask where it stops in that zone
  if(which(!is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))==2&
          # lowest latitude in the zone is NOT higher than highest latitude of range
     !(ext(mask)[3]>range_limits[range_limits$species==speciesnames[j],]$ymax)){
    mask<-terra::mask(mask,ext(ext(mask)[1:3],range_limits[range_limits$species==speciesnames[j],]$ymax), touches=T, updatevalue=NA)
  }
        # if the limit does not reach the zone, mask the entire zone
  if(which(!is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))==2 &
          # lowest latitude in the zone IS higher than highest latitude of range
     ext(mask)[3]>range_limits[range_limits$species==speciesnames[j],]$ymax){
    values(mask)<-NA
  }
}
  # if the species has no latitude limits in the study region, don't mask anything
if(all(is.na(range_limits[range_limits$species==speciesnames[j],c("ymin","ymax")]))){
    mask<-mask
}
    
    masks3[[i]]<-mask
  }
    #check NA count
global(masks3[[3]],fun="isNA")

  
  
# e) select n random background points in each zone that equal the total number of nest records that zone
    #first assign each nest to a geographic zone (base on which correll raster zone layer they overlap)
  nests$Region<-NA
      #for each zone...
  for(i in 1:length(masks3)){
      #nests that are within the zone, assign them that zone's region number
    nests[which(!is.nan(terra::extract(masks3[[i]], nests,cells=T)[,3])),"Region"]<-i
  }
      #fill in region for any nests that plot out of the mask based on regions that other records in their site were assigned
  nests<-ungroup(mutate(group_by(nests,site_cd),Region=max(Region,na.rm=T)))
  
  # total number of nests in each zone
  n_nests<-st_drop_geometry(summarise(group_by(nests,Region),n=n()))
  # increase the count by a little bit (1%) in case some nest locations are on the raster edge
  n_nests$n<-round(n_nests$n+(n_nests$n*0.01))
  
  #if a zone is within the species range, but there are no records in that zone, give it the minimum observation count across zones.
  temp<-list()
  for(i in 1:length(masks3)){
  if(length(cells(masks3[[i]]))>0 & !(i%in%n_nests$Region)){
    temp[[i]]<-data.frame(Region=i,n=min(n_nests$n))
  }
  }
  new_regions<-do.call("rbind",temp)
  
  #remove the extra bg points added from the other zones to keep the bg point balanced with the presence count (if there are a lot).
  if((length(new_regions)*min(n_nests$n))>18){
  n_nests<-n_nests%>%
    mutate(n=n-((length(new_regions)*min(n))/nrow(n_nests)))
  }
  
  n_nests2<-rbind(n_nests,new_regions)
  
  

  #for each region...
  for(i in n_nests2$Region){
    mask<-masks3[[i]]
  # set parameters: 
    # adjustment that accounts for the total non-NA cells to sample from (number of available cells divided by number of needed points)
    #tf<-floor(ncell(mask[cells(mask)])/n_nests[n_nests$Region==i,]$n)
    # seed to ensure the same random sample each time
    set.seed(1963)
  # select N random background points
  bg[[i]] <- as.data.frame(as.matrix(raptr::randomPoints(mask=mask, n_nests2[n_nests2$Region==i,]$n, prob = T)))%>% 
    mutate(Region=i,
          bp="b") #mark as a background location
  }
  # bind all the background points across zones together
  bg_all<-do.call("rbind",bg)

  
#write coordinates to csv
write.csv(bg_all%>%mutate(Year=NA,site=NA,fate=NA,id=paste0("b",c(1:nrow(.)),"r",Region)),
          paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"),row.names=F)

#}
}
