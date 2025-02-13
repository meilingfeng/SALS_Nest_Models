library(tidyverse)
library(sf)
library(terra)#updated version of raster package
library(tidyterra)
library(dismo)
library(gbm)
#source("C:/Users/10788/Desktop/SaltMarsh/SHARP/Functions/gridSample_sf.R")
#https://rspatial.org/sdm/1_sdm_introduction.html

### Set up
# -------------------------------------------
##  file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## all species to run through
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")


## list of predictor surface files for each zone
load(paste0(path_out,"predictor_files_all_zones_30m.rds"))
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 


## Load nest observation datasets with predictor variables
for (j in 1:length(speciesnames)){
  # predictors at nest location point
  dat_point<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_vars_local.csv"))%>%
    dplyr::select(id,site,Year,bp,tideres)
  
  
  # buffered predictors - measures summarized within a 30 m buffer of nest
  dat_buff<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_vars_buff15.csv"))%>%
    dplyr::select(id,site,Year,bp,uvvr_mean,ndvi,pca,HIMARSH,LOMARSH,uvvr_diff,elevation)
  
  #combine point and buffer variables
  dat<-left_join(dat_point,dat_buff,by=c("id","site","Year","bp"))%>%
    mutate(presence=ifelse(bp=="p",1,0))%>%
    dplyr::select(-bp)
  
  #add fate variable, geographic region/zone, latitude
  fates<-read.csv(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_nests2024_regions.csv"))%>%
    dplyr::select(id,fate,Region)%>%
    mutate(fate=case_when(
      fate=="FLEDGED"~1,
      fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED")~0,
      is.na(fate)|!(fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED"))~NA))%>%
    distinct(id,.keep_all = T)
  lat<-st_read(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))%>%
    st_transform("EPSG:4269")%>%
    dplyr::mutate(latitude = sf::st_coordinates(.)[,2],
                  longitude = sf::st_coordinates(.)[,1])%>%
    dplyr::select(id,latitude,longitude)
  fates<-fates%>%
    left_join(lat,by="id")
  
  bg_regions<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"))%>%
    st_as_sf(coords=c("x","y"),crs="EPSG:26918")%>%
    st_transform("EPSG:4269")%>%
    dplyr::mutate(latitude = sf::st_coordinates(.)[,2],
                  longitude = sf::st_coordinates(.)[,1])%>%
    dplyr::select(id,Region,latitude,longitude)%>%
    mutate(fate=NA)

  fates<-rbind(fates,bg_regions)
  
  
  dat<-dat%>%
    left_join(fates,by="id")%>%
    dplyr::select(-geometry)
  
  # replace NA land cover proportions with 0 
  dat[,c("HIMARSH","LOMARSH")]<-dat[,c("HIMARSH","LOMARSH")]%>%
    replace(is.na(.),0)
  
  

  
### Add some additional nuisance variables for GLMs
#--------------------------------------------------------------------------
  
## 1. Days since the last spring tide the nest was initiated
  # More time since last moon = less time to fledge before the nest highest high tide = less successful nests 
time_since_tide<-read.csv(paste0(path_out,"Intermediate_outputs/Tides/time_since_tide_",speciesnames[j],".csv"))

dat2<-dat%>%
  left_join(time_since_tide,by="id")

  #nest records missing initiation dates (75)- nests without monitoring data
test<-dat2%>%
  filter(presence==1&is.na(init_date))
nrow(test)


## 2. Day of the year
  # Expect nest site selection to change as vegetation grows in.
dat3<-dat2%>%
  mutate(doy=yday(init_date))

  #fill in day of year for background points by assigning each bp to a nest initiation date
date_order<-list()

    # draw random index numbers for all the nest initiation dates in each region
regions<-sort(unique(dat3[dat3$presence==1,]$Region))
for(i in 1:length(regions)){
  date_order[[i]]<-sample(which(dat3$presence==1&!is.na(dat3$init_date)&dat3$Region==regions[[i]]),
                          nrow(dat3[dat3$presence==0 & dat3$Region==regions[[i]], ]),
                          replace = T)
  
  dates<-vector()
  
  for(k in 1:length(date_order[[i]])){
  dates[k]<-dat3$init_date[date_order[[i]][k]]
  }
  
  dat3[dat3$presence==0 & dat3$Region==regions[[i]], ]$init_date<-dates
}


  # for background points in regions without nest observations, assign them a random initiation date from the nest data
dates2<-sample(which(dat3$presence==1&!is.na(dat3$init_date)),
       nrow(dat3[dat3$presence==0&is.na(dat3$init_date),]),
       replace = F)
dat3[dat3$presence==0&is.na(dat3$init_date),]$init_date <-dat3$init_date[dates2]


dat3<-dat3%>%
  mutate(doy=yday(init_date))


  
### Filter for complete records (necessary for GLMs)
#--------------------------------------------------------------------------------
  # filter for remote sensing data
  dat_comp<-dat3[complete.cases(dat3[,all_terms]),]
  nrow(dat3[!(dat3$id%in%dat_comp$id),]) #removes 380 background points


  #remove nest records without initiation dates. 
  dat_comp2<-dat_comp[!is.na(dat_comp$doy),]
  nrow(dat_comp)-nrow(dat_comp2) # removes 75 nest records (only 65 have nest fates)
  
  
  
  
###  Divide observations into nest presence and nest survival datasets
#----------------------------------------------------------------------------------
  #for GLM methods (N=5625 presence data (2949 nest locations) and N= 2613 fate data) 
  #(only nest records removed are those without monitoring data/initiation dates)
  pres_dat<-dat_comp2
  sum(pres_dat$presence)
  surv_dat<-dat_comp2[!is.na(dat_comp2$fate),] #filter just the nests with known fates
  surv_dat<-surv_dat[!(is.na(surv_dat$time_since_tide)),]
  
  #for the ML methods (N=6080 presence data (3024 nest locations) and N= 2678 fate data)
  pres_dat_ml<-dat3
  sum(pres_dat_ml$presence)
  surv_dat_ml<-dat3[!is.na(dat3$fate),]
  
  
  
  
  ## look at the prevalence of:
  
  # nest successes vs fails 
  # complete data (GLMs) (For SALS- 1660:953, 0.36 prevalence of success, total records =2613)
  addmargins(table(surv_dat$fate));sum(surv_dat$fate)/nrow(surv_dat)
  # all data (ML) (For SALS 1703:975, 0.36 prevalence of success, total records =2678)
  addmargins(table(surv_dat_ml$fate));sum(surv_dat_ml$fate)/nrow(surv_dat_ml)
  
  # nest presence to background 
  # complete data (GLMs) (For SALS 2676b:2949p, 0.52)
  addmargins(table(pres_dat$presence)); sum(pres_dat$presence)/nrow(pres_dat)
  #all data (ML) (For SALS 3056b:3024p, 0.50)
  addmargins(table(pres_dat_ml$presence)); sum(pres_dat_ml$presence)/nrow(pres_dat_ml)
  
  

  
### Adjust class balance (if needed)
#---------------------------------------------------------------------------------------------

# Spatially thin to balance majority class in the ML dataset if they are unbalanced
#in general, balanced datasets are more important than maintaining true prevalence when discriminating between sites is the goal and ML approach is to be used (Steen; Barbet)
# for GLMs data availability is comparatively more important, especially for rare species (prevalence <0.1) (Steen; Barbet), so only worry about the ML dataset
#make sure we maintain at least 100 records of each class
if(sum(surv_dat_ml$fate)/nrow(surv_dat_ml)<0.4&nrow(surv_dat_ml[surv_dat_ml$fate==1,])>100){
  res_list<-c(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
    # create a Raster Layer with the extent of nest locations
    bg<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"))
    nest_extent<-vect(st_as_sf(bg,coords=c("x","y"),crs="EPSG:26918"))
    points<-vect(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))%>%
      # project to UTM so we can specify resolution in meters
      project("EPSG:26918")%>%
      tidyterra::mutate(fate=case_when(
        fate=="FLEDGED"~1,
        fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED")~0,
        is.na(fate)|!(fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED"))~NA))
    
    
    for(i in 1:length(res_list)){
    # set the resolution of the cells to 30 m
    res_thin<- res_list[i]
    r <- rast(extent=ext(nest_extent), resolution = res_thin, crs="EPSG:26918")
    
    # thin nest fails down to 1 nest per 30 meters
    fails<-points%>%
      tidyterra::filter(fate==0)
    
    non_fails<-points%>%
      tidyterra::filter(fate!=0|is.na(fate))
    
    set.seed(123)
    nest_thin <- spatSample(fails, size=1, strata=r)#strata are the 30m raster cells. Sample 1 point per cell.
    
    fails_ss<- filter(fails,id%in%nest_thin$id)# select only the thinned points
    
    #number of nest records removed N=622
    nrow(fails)-nrow(fails_ss)
    
    #recombine subset fails and fledges
    nests_ss<-rbind(fails_ss,non_fails)
    
    
    #balance survival data by filtering for the thinned record ids
    surv_dat_ml2<-surv_dat_ml%>%
      # all the successes and any fails that made it through the thinning
      filter(fate==1|id%in%fails_ss$id)
    
    # get close to 0.5 prevelance
    if(sum(surv_dat_ml2$fate)/nrow(surv_dat_ml2)>=0.4) break()
    }
    
    surv_dat_ml<-surv_dat_ml2
}

table(surv_dat_ml$fate)


if(sum(surv_dat_ml$fate)/nrow(surv_dat_ml)>0.6&nrow(surv_dat_ml[surv_dat_ml$fate==0,])>100){
  res_list<-c(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
    # create a Raster Layer with the extent of nest locations
    bg<-read.csv(paste0(path_out,"Intermediate_outputs/background_points_",speciesnames[j],"_30m_01_31_2025.csv"))
    nest_extent<-vect(st_as_sf(bg,coords=c("x","y"),crs="EPSG:26918"))
    points<-vect(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))%>%
      # project to UTM so we can specify resolution in meters
      project("EPSG:26918")%>%
      tidyterra::mutate(fate=case_when(
        fate=="FLEDGED"~1,
        fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED")~0,
        is.na(fate)|!(fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED"))~NA))
    
    
    for(i in 1:length(res_list)){
        # set the resolution of the cells to 30 m
        res_thin<- res_list[i]
    r <- rast(extent=ext(nest_extent), resolution = res_thin, crs="EPSG:26918")
    
    # thin nest successes down to 1 nest per 30 meters
    fails<-points%>%
      tidyterra::filter(fate==0)
    
    non_fails<-points%>%
      tidyterra::filter(fate!=0|is.na(fate))
    
    set.seed(123)
    nest_thin <- spatSample(non_fails, size=1, strata=r)#strata are the 30m raster cells. Sample 1 point per cell.
    
    non_fails_ss<- filter(non_fails,id%in%nest_thin$id)# select only the thinned points
    
    #number of nest records removed N=622
    nrow(non_fails)-nrow(non_fails_ss)
    
    #recombine subset fails and fledges
    nests_ss<-rbind(fails,non_fails_ss)
    
    
    #balance survival data by filtering for the thinned record ids
    surv_dat_ml2<-surv_dat_ml%>%
      # all the fails and any successes that made it through the thinning
      filter(fate==0|id%in%non_fails_ss$id)
    

    # get close to 0.5 prevelance
    if(sum(surv_dat_ml2$fate)/nrow(surv_dat_ml2)<=0.6) break()
    }
    surv_dat_ml<-surv_dat_ml2
}
table(surv_dat_ml$fate)
 
  
  ## Summary of data availability after filtering and thinning/balancing
  #-----------------------------------------------------------------------------------
  nrow(surv_dat_ml)
  #There are 2271 nest fate observations after spatial thinning
  
  
  
  
  ## Divide into testing and training data (80 train/20 test) - k-fold data partitioning
  #---------------------------------------------------------------------------------------------
  
  # if doing presence absence, split both. If doing presence only (bioclim), only split the presence data, background data only used for testing
  
  #number of partitions to make. eg, 5 splits into 5 groups each of which will be used once as a test dataset (20% testing 80% training)
  k<-5
  set.seed(123)
  pres_dat$group<-kfold(pres_dat,k)
  surv_dat$group<-kfold(surv_dat,k)
  pres_dat_ml$group<-kfold(pres_dat_ml,k)
  surv_dat_ml$group<-kfold(surv_dat_ml,k)
  
  
  #check balance in each split dataset
  table(pres_dat$presence,pres_dat$group) # stays mostly well balanced
  table(surv_dat$fate,surv_dat$group)
  table(pres_dat_ml$presence,pres_dat_ml$group) # stays mostly well balanced
  table(surv_dat_ml$fate,surv_dat_ml$group)
  
  
  # set a constant response variable name
  pres_dat$y<-pres_dat$presence
  surv_dat$y<-surv_dat$fate
  pres_dat_ml$y<-pres_dat_ml$presence
  surv_dat_ml$y<-surv_dat_ml$fate

  

  write.csv(pres_dat,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_pres_GLM_dat.csv"),row.names = F)
  write.csv(surv_dat,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_fate_GLM_dat.csv"),row.names = F)
  write.csv(pres_dat_ml,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_pres_ML_dat.csv"),row.names = F)
  write.csv(surv_dat_ml,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_fate_ML_dat.csv"),row.names = F)
}
