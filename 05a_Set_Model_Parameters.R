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
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## 2. parameters
# what prediction surface resolution?
reso<-30

# use random background points or veg plots as absences?
#ab_type<-"v"
ab_type<-"b" 

# veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))


## 3. Load observations with predictors
# list of predictor surface files for each zone
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 

speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
for (j in 1:length(speciesnames)){
  # point predictors - UVVR, tidal restriction
  dat1<-read.csv(paste0(path_out,"Final_outputs/",speciesnames[j],"_nest_vars_local.csv"))%>%
    filter(bp%in%c("p",ab_type))%>%
    mutate(presence=ifelse(bp=="p",1,0))%>%
    left_join(veg_codes,by="veg_class")%>%
    dplyr::select(-pca,-ndvi,-cor_txt,-ent_txt)
  
  # buffered predictors - high marsh proportion, avg texture, ndvi, pca, and elevation
  dat<-read.csv(paste0(path_out,"Final_outputs/",speciesnames[j],"_nest_vars_buff15.csv"))%>%
    dplyr::select(id,HIMARSH,LOMARSH,ndvi,pca,elevation)%>% #remove UPLND
    right_join(dat1,by="id")%>%
    mutate(veg_code=as.factor(veg_code),
           #create binary variable for if nests intersect High Marsh habitat
           Highmarsh=as.factor(ifelse(veg_class=="HIMARSH",1,0)))
  
  # fill in 0's for land cover proportions with NAs 
  dat[,c("HIMARSH","LOMARSH")]<-dat[,c("HIMARSH","LOMARSH")]%>%
    replace(is.na(.),0)
  
  
  # Remove observation records with missing predictor values (regression modeling methods cannot have missing values)
  dat_comp<-dat[complete.cases(dat[,all_terms]),]
  dat[!(dat$id%in%dat_comp$id),]
  
  
  # Add some additional nuisance variables
  #--------------------------------------------------------------------------
  # 1. Time nest was initiated since the last new moon - estimate of how late in the flood free window the nest was started
  # later nests, more time since last moon, have less time to fledge before the nest highest high tide. 
  moon<-read.csv(paste0(dat_path,"Environmental Predictors/new_moon_dates.csv"))%>%
    mutate(new.moon.date=mdy(new.moon.date))%>%
    arrange(new.moon.date)%>%
    mutate(lag=lead(new.moon.date))%>%
    filter(!is.na(lag))
  
  #create a list of potential dates that come after each moon date to join to the nest data
  moon.dates<-list()
  for(i in 1:nrow(moon)){
    moon.dates[[i]]<-seq(moon[i,]$new.moon.date,moon[i,]$lag,by="days")
  }
  
  #add to each new.moon.date
  moon.dat<-expand.grid(new.moon.date=moon$new.moon.date, date=do.call("c",moon.dates))
  
  moon.dat<-moon.dat%>%
    # subtract the new moon date from the survey date to get days since the new moon
    mutate(time_since_moon=date-new.moon.date)%>%
    # moon cycle is 29 days, so time should be more than 0 but less than 29 days
    filter(time_since_moon>=0&time_since_moon<=29)%>%
    distinct(.keep_all = T)%>%
    mutate(time_since_moon=as.numeric(time_since_moon))%>%
    dplyr::select(-new.moon.date)
  
  #get all possible first egg dates using nest monitoring data
  mon<-read.csv(paste0(dat_path,"Demographic Database/Nest_Monitoring_May_24.csv"))
  #replace all missing chick age values with NA
  mon[,10:25]<-replace(mon[,10:25],mon[,10:25]==""|mon[,10:25]=="NOT REC"|mon[,10:25]=="na",NA)
  
  #Calculate hatch dates based on chick age for records missing first egg dates.
  mon2<-mon%>%
    #filter the NA ages out
    filter(if_any(contains("Chick.Age"),~!is.na(.x)))%>%
    #keep only one observation per nest
    distinct(SHARP.Nest,.keep_all = T)%>%
    rowwise()%>%
    #calculate the oldest chick age for the nest
    mutate(max_age=as.numeric(max(c_across(contains("Chick.Age")),na.rm = T)),
           # estimate the hatch date by subtracting the oldest chick age from the survey date
           hatch.date=mdy(Visit.Date)-days(max_age))%>%
    dplyr::select(id=SHARP.Nest,hatch.date)%>%
    ungroup()
  
  #also get all records with hatch dates from nest fates data
  add.hatch<-read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2020.csv"))%>%
    dplyr::select(id=SHARPNestID,hatch.date=EstHatchDate)%>%
    filter(!(hatch.date%in%c("Unknown","UNK","NOT REC",""))&!is.na(hatch.date))%>%
    mutate(hatch.date=mdy(hatch.date))
  
  #remove any age estimated hatch dates that already have estimate dates in the fates data
  mon3<-mon2%>%
    filter(!(id%in%c(add.hatch$id)))
  
  #join all the estimates hatch dates
  all.add.hatch<-rbind(mon3,add.hatch)%>%
    # subtract 10 incubation days from the hatch date, this is the day the last egg was laid per Ruskin 2016 (Field 2018 Ecography paper)
    mutate(inc.date=hatch.date-days(10))
  
  # calculate total eggs
  mon<-mon%>%
    rowwise()%>%
    mutate(Eggs=as.numeric(Eggs),
           Nestlings=as.numeric(Nestlings))%>%
    mutate(total.eggs=sum(c(Eggs,Nestlings),na.rm=T))%>%
    ungroup()
  
  max.eggs<-read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2020.csv"))%>%
    dplyr::select(id=SHARPNestID,total.eggs=MaxNumEggs)%>%
    filter(!(total.eggs%in%c("UNKNOWN","UNK","NOT REC",""))&!is.na(total.eggs))%>%
    mutate(total.eggs=as.numeric(total.eggs))
  
  #remove any total egg counts that already have estimates in the fates data
  mon4<-mon%>%
    filter(!(SHARP.Nest%in%c(max.eggs$id)))%>%
    dplyr::select(id=SHARP.Nest,total.eggs)%>%
    distinct(id,.keep_all = T)
  #combine
  all.add.firstegg<-rbind(mon4,max.eggs)%>%
    right_join(all.add.hatch,by="id")%>%
    #estimated first egg date is last egg date-(N eggs-1)
    mutate(date2=inc.date-days(total.eggs-1))%>%
    dplyr::select(id,date2)
  
  
  dat_comp<- read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2020.csv"))%>%
    dplyr::select(id=SHARPNestID, date=EstFirstEggDate)%>%
    mutate(date=mdy(date))%>%
    #first add estimated first egg dates from fates data to the cleaned nesting data
    right_join(dat_comp,by="id")
  
  #how many observations does calculating first egg date fix? (323)
  missing<-dat_comp[is.na(dat_comp$date),]
  nrow(missing[missing$id%in%all.add.firstegg$id,])
  
  dat_comp2<- dat_comp%>%
    #then add the calculated first egg dates for those that are missing
    left_join(all.add.firstegg,by="id")%>%
    mutate(date=coalesce(date,date2))%>%#coalesce takes values from date unless they are NA, then takes from date2
    dplyr::select(-date2)
  
  dat_comp2<-dat_comp2%>%
    left_join(moon.dat,by="date")%>%
    distinct(id,.keep_all = T)
  
  # 2. Day of the year. Expect nest site selection to change as vegetation grows in.
  dat_comp2<-dat_comp2%>%
    mutate(doy=yday(date))
  
  
  
  
  
  # Check for complete records again
  #-------------------------------------------------------------------------------
  #filter only complete records, those with first egg dates and background points
  dat_comp3<-dat_comp2[!is.na(dat_comp2$doy)|dat_comp2$presence==0,]
  dat_comp3<-dat_comp3%>%dplyr::filter(!is.na(region))
  
  
  ### Check Class Balance
  #---------------------------------------------------------------------------------------------
  # balancing can help machine learning models 
  
  
  # randomly select the same number of background points as the filtered, complete nest records
  # sample the same sample size as each region
  reg_size<-as.vector(table(dat_comp3[dat_comp3$presence==1,]$region))
  regions<-sort(unique(dat_comp3[dat_comp3$presence==1,]$region))
  temp1<-list()
  
  set.seed(123)
  for(i in 1:length(reg_size)){
    temp1[[i]]<-dat_comp3[dat_comp3$region==regions[[i]]&dat_comp3$presence==0,]
    temp1[[i]]<-temp1[[i]][sample(1:nrow(temp1[[i]]), reg_size[i]), ]
  }
  temp1[[length(reg_size)+1]]<-dat_comp3[dat_comp3$presence==1,]
  dat_comp4 <- do.call("rbind", temp1)
  dat_comp4<-dplyr::filter(dat_comp4,!is.na(id))
  
  #fill in day of year for background points by assigning each bp to a nest observation date
  date_order<-list()
  #drawn random index numbers for all the dates in each region
  for(i in 1:length(reg_size)){
    date_order[[i]]<-sample(1:reg_size[i],reg_size[i],replace = F)
    dat_comp4[dat_comp4$presence==0 & dat_comp4$region==regions[[i]], ]$date<-dat_comp4[dat_comp4$presence==1 & dat_comp4$region==regions[[i]], ]$date[date_order[[i]]]
  }
  
  dat_comp4<-dat_comp4%>%
    mutate(doy=yday(date))
  
  
  ##  divide observations into nest presence and nest survival datasets
  pres_dat<-dat_comp4
  surv_dat<-dat_comp4[!is.na(dat_comp4$fate),] #filter just the nests with known fates
  surv_dat<-surv_dat[!(is.na(surv_dat$time_since_moon)),]#removes one 2021 observation
  
  
  ## look at the prevalence of:
  
  # nest fledges to nest fails
  addmargins(table(surv_dat$fate));sum(surv_dat$fate)/nrow(surv_dat)
  
  # prevalence of presence to background
  addmargins(table(pres_dat$presence)); sum(pres_dat$presence)/nrow(pres_dat)
  

  
  
  
  # Thin to balance nest failures if there are a lot of them (zeros)
  if(sum(surv_dat$fate)/nrow(surv_dat)>0.59|sum(surv_dat$fate)/nrow(surv_dat)<0.41){
    # create a RasterLayer with the extent of nest locations
    points<-st_read(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_nests_2010_2020_dist_err_removed.shp"))%>%

      # project to UTM so we can specify resolution in meters
      project("EPSG:26918")
      # set the resolution of the cells to 30 m
      res<- 30
      r <- rast(extent=ext(points), resolution = res, crs="EPSG:26918")
      
      # thin nest fails down to 1 nest per 30 meters
      fails<-points%>%
        tidyterra::filter(fate==0)
      
      non_fails<-points%>%
        tidyterra::filter(fate!=0|is.na(fate))
      
      nest_thin <- spatSample(fails, size=1, strata=r)#strata are the 30m raster cells. Sample 1 point per cell.
      
      fails_ss<- filter(fails,id%in%nest_thin$id)# select only the thinned points
      
      #number of nest records removed N=562
      nrow(fails)-nrow(fails_ss)
      
      #recombine subset fails and fledges
      nests_ss<-rbind(fails_ss,non_fails)
      
      
      #balance survival data by filtering for the thinned record ids
      surv_dat<-surv_dat%>%
        # all the successes and any fails that made it through the thinning
        filter(fate==1|id%in%fails_ss$id)
      
      table(surv_dat$fate);sum(surv_dat$fate)/nrow(surv_dat) # get close to 0.5 prevelance 
  }
  
  
  # thin the veg data if using as background points
  if(ab_type=="v"){
    
    # create a RasterLayer with the extent of veg locations
    points<-st_read(paste0(path_out,"Final_outputs/veg_locations/veg_locations_12_29_22.shp"))%>%
      filter(PontTyp=="Random")
    r <- rast(points)
    
    # project to UTM so we can specify resolution in meters
    r<-terra::project(r,"EPSG:26918")
    veg_m<-st_transform(points,"EPSG:26918")
    
    # set the resolution of the cells to 30 m
    res(r) <- 30
    # expand (extend) the extent of the RasterLayer a little
    r <- extend(r, ext(r)+10)
    
    # thin veg points down to 1 point per 30 meters
    fails<-veg_m
    veg_thin <- gridSample(fails, r, n=1)
    
    fails$Easting<-sf::st_coordinates(fails)[,1]
    fails$Northing = sf::st_coordinates(fails)[,2]
    veg_ss<- filter(fails,(Easting%in%veg_thin[,1]&Northing%in%veg_thin[,2]))%>%
      distinct(geometry,.keep_all = T)
    
    #number of nest records removed N=293
    nrow(fails)-nrow(veg_ss)
    
    #balance survival data by filtering for the thinned record ids
    pres_dat<-pres_dat%>%
      filter(id%in%veg_ss$veg_id | bp=="p")
    
    table(pres_dat$presence);sum(pres_dat$presence)/nrow(pres_dat) 
    
  }
  
  
  
  ## Summary of data availability after filtering and thinning/balancing
  #-----------------------------------------------------------------------------------
  nrow(pres_dat) 
  #There are 3230 nest observations
  nrow(surv_dat)
  #There are 1528 nest fate observations (now 1528)
  
  length(unique(pres_dat$site))
  #32 sites
  
  
  
  
  ## Divide into testing and training data (80 train/20 test) - k-fold data partitioning
  #---------------------------------------------------------------------------------------------
  
  # if doing presence absence, split both. If doing presence only (bioclim), only split the presence data, background data only used for testing
  
  #number of partitions to make. eg, 5 splits into 5 groups each of which will be used once as a test dataset (20% testing 80% training)
  k<-5
  set.seed(123)
  pres_dat$group<-kfold(pres_dat,k)
  surv_dat$group<-kfold(surv_dat,k)
  
  
  #check balance in each split dataset
  table(pres_dat$presence,pres_dat$group) # stays mostly well balanced
  table(surv_dat$fate,surv_dat$group)
  
  
  
  # set a constant response variable name
  pres_dat$y<-pres_dat$presence
  
  surv_dat$y<-surv_dat$fate
  
  
  table(pres_dat$y)
  table(surv_dat$y)
  write.csv(pres_dat,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_pres_dat.csv"),row.names = F)
  write.csv(surv_dat,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_surv_dat.csv"),row.names = F)
  
}
