
library(tidyverse)
library(sf)
library(tmap)
library(terra)
library(tidyterra)
library(exactextractr)



## Set file path to data
# -------------------------------------------
dat_path<-"G:/My Drive/Research/SHARP/Data/"


## Load in data
# -------------------------------------------
#Set all coordinate systems to ESPG: 26918 (NAD 1983 UTM zone 18N)

## Environmental Predictor 1: Marsh Vegetation Classes 
  #list raster files
  file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

  #read as raster layers
  vg_cls<-map(file_list,rast)

  #set coordinate systems
  for(i in 1:length(vg_cls)){
  crs(vg_cls[[i]])<-"EPSG:26918"
  }

  #get cell area
  area<-prod(res(vg_cls[[1]]))

  #vegetation class definitions
  veg_class<-data.frame(value=c(1:2,4:9),
                        veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))


  

## Enviromental Predictor 2: Unvegetated-Vegetated ratio
  # transformed uvvr coordinate system from WGS to NAD83 in ArcPro
  # just want band 4 (UVVR). Band 1 is unveg, 2 is veg, and 3 is proportion of water.
  
  # A. UVVR mean across 2014-2018
  uvvr<-rast(paste0(dat_path,"UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything()) 
  
  # B. UVVR change (2014-2018)
  uvvr14<-rast(paste0(dat_path,"UVVR/uvvr_14_utm18_mean_ext.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_14'},.cols = everything())
    # values above 2 are not accurate, set to NA
    uvvr14[uvvr14>2]<-NA
  uvvr18<-rast(paste0(dat_path,"UVVR/uvvr_18_utm18_mean_ext.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_18'},.cols = everything())
    uvvr18[uvvr18>2]<-NA
    # subtract 2018 from 2014 to get change
  uvvr_diff<-uvvr18-uvvr14
  uvvr_diff<-uvvr_diff%>%
    dplyr::rename_with(function(x){x<-'uvvr_diff'},.cols = everything())


  
  
## Response Variable: Probability of nest success
  #load original data with all nest fate variables and observations
  nest_fates<-read.csv(paste0(dat_path,"NestFates_2001-2020.csv"))%>%
  #select just the unique nest observation ID and it's fate
    dplyr::select(id=SHARPNestID,fate=UltimateNestFate)%>%
  #remove nests missing fates (~500)
  filter(!(fate%in%c("UNKNOWN", "UNKNOWN IF FLEDGED OR FAILED", "NA","NOT MONITORED","NEVER HAD EGGS","NEVER HAD")))%>%
  #create binary nest success variable for remaining records
  mutate(fate=as.factor(ifelse(fate=="FLEDGED",1,0)))

  #load point layer of spatially projected nest records
  nests<-st_read(paste0(dat_path,"Nest Locations/nest_locations_12_9_22.shp"))%>%
  #filter records to just SALS
    filter(Species=="SALS"&crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  #remove unneeded variables
    select(-c("Lat","Long","crd_typ","mssng_l_","mssng_s_","mssng_c","fate"))%>%
  #project to UTM 18
    st_transform("EPSG:26918")
  nests_buff<- nests%>%
    #buffer nest points by 15 meters to match largest resolution of environmental data (30m resolution)
    st_buffer(dist = 15)%>%
    #join information on nest success to the filtered records in the nest point layer 
    inner_join(nest_fates,by="id")%>%
      distinct(.keep_all = TRUE)

  

  
## summarize environmental data at nest site points
# -----------------------------------------------------

  ## 1. Marsh Vegetation Class at each nest
  
  #Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
  zones<-c(1:8)
  #Create an empty list to hold a summary of vegetation in nests for each regional layer (will hold a list of summary tables, one for each region)
  out_list<-list()
  
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(vg_cls)) {
    #select zone layer and rename raster values
    layer<-dplyr::rename_with(vg_cls[[i]], function(x){x<-'value'},.cols = everything())
    
    #extract the raster cell value at each point (ID represents order of observations in points layer)
    out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
      mutate(region=zones[i])%>%
      #join additional attributes
      st_as_sf()%>%
      st_drop_geometry()%>%
      left_join(veg_class,by="value")%>%
      filter(!is.na(value))%>%
      dplyr::select(-value)
  }
  
  #empty region dataframes indicate no SALS nests in that region, remove them
  out_list2<-out_list[-c(2,6:8)]
  
  #combine vegetation values for nests across all regions into 1 dataframe
  veg_prop<-do.call("rbind",out_list2)

  
  
  
  ## 2. mean UVVR ratio at each nest
  uvvr_mean<-terra::extract(uvvr, vect(nests), bind=T)%>%
    st_as_sf()%>%
    st_drop_geometry()%>%
    dplyr::select(id,uvvr_mean)%>%
  #UVVR values greater than 2 shouldn't be trusted, set to NA
    mutate(uvvr_mean=ifelse(uvvr_mean>2,NA,uvvr_mean))
  
  
  ## 3. Change in UVVR at each nest
  uvvr_diff2<-terra::extract(uvvr_diff, vect(nests), bind=T)%>%
    st_as_sf()%>%
    st_drop_geometry()%>%
    dplyr::select(id,uvvr_diff)

  
  
  ## join mean UVVR and Proportion of vegetation classes in each nest buffer into 1 table
  final_dat_local<-left_join(veg_prop,uvvr_mean, by='id')%>%
    left_join(uvvr_diff2,by='id')%>%
    left_join(nest_fates,by='id')%>%
    dplyr::select(id,region,site_cd,Year,fate, veg_class,uvvr_mean,uvvr_diff)%>%
    rename(site=site_cd,yr=Year)
  
  write.csv(final_dat_local,paste0(dat_path,"nest_fate_model_local_dataset.csv"),row.names = F)
  
  

  
  


## summarize environmental data in each buffer zone (15m around each nest = 30m to match coarsest resolution of environmental data -UVVR)
# ------------------------------------------------------
## 1. Marsh Vegetation Classes
  
  #Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
  zones<-c(1:8)
  #Create an empty list to hold a summary of vegetation in nests for each regional layer (will hold a list of summary tables, one for each region)
  out_list<-list()

  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(vg_cls)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(vg_cls[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                         summarize_df=T,include_cols='id')%>%
  #convert cell count to area by multiplying count by cell area. 
  #then calculate the proportional area of each veg class
      group_by(id)%>%
      mutate(area=n*area,
           prop_area=round(area/sum(area),digits=2),
           region=zones[i])%>%
      ungroup()%>%
  #join additional attributes
    left_join(nests%>%st_drop_geometry(),by='id')%>%
    left_join(veg_class,by="value")%>%
    dplyr::select(-c("value","area","n"))%>%
    pivot_wider(names_from= "veg_class",values_from = "prop_area")
  }

  #empty region dataframes indicate no SALS nests in that region, remove them
  out_list2<-out_list[-c(2,6:8)]

  # align columns across dataframes 
  nms <- veg_class$veg_class   # Vector of columns you want in the data.frames (the names of each vegetation class)
  nms<-c(nms,"NA") #NA indicates outside of marsh area (no vegetation class within nest buffer)

  for(i in 1:length(out_list2)){
  
    Missing <- setdiff(nms, names(out_list2[[i]]))  # Find names of missing columns by comparing dataframe columns to vector of desired column names
   out_list2[[i]][Missing] <- 0                    # Add missing columns, fill observations with '0's (0% cover)
  
  }
  
  #combine vegetation values for nests across all regions into 1 dataframe
  veg_prop<-do.call("rbind",out_list2)
  #rename the "NA" vegetation variable column to "MISSING"
  colnames(veg_prop)[28]<-"MISSING"
  #remove nests outside the marsh area (with 100% of their nest buffer missing data)
  veg_prop2<-veg_prop%>%
    filter(MISSING!=1|is.na(MISSING))


  

## 2. Calculate mean UVVR within each nest buffer
  uvvr_mean<-exact_extract(uvvr, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                          summarize_df=T,include_cols='id')
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
  uvvr_mean<-group_by(uvvr_mean,id)%>%
    mutate(weighted=n*value)%>%
    summarise(uvvr_mean=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()%>%
    #UVVR values greater than 2 shouldn't be trusted, set to NA
    mutate(uvvr_mean=ifelse(uvvr_mean>2,NA,uvvr_mean))
  
  
  
  
## 3. Calculate UVVR change within each nest buffer
  uvvr_diff2<-exact_extract(uvvr_diff, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                           summarize_df=T,include_cols='id')
  #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
  uvvr_mean<-group_by(uvvr_diff,id)%>%
    mutate(weighted=n*value)%>%
    summarise(uvvr_diff=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()

  

## join UVVR and Proportion of vegetation classes in each nest buffer into 1 table
final_dat<-left_join(veg_prop2,uvvr_mean, by='id')%>%
  left_join(uvvr_diff2,by='id')
  dplyr::select(id,region,site_cd,Year,fate.y, HIMARSH,LOMARSH,POOL,PHRG,MISSING,STRM,MUD,UPLND,TERRBRD,uvvr_mean,uvvr_diff)%>%
  rename(fate=fate.y,site=site_cd,yr=Year)

write.csv(final_dat,paste0(dat_path,"nest_fate_model_buff15m_dataset.csv"),row.names = F)



#Other variables:
#NAIP raw data- create texture layer
#elevation


library(GLCMTextures)

ndvi3<-rast(paste0("C:/Users/mefen/Downloads/3_NDVI.tif"))
crs(ndvi3)

# 1. Raster Quanitization
ndvi3_rq<- quantize_raster(r = ndvi3, n_levels = 16, method = "equal range")
