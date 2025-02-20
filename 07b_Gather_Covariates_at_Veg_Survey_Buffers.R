library(tidyverse)
library(terra)
library(sf)
library(tidyterra)
library(exactextractr)



### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


# load vegetation survey data
veg_data<-read.csv(paste0(path_out,"Intermediate_outputs/Survey_Vegetation/processed_rapid_veg.csv"))%>%
  #need a unique identifer since id's are shared across multiple visit years
  rownames_to_column()
  

# Load all remote sensing data
source("04a_Format_Load_Predictors.R")


### Gather covariates at veg points
#----------------------------------------------------------------------------
# Latitude
veg_buff<-st_as_sf(veg_data,coords=c("Long","Lat"),crs="EPSG:4269")%>%
  st_transform("EPSG:26918")%>%
  cbind(veg_data%>%dplyr::select(Lat))%>%
  st_buffer(15)

#Year, Day of year
veg_buff<-veg_buff%>%
  rename(Year=year)
veg_buff<-veg_buff%>%
  mutate(doy=yday(ymd(Date)))

# Time since tide (set this to the time that maximizes nesting from the models) (assume nests that nest immediate after the tide are most successful)
veg_buff<-veg_buff%>%
  mutate(time_since_tide=1)

# Remote sensing covariates
  ## Summarize predictors in buffers
  #*Exact extract only returns values for points that intersect the raster
  ## a. Marsh Vegetation Classes (percent area)
  #-------------------------------
#list for the 8 region outputs
out_list<-list()

#constants for area calculations
area_pix<-round(res(vg_cls[[1]])[1]^2)
area_buff<-as.numeric(st_area(veg_buff)[1])#30m diameter circle (15m buff)

#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(vg_cls[[i]], st_transform(veg_buff,crs(vg_cls[[1]])), #set coord system of nest points to raster layer
                               function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='rowname')%>%
    #convert cell count to area by multiplying count by cell area. 
    #then calculate the proportional area of each veg class
    mutate(area=n*area_pix,
           prop_area=round(area/area_buff,digits=3),
           #this will generate accurate region assignments for locations (use these)
           region=zones[i])%>%
    ungroup()%>%
    #join additional attributes
    left_join(veg_class,by="value")%>%
    dplyr::select(-c("value","area","n"))%>%
    pivot_wider(names_from= "veg_class",values_from = "prop_area")%>%
    # if a location only has "NA" as a cover class, remove. Those locations are outside the marsh area for this region.
    filter(!if_all(contains(veg_class$veg_class),is.na))%>%
    # if the table has an NA column, remove
    dplyr::select(rowname,region,contains(veg_class$veg_class))%>% 
    #assign 0's to veg classes that are absent from nest buffers in marsh area
    replace(is.na(.), 0)
}

#empty region dataframes indicate no plots in that region

# Make sure locations in all regions have all veg classes sampled 
# Vector of columns you want in the data.frames (the names of each vegetation class)
nms <- veg_class$veg_class   
# if a DF is missing a veg class, it didn't appear at any locations in that region. Add the class and fill with zero % cover.
for(i in 1:length(out_list)){
  # Find names of missing columns by comparing dataframe columns to vector of desired column names  
  Missing <- setdiff(nms, names(out_list[[i]]))  
  # Add missing columns, fill observations with '0's (0% cover)
  out_list[[i]][Missing] <- 0                    
  
}

#combine vegetation values for plots across all regions into 1 dataframe
veg_prop<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(rowname,.keep_all=T)



  ## b. Mean UVVR 
  #-------------------------------
uvvr_mean<-exact_extract(uvvr, st_transform(veg_buff,crs(uvvr)), #set coord system of veg points to raster layer
                         function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                         summarize_df=T,include_cols='rowname')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each plot buffer(id)
uvvr_mean2<-group_by(uvvr_mean,rowname)%>%
  #remove locations within raster layer extent but not in marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
  mutate(n=n/sum(n,na.rm=T),
         weighted=n*value)%>%
  summarise(uvvr_mean=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()


## c. UVVR change 
#-------------------------------
uvvr_diff2<-exact_extract(uvvr_diff, st_transform(veg_buff,crs(uvvr_diff)), #set coord system of veg points to raster layer
                          function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                          summarize_df=T,include_cols='rowname')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each plot buffer(id)
uvvr_diff2<-group_by(uvvr_diff2,rowname)%>%
  #remove locations within raster layer extent but not in marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
  mutate(n=n/sum(n,na.rm=T),
         weighted=n*value)%>%
  summarise(uvvr_diff=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()



## j. Tidal Restrictions
#--------------------------------
tideres2<-exact_extract(tideres, st_transform(veg_buff,crs(tideres)),
                        function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                        summarize_df=T,include_cols='rowname')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each plot buffer(id)
tideres2<-group_by(tideres2,rowname)%>%
  #remove locations within raster layer extent but not in marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
  mutate(n=n/sum(n,na.rm=T),
         weighted=n*value)%>%
  summarise(tideres=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()


## d. NDVI
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each veg buffer (coverage fraction)
  out_list[[i]]<-exact_extract(ndvi[[i]], st_transform(veg_buff,crs(ndvi[[1]])), function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='rowname')%>%
    group_by(rowname)%>%
    #remove locations within raster layer extent but not in this region's marsh area.
    filter(!(is.nan(value)|is.na(value)))%>%
    #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each veg buffer(rowname)
    mutate(n2=n/sum(n,na.rm=T),
           weighted=n2*value)%>%
    summarise(ndvi=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no plots in that region
#combine values for nests across all regions into 1 dataframe
ndvi_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(rowname,.keep_all=T)


## e. PCA
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(pca)) {
  #summarize the number of cells weighted by proportion of coverage within each veg buffer (coverage fraction)
  out_list[[i]]<-exact_extract(pca[[i]], st_transform(veg_buff,crs(pca[[1]])), function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='rowname')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each veg buffer(id)
    group_by(rowname)%>%
    #remove locations within raster layer extent but not in this region's marsh area.
    filter(!(is.nan(value)|is.na(value)))%>%
    mutate(n=n/sum(n,na.rm=T),
           weighted=n*value)%>%
    summarise(pca=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no plots in that region
#combine values for nests across all regions into 1 dataframe
pca_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(rowname,.keep_all=T)



## f. Elevation
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(dem)) {
  #summarize the number of cells weighted by proportion of coverage within each veg buffer (coverage fraction)
  out_list[[i]]<-exact_extract(dem[[i]], st_transform(veg_buff,crs(dem[[1]])), function(df) summarize(group_by(df,value,rowname),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='rowname')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each veg buffer(id)
    group_by(rowname)%>%
    #remove locations within raster layer extent but not in this region's marsh area.
    filter(!(is.nan(value)|is.na(value)))%>%
    mutate(n=n/sum(n,na.rm=T),
           weighted=n*value)%>%
    summarise(elevation=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no plots in that region
#combine values for nests across all regions into 1 dataframe
dem_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(rowname,.keep_all=T)


## g. join UVVR and Proportion of vegetation classes in each nest buffer into 1 table
#-------------------------------
# joing back to the original nest location table. Any locations without information in the raster extraction tables
# should come back as NA- these are true missing values (do not replace with 0)
final_dat<-left_join(veg_buff,veg_prop,by='rowname')%>%
  left_join(uvvr_mean2, by='rowname')%>%
  left_join(uvvr_diff2,by='rowname')%>%
  left_join(ndvi_buff,by='rowname')%>%
  left_join(pca_buff,by='rowname')%>%
  left_join(dem_buff,by='rowname')%>%
  left_join(tideres2,by='rowname')%>%
  distinct(rowname,.keep_all=T)%>%
  st_drop_geometry()

  
  write.csv(final_dat,paste0(path_out,"Intermediate_outputs/Survey_Vegetation/veg_buff_RS_covariates_15buff.csv"),row.names = F)
  