
library(sf)
library(tidyverse)
library(tidyterra)
library(terra)#updated version of raster package
library(exactextractr)

#####################################################################################################################################
## Summarize environmental data in each buffer zone - 30m since this is our prediction resolution
######################################################################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

# Load all nest locations and environmental data
if(!exists("veg_class")){
  source("04a_Format_Load_Predictors.R")
}

speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")

for (j in 1:length(speciesnames)){
nests_buff<-st_read(paste0(path_out,"Intermediate_outputs/Nest_locations/",speciesnames[j],"_nest_pres_bg_buff.shp"))
  
#list for the 8 region outputs
out_list<-list()

#constants for area calculations
area_pix<-round(res(vg_cls[[1]])[1]^2)
area_buff<-as.numeric(st_area(nests_buff)[1])#30m diameter circle (15m buff)

## Summarize predictors in buffers
#---------------------------------------------------
#*exact_extract returns cells from the raster that are covered by the polygons. It won't return values for polygons that don't overlap the RS data.
#
## a. Marsh Vegetation Classes
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(vg_cls[[i]], st_transform(nests_buff,crs(vg_cls[[1]])), #set coord system of nest points to raster layer
                               function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
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
    dplyr::select(id,region,contains(veg_class$veg_class))%>% 
    #assign 0's to veg classes that are absent from nest buffers in marsh area
    replace(is.na(.), 0)
}

#empty region dataframes indicate no nests in that region

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

#combine vegetation values for nests across all regions into 1 dataframe
veg_prop<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(id,.keep_all=T)



## b. Mean UVVR 
#-------------------------------
uvvr_mean<-exact_extract(uvvr, st_transform(nests_buff,crs(uvvr)), #set coord system of nest points to raster layer
                         function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                         summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_mean2<-group_by(uvvr_mean,id)%>%
  #remove locations within raster layer extent but not in marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
  mutate(n=n/sum(n,na.rm=T),
         weighted=n*value)%>%
  summarise(uvvr_mean=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()


## c. UVVR change 
#-------------------------------
uvvr_diff2<-exact_extract(uvvr_diff, st_transform(nests_buff,crs(uvvr_diff)), #set coord system of nest points to raster layer
                          function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                          summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_diff2<-group_by(uvvr_diff2,id)%>%
  #remove locations within raster layer extent but not in marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
  mutate(n=n/sum(n,na.rm=T),
         weighted=n*value)%>%
  summarise(uvvr_diff=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()



## j. Tidal Restrictions
#--------------------------------
tideres2<-exact_extract(tideres, st_transform(nests_buff,crs(tideres)),
                        function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                        summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
tideres2<-group_by(tideres2,id)%>%
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
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(ndvi[[i]], st_transform(nests_buff,crs(ndvi[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
      group_by(id)%>%
      #remove locations within raster layer extent but not in this region's marsh area.
      filter(!(is.nan(value)|is.na(value)))%>%
      #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each nest buffer(id)
      mutate(n2=n/sum(n,na.rm=T),
           weighted=n2*value)%>%
      summarise(ndvi=round(sum(weighted,na.rm=T),digits=5))%>%
      ungroup()
}

#empty region dataframes indicate no nests in that region
#combine values for nests across all regions into 1 dataframe
ndvi_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(id,.keep_all=T)


## e. PCA
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(pca)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(pca[[i]], st_transform(nests_buff,crs(pca[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
  #remove locations within raster layer extent but not in this region's marsh area.
  filter(!(is.nan(value)|is.na(value)))%>%
    mutate(n=n/sum(n,na.rm=T),
           weighted=n*value)%>%
    summarise(pca=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no nests in that region
#combine values for nests across all regions into 1 dataframe
pca_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(id,.keep_all=T)



## f. Elevation
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(dem)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(dem[[i]], st_transform(nests_buff,crs(dem[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    #remove locations within raster layer extent but not in this region's marsh area.
    filter(!(is.nan(value)|is.na(value)))%>%
    mutate(n=n/sum(n,na.rm=T),
           weighted=n*value)%>%
    summarise(elevation=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no nests in that region
#combine values for nests across all regions into 1 dataframe
dem_buff<-do.call("rbind",out_list)%>%
  # some points overlap marsh area in two regions (rasters overlap). Just pick the values from one region.
  distinct(id,.keep_all=T)


## g. join UVVR and Proportion of vegetation classes in each nest buffer into 1 table
#-------------------------------
# joing back to the original nest location table. Any locations without information in the raster extraction tables
  # should come back as NA- these are true missing values (do not replace with 0)
final_dat<-left_join(nests_buff,veg_prop,by='id')%>%
  left_join(uvvr_mean2, by='id')%>%
  left_join(uvvr_diff2,by='id')%>%
  left_join(ndvi_buff,by='id')%>%
  left_join(pca_buff,by='id')%>%
  left_join(dem_buff,by='id')%>%
  left_join(tideres2,by='id')%>%
  dplyr::select(id,bp,region,site,Year,fate, 
                HIMARSH,LOMARSH,POOL,PHRG,STRM,MUD,UPLND,TERRBRD,
                uvvr_mean,uvvr_diff, tideres,
                ndvi,pca,elevation)%>%
  st_drop_geometry()

#save(uvvr_mean,uvvr_diff2,ndvi_buff,pca_buff,entro_buff,corr_buff, file=paste0(path_out,"Final_outputs/4c_final_files.rds"))
write.csv(final_dat,paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[j],"_nest_vars_buff15.csv"),row.names = F)

}
