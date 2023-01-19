library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(tmap)
library(exactextractr)
library(GLCMTextures)
#run test, get full zone 1 NDVI into cluster as well as UVVR data, finish loading pca

# test data
#uvvr<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=60)%>%
#  dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything())
#values(uvvr)<-runif(ncell(uvvr),0,2)


#ndvi<-list()
#ndvi[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
#values(ndvi[[1]])<-runif(ncell(ndvi[[1]]),0,1)
#ndvi[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
#values(ndvi[[2]])<-runif(ncell(ndvi[[2]]),0,1)

#pca<-list()
#pca[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
#values(pca[[1]])<-runif(ncell(pca[[1]]),0,1)
#pca[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
#values(pca[[2]])<-runif(ncell(pca[[2]]),0,1)

#t_veg<-list()
#t_veg[[1]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=40)
#values(t_veg[[1]])<-round(runif(ncell(t_veg[[1]]),1,9))
#t_veg[[2]]<-rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=41, ymax=60)
#values(t_veg[[2]])<-round(runif(ncell(t_veg[[2]]),1,9))
#vg_cls<-t_veg

## Set file path to data
# -------------------------------------------
#dat_path<-"D:/Research/SHARP/Data/"
dat_path<-"/home/FCAM/mlfeng/Data/"



## 1. Format Nest Observations
# -------------------------------------
# load nest observation shapefile
nests<-st_read(paste0(dat_path,"Nest_Locations/nest_locations_01_3_23.shp"))%>%
  # filter records to just SALS or species of interest
  filter(Species=="SALS"&
  # also filter records missing coordinate information or that have coordinate errors
           crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  # convert all coordinates to degrees (NAD83) and create Long and Lat columns
  st_transform("EPSG:4269")%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
           #create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                         fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(name=id,latitude=Lat,longitude=Long,fate,Year,site=site_cd)

  # output table file
#write.csv(st_drop_geometry(nests),paste0(dat_path,"Nest_Coords_fates_SALS_01_5_23.csv"),row.names=F)
nests<-rename(nests,id=name)

nests_buff<- nests%>%
  #buffer nest points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)


## 2. Sample environmental data at nests using their original resolutions
# --------------------------------------------------------------------------

## 2. a. Load in data
#---------------------


## Environmental Predictor 1: Marsh Vegetation Classes

#Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
zones<-c(1:8)

#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

#read as raster layers
vg_cls<-map(file_list,rast)

#set coordinate systems
  # Set all coordinate systems to ESPG:4269 (NAD 1983 geographic coord system)
#for(i in 1:length(vg_cls)){
#vg_cls[[i]] <- terra::project(vg_cls[[i]],"EPSG:4269")
#}

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
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything()) 
# values above 2 are not accurate, set to NA
uvvr[uvvr>2]<-NA

# B. UVVR change (2014-2018)
uvvr14<-rast(paste0(dat_path,"UVVR/uvvr_14_utm18_mean_ext.tif"))%>%
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_14'},.cols = everything())
# values above 2 are not accurate, set to NA
uvvr14[uvvr14>2]<-NA
uvvr18<-rast(paste0(dat_path,"UVVR/uvvr_18_utm18_mean_ext.tif"))%>%
  #terra::project("EPSG:4269")%>%
  dplyr::rename_with(function(x){x<-'uvvr_18'},.cols = everything())
uvvr18[uvvr18>2]<-NA
# subtract 2018 from 2014 to get change
uvvr_diff<-uvvr18-uvvr14
uvvr_diff<-uvvr_diff%>%
  dplyr::rename_with(function(x){x<-'uvvr_diff'},.cols = everything())



## Enviromental Predictor 3: NAIP

#list raster files
file_list1<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "NDVI.TIF$",full.names=T)))
file_list2<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "PCA.tif$",full.names=T)))

#read as raster layers
ndvi<-map(file_list1,rast)
pca<-map(file_list2,rast)

#set values below 0 to 0 - barren land or water
for(i in 1:length(ndvi)){
  dat<-ndvi[[i]]
  dat[dat<0]<-0
  ndvi[[i]]<-dat
}


## 3. Calculate texture surfaces
# --------------------------------------------

## 3. a. Raster Quanitization 
#----------------------------
  # use equal range method when comparing across several datasets (can set gloabl min and max range using max_val and min_val) - splits data into equal ranges
  # alternative option is "equal prob" which splits by quantiles, used in original paper (Haralick and Shanmugam 1973)
  # NDVI of live plants ranges 0:1, set range and breaks at 0.1 intervals
txt_homo<-txt_entro<-txt_corr<-list()
for (i in 1:length(ndvi)){
  ndvi_rq<- quantize_raster(r = ndvi[[i]], n_levels = 10, method = "equal range")
  
## 3. b. calc texture metric surfaces
#----------------------------
  # use window size of 3x3m 
  # calculate across multiple shifts (invariant)
  txt<- glcm_textures(ndvi_rq, w = c(3,3), n_levels =10, 
                      metrics = c("glcm_homogeneity", 
                                  "glcm_entropy", "glcm_correlation"),
                      quantization = "none", shift = list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)))
  
  txt_homo[[i]]<-txt[[1]]
  txt_entro[[i]]<-txt[[2]]
  txt_corr[[i]]<-txt[[3]]
  
}



## 4. summarize at nest points
# ---------------------------------------------
#test points
# grab 20 cell index numbers at random
#samp <- round(runif(20,1,ncell(uvvr)))

# their location coordinates
#samplocs <- xyFromCell(uvvr, samp)

# convert to a data frame
#samp <- as.data.frame(cbind(samplocs, samp))
#names(samp) <- c('x', 'y', 'id')
#samp<-mutate(samp,
#             Year=round(runif(20,2011,2020)),
#             site="AT",
#             fate=1)
# spatial object
#nests <- st_as_sf(as.data.frame(samp), coords = c('x', 'y'), crs = crs(uvvr))
#nests_buff<-nests%>%
#  st_buffer(dist=15)

## 4. a. Marsh Vegetation Classes
#-------------------------------

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
    dplyr::select(-value,-Year,-site,-fate)
}

#empty region dataframes indicate no SALS nests in that region

#combine vegetation values for nests across all regions into 1 dataframe
veg_prop<-do.call("rbind",out_list)



## 4. b. mean UVVR ratio at each nest
#-------------------------------
uvvr_mean<-terra::extract(uvvr, vect(nests), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_mean)


## 4. c. Change in UVVR at each nest
#-------------------------------
uvvr_diff2<-terra::extract(uvvr_diff, vect(nests), bind=T)%>%
  st_as_sf()%>%
  st_drop_geometry()%>%
  dplyr::select(id,uvvr_diff)#=uvvr_mean)#For test



## 4. d. NDVI 
#-------------------------------

for(i in 1:length(ndvi)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(ndvi[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,ndvi=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
ndvi2<-do.call("rbind",out_list)



## 4. e. PCA
#-------------------------------

for(i in 1:length(pca)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(pca[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,pca=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
pca2<-do.call("rbind",out_list)



## 4. f. Homogeneity (of NDVI)
#-------------------------------

for(i in 1:length(txt_homo)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_homo[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,hom_txt=value)
}

#empty region dataframes indicate no SALS nests in that region

#combine NDVI values for nests across all regions into 1 dataframe
txt_homo2<-do.call("rbind",out_list)


## 4. g. Entropy (of NDVI)
#-------------------------------

for(i in 1:length(txt_entro)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_entro[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,ent_txt=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
txt_entro2<-do.call("rbind",out_list)




## 4. h. Correlation (of NDVI)
#-------------------------------

for(i in 1:length(txt_corr)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_corr[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(id,cor_txt=value)
}

#empty region dataframes indicate no SALS nests in that region
#combine NDVI values for nests across all regions into 1 dataframe
txt_corr2<-do.call("rbind",out_list)



## 4. i. join mean UVVR and Proportion of vegetation classes in each nest buffer into 1 table
#-------------------------------
final_dat_local<-left_join(ndvi2,txt_homo2, by='id')%>%
  left_join(txt_entro2,by='id')%>%
  left_join(txt_corr2,by='id')%>%
  dplyr::select(id,ndvi,hom_txt,ent_txt,cor_txt)%>%
  distinct(.keep_all = T)

write.csv(final_dat_local,paste0(dat_path,"SALS_nest_vars_local.csv"),row.names = F)





## 5. summarize environmental data in each buffer zone (15m around each nest = 30m to match coarsest resolution of environmental data -UVVR)
# -----------------------------------------------------------------------------------------------------------------

## 5.a. Marsh Vegetation Classes
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(vg_cls[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell count to area by multiplying count by cell area. 
    #then calculate the proportional area of each veg class
    group_by(id)%>%
    distinct(.keep_all = T)%>%
    mutate(area=n*area,
           prop_area=round(area/sum(area),digits=2),
           region=zones[i])%>%
    ungroup()%>%
    #join additional attributes
    left_join(nests%>%st_drop_geometry(),by='id')%>%
    left_join(veg_class,by="value")%>%
    dplyr::select(-c("value","area","n"))%>%
    distinct(.keep_all=T)%>%
    pivot_wider(names_from= "veg_class",values_from = "prop_area")
}

#empty region dataframes indicate no SALS nests in that region
# align columns across dataframes 
nms <- veg_class$veg_class   # Vector of columns you want in the data.frames (the names of each vegetation class)
nms<-c(nms,"NA") #NA indicates outside of marsh area (no vegetation class within nest buffer)

for(i in 1:length(out_list)){
  
  Missing <- setdiff(nms, names(out_list[[i]]))  # Find names of missing columns by comparing dataframe columns to vector of desired column names
  out_list[[i]][Missing] <- 0                    # Add missing columns, fill observations with '0's (0% cover)
  
}

#combine vegetation values for nests across all regions into 1 dataframe
veg_prop<-do.call("rbind",out_list)
#rename the "NA" vegetation variable column (last column) to "MISSING"
colnames(veg_prop)[ncol(veg_prop)]<-"MISSING"





## 5.b. Mean UVVR 
#-------------------------------
uvvr_mean<-exact_extract(uvvr, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                         summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_mean<-group_by(uvvr_mean,id)%>%
  mutate(weighted=n*value)%>%
  summarise(uvvr_mean=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()


## 5.c. UVVR change 
#-------------------------------
uvvr_diff2<-exact_extract(uvvr_diff, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                          summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_diff2<-group_by(uvvr_diff2,id)%>%
  mutate(weighted=n*value)%>%
  summarise(uvvr_diff=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()



## 5.d. NDVI
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(ndvi[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(ndvi=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
ndvi_buff<-do.call("rbind",out_list)


## 5.e. PCA
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(pca)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(pca[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(pca=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
pca_buff<-do.call("rbind",out_list)




## 5.f. Homogeneity
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_homo[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(hom_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
homo_buff<-do.call("rbind",out_list)



## 5.g. Entropy
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_entro[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(ent_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
entro_buff<-do.call("rbind",out_list)



## 5.h. Correlation
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_corr[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(cor_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
corr_buff<-do.call("rbind",out_list)




## 5.i. join UVVR and Proportion of vegetation classes in each nest buffer into 1 table
#-------------------------------
final_dat<-left_join(veg_prop,uvvr_mean, by='id')%>%
  left_join(uvvr_diff2,by='id')%>%
  left_join(ndvi_buff,by='id')%>%
  left_join(pca_buff,by='id')%>%
  left_join(homo_buff,by='id')%>%
  left_join(entro_buff,by='id')%>%
  left_join(corr_buff,by='id')%>%
dplyr::select(id,region,site,Year,fate, 
              HIMARSH,LOMARSH,POOL,PHRG,MISSING,STRM,MUD,UPLND,TERRBRD,
              uvvr_mean,uvvr_diff,
              ndvi,pca,hom_txt,ent_txt,cor_txt)%>%
  distinct(.keep_all = T)

write.csv(final_dat,paste0(dat_path,"SALS_nest_vars_buff15.csv"),row.names = F)





## 6. Create prediction surface
# ---------------------------------------------------------------------------------------------
# align extent and resolution of all environmental rasters
# use a prediction surface resolution of 3 meters (matches the Correll data, plus middle ground between NAIP and UVVR)
# Set all coordinate systems to ESPG:4269 (NAD 1983 geographic coord system)
# Use Atlantic coast marsh extent (refer to Correll marsh layer) for each zone (n = 8)- predict by zone to reduce computation time



## 6.a. Create an identical list of templates with cell ID numbers based on marsh veg zone layers
# -------------------
cellid<-list()
for(i in 1:length(vg_cls)){
  dat<-vg_cls[[i]]
  values(dat)<-as.numeric(paste0(i,1:ncell(dat)))
  cellid[[i]]<-dat
}


## 6.b. Disaggregate UVVR data - 30m to 3m
# --------------------

  # UVVR mean
    # create separate files for each zone
uvvr_temp<-list()
for (i in 1:length(vg_cls)) {
uvvr_temp[[i]]<-crop(extend(uvvr,vg_cls[[i]]),vg_cls[[i]])%>%
           terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
writeRaster(uvvr_temp[[i]],filename=paste0(dat_path,"UVVR/UVVR_annual_mean/Z",i,"_uvvr_mean_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT8U'))

}

  # UVVR change
    # create separate files for each zone
uvvr_diff_temp<-list()
for (i in 1:length(vg_cls)) {
  uvvr_diff_temp[[i]]<-crop(extend(uvvr_diff,vg_cls[[i]]),vg_cls[[i]])%>%
    resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(uvvr_diff_temp[[i]],filename=paste0(dat_path,"UVVR/UVVR_change_14_18/Z",i,"_uvvr_diff_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT8S'))
  
}




## 6.c. Aggregate NAIP data - 1m to 3m
# ----------------------
  # NDVI
    # create separate files for each zone
ndvi_temp<-list()
for(i in 1:length(ndvi)) {
  ndvi_temp[[i]]<-crop(extend(ndvi[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(ndvi_temp[[i]],filename=paste0(path_out,"Z",i,"_NDVI_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


  # PCA
    # create separate files for each zone
pca_temp<-list()
for(i in 1:length(pca)) {
  pca_temp[[i]]<-crop(extend(pca[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(pca_temp[[i]],filename=paste0(path_out,"Z",i,"_PCA_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


  # Homogeneity
    # create separate files for each zone
homo_temp<-list()
for(i in 1:length(txt_homo)) {
  homo_temp[[i]]<-crop(extend(txt_homo[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(homo_temp[[i]],filename=paste0(path_out,"Z",i,"_homo_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


  # Entropy
    # create separate files for each zone
entro_temp<-list()
for(i in 1:length(txt_entro)) {
  entro_temp[[i]]<-crop(extend(txt_entro[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(entro_temp[[i]], filename=paste0(path_out,"Z",i,"_entro_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}



  # Correlation
    # create separate files for each zone
corr_temp<-list()
for(i in 1:length(txt_corr)) {
  corr_temp[[i]]<-crop(extend(txt_corr[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
    #write to file
  writeRaster(corr_temp[[i]], filename=paste0(path_out,"Z",i,"_corr_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}




## 6.d. combine all values on prediction surface into a single dataframe by cell ID
# --------------------
dat_all<-list()
mat<-list()
for(i in 1:length(vg_cls)){
dat_all[[i]]<-c(cellid[[i]],vg_cls[[i]],uvvr_temp[[i]],uvvr_diff_temp[[i]],
                ndvi_temp[[i]],pca_temp[[i]],homo_temp[[i]],entro_temp[[i]],corr_temp[[i]])
names(dat_all[[i]])<-c("id","vg_clss","uvvr_mean","uvvr_diff","ndvi","pca","hom_txt","ent_txt","cor_txt")
mat[[i]]<-as.data.frame(as.matrix(dat_all[[i]]))
}
mat_final<-do.call("rbind",mat)
