


###*** first test with z 4 and 5 on local
###*then upload all NDVI and PCA files to cluster and run there
library(tidyverse)
library(GLCMTextures)
library(terra)
library(tidyterra)

## Set file path to data
# -------------------------------------------
dat_path<-"G:/My Drive/Research/SHARP/Data/"
#dat_path<-"/home/FCAM/mlfeng/environmental_dat/"

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

plot(ndvi[[1]])
## Calculate texture surfaces
# ------------------------------
# 1. Raster Quanitization 
# use equal range method when comparing across several datasets (can set gloabl min and max range using max_val and min_val) - splits data into equal ranges
# alternative option is "equal prob" which splits by quantiles, used in original paper (Haralick and Shanmugam 1973)
# NDVI of live plants ranges 0:1, set range and breaks at 0.1 intervals
txt_homo<-txt_entro<-txt_corr<-list()
for (i in 1:length(ndvi)){
ndvi_rq<- quantize_raster(r = ndvi[[1]], n_levels = 10, method = "equal range")

# 2. calc texture metric surfaces
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



## summarize NAIP data at nest site points
# -----------------------------------------------------

#load point layer of spatially projected nest records
nests<-st_read(paste0(dat_path,"Nest_Locations/nest_locations_01_3_23.shp"))%>%
  #filter records to just SALS
  filter(Species=="SALS"&crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  #remove unneeded variables
  select(-c("Lat","Long","crd_typ","mssng_l_","mssng_s_","mssng_c","fate"))%>%
  #project to UTM 18
  st_transform("EPSG:26918")%>%
  distinct(.keep_all = TRUE)

nests_buff<- nests%>%
  #buffer nest points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)


# 8 distinct geographic regions along the Eastern coast
zones<-c(1:8)
#Create an empty list to hold a summary of NAIP in nests for each regional layer (will hold a list of summary tables, one for each region)
out_list<-list()

#for each regional zone (1-8) along the east coast...

## 1. NDVI 

for(i in 1:length(ndvi)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(ndvi[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(-value)
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list
#combine NDVI values for nests across all regions into 1 dataframe
ndvi2<-do.call("rbind",out_list2)



## 2. PCA

for(i in 1:length(pca)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(pca[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(-value)
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list
#combine NDVI values for nests across all regions into 1 dataframe
pca2<-do.call("rbind",out_list2)



## 3. Homogeneity (of NDVI)

for(i in 1:length(txt_homo)) {
#select zone layer and rename raster values
layer<-dplyr::rename_with(txt_homo[[i]], function(x){x<-'value'},.cols = everything())

#extract the raster cell value at each point (ID represents order of observations in points layer)
out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
  #join additional attributes
  st_as_sf()%>%
  st_drop_geometry()%>%
  filter(!is.na(value))%>%
  dplyr::select(-value)
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine NDVI values for nests across all regions into 1 dataframe
txt_homo2<-do.call("rbind",out_list2)


## 4. Entropy (of NDVI)

for(i in 1:length(txt_entro)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_entro[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(-value)
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine NDVI values for nests across all regions into 1 dataframe
txt_entro2<-do.call("rbind",out_list2)




## 5. Correlation (of NDVI)

for(i in 1:length(txt_corr)) {
  #select zone layer and rename raster values
  layer<-dplyr::rename_with(txt_corr[[i]], function(x){x<-'value'},.cols = everything())
  
  #extract the raster cell value at each point (ID represents order of observations in points layer)
  out_list[[i]]<-terra::extract(layer, vect(nests), bind = T)%>%
    #join additional attributes
    st_as_sf()%>%
    st_drop_geometry()%>%
    filter(!is.na(value))%>%
    dplyr::select(-value)
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine NDVI values for nests across all regions into 1 dataframe
txt_corr2<-do.call("rbind",out_list2)



## join mean UVVR and Proportion of vegetation classes in each nest buffer into 1 table
final_dat_local<-left_join(ndvi2,txt_homo2, by='id')%>%
  left_join(txt_entro2,by='id')%>%
  left_join(txt_corr2,by='id')%>%
  dplyr::select(id,ndvi,txt_homo,txt_entro,txt_corr)

write.csv(final_dat_local,paste0(dat_path,"SALSnests_NAIP_vars_local.csv"),row.names = F)




## summarize environmental data in each buffer zone (15m)
# ------------------------------------------------------
## 1. NDVI

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

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine values for nests across all regions into 1 dataframe
ndvi_buff<-do.call("rbind",out_list2)


## 2. PCA

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

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine values for nests across all regions into 1 dataframe
pca_buff<-do.call("rbind",out_list2)




## 3. Homogeneity
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_homo[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(homo=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine values for nests across all regions into 1 dataframe
homo_buff<-do.call("rbind",out_list2)



## 4. Entropy
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_entro[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(entro=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine values for nests across all regions into 1 dataframe
entro_buff<-do.call("rbind",out_list2)



## 5. Correlation
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_corr[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(corr=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region, remove them
#out_list2<-out_list[-c(2,6:8)]
out_list2<-out_list

#combine values for nests across all regions into 1 dataframe
corr_buff<-do.call("rbind",out_list2)





## join NAIP variables in each nest buffer into 1 table
final_dat<-left_join(ndvi_buff,homo_buff, by='id')%>%
  left_join(entro_buff,by='id')%>%
  left_join(corr_buff,by='id')%>%
dplyr::select(id,ndvi,homo,entro,corr)

write.csv(final_dat,paste0(dat_path,"SALSnests_NAIP_vars_15mbuff.csv"),row.names = F)







## Aggregate 1 m NAIP to 3 m resolution aligned with Correll Veg layer
# ----------------------------------------------------------------------

# list marsh extent raster files for each zone (use as a template)
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))
# read as raster layers
vg_cls<-map(file_list,rast)



for(i in 1:length(ndvi)) {
  ndvi_agg<-crop(extend(ndvi[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
  #write to file
  writeRaster(ndvi_agg,filename=paste0(dat_path,"Z",i,"_NDVI_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}

for(i in 1:length(pca)) {
  pca_agg<-crop(extend(pca[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
  #write to file
  writeRaster(pca_agg,filename=paste0(dat_path,"Z",i,"_NPCA_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


for(i in 1:length(txt_homo)) {
  homo_agg<-crop(extend(txt_homo[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
  #write to file
  writeRaster(homo_agg,filename=paste0(dat_path,"Z",i,"_homo_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


for(i in 1:length(txt_entro)) {
  entro_agg<-crop(extend(txt_entro[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
  #write to file
  writeRaster(entro_agg, filename=paste0(dat_path,"Z",i,"_entro_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}


for(i in 1:length(txt_corr)) {
  corr_agg<-crop(extend(txt_corr[[i]],vg_cls[[i]]),vg_cls[[i]])%>%
    terra::resample(vg_cls[[i]],method="bilinear")
  #write to file
  writeRaster(corr_agg, filename=paste0(dat_path,"Z",i,"_corr_3m.tif"), overwrite=TRUE, wopt= list(gdal=c("COMPRESS=NONE"), datatype='FLT4S'))
}



