library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(tmap)
library(exactextractr)
library(GLCMTextures)

## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"



## 1. Format Nest Observations
# -------------------------------------
# load nest observation shapefile

nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  # filter records to just SALS or species of interest
  filter(Species=="SALS"&
           # also filter records missing coordinate information or that have coordinate errors
           crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  # convert all coordinates to degrees (NAD83) and create Long and Lat columns
  st_transform("EPSG:4269")%>% #"EPSG:26918"
  mutate(bp="p", #mark as a nest presence location
         Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         #create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                        fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(name=id,latitude=Lat,longitude=Long,fate,Year,site=site_cd,bp)%>%
  distinct(.keep_all = TRUE)

# output csv file for all nest coordinates and their fate, if recorded
if(!file.exists(paste0(path_out,"Final_outputs/Nest_Coords_fates_SALS_01_5_23.csv"))){
write.csv(st_drop_geometry(nests),paste0(path_out,"Final_outputs/Nest_Coords_fates_SALS_01_5_23.csv"),row.names=F)
}


# nest points sf object
nests<-rename(nests,id=name)

#nest buffers sf object
nests_buff<- nests%>%
  #buffer nest points by 5 meters 
  st_buffer(dist = 5)%>%
  distinct(.keep_all = TRUE)






## 2. Add random veg and background points to nest data
# --------------------------------------------------------------------------
#load veg and background points
source("03_background_selection.R")

#filter random background points that are within 5 meters of a nest
bg_points2<-bg_points[!(st_intersects(bg_points, nests_buff) %>% lengths > 0),] #https://stackoverflow.com/questions/57014381/how-to-filter-an-r-simple-features-collection-using-sf-methods-like-st-intersect

#make sure columns align
names(nests)
names(veg)
names(bg_points)

#add points to nest data
nests<-rbind(nests,bg_points2,veg)%>%
  distinct(.keep_all = TRUE)


#nest and background buffers sf object
nests_buff<- nests%>%
  #buffer nest and background points by 15 meters to match largest resolution of environmental data (30m resolution)
  st_buffer(dist = 15)%>%
  distinct(.keep_all = TRUE)







## 3. Load environmental data at their original resolutions
# --------------------------------------------------------------------------


## Environmental Predictor 1: Marsh Vegetation Classes

#Marsh vegetation data comes in individual rasters for 8 distinct geographic regions along the Eastern coast
zones<-c(1:8)

#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

#read as raster layers
vg_cls<-map(file_list,rast)

#vegetation class definitions
veg_class<-data.frame(value=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))




## Environmental Predictor 2: Unvegetated-Vegetated ratio
# transformed uvvr coordinate system from WGS to NAD83 in ArcPro
# just want band 4 (UVVR). Band 1 is unveg, 2 is veg, and 3 is proportion of water.

# Read raster if it already exists. Otherwise process the UVVR layer.
if(file.exists(paste0(path_out,"Final_outputs/UVVR/uvvr_mean_noOutlier.tif"))){
  uvvr<-  rast(paste0(path_out,"Final_outputs/UVVR/uvvr_mean_noOutlier.tif"))
  uvvr_diff<-  rast(paste0(path_out,"Final_outputs/UVVR/uvvr_diff_noOutlier.tif"))
}
if(!file.exists(paste0(path_out,"Final_outputs/UVVR/uvvr_mean_noOutlier.tif"))){
  
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

  
  writeRaster(uvvr,filename = paste0(path_out,"Final_outputs/UVVR/uvvr_mean_noOutlier.tif"),overwrite=T)
  writeRaster(uvvr_diff,filename = paste0(path_out,"Final_outputs/UVVR/uvvr_diff_noOutlier.tif"),overwrite=T)
  
}


## Environmental Predictor 3: NAIP
#if files do not exist
if(!(file.exists(paste0(path_out,"Final_outputs/Correll_NAIP/Z1_zeroed_NDVI.tif")))){
#list raster files
file_list1<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))
file_list2<-unlist(map(paste0(path_out,"Final_outputs/Correll_NAIP/"),~list.files(.,pattern = "PC1.tif$",full.names=T)))

#read as raster layers
ndvi<-map(file_list1,rast)
pca<-map(file_list2,rast)

#set values below 0 NDVI to 0 - barren land or water
for(i in 1:length(ndvi)){
  dat<-ndvi[[i]]
  dat[dat<0]<-0
  ndvi[[i]]<-dat
}

for(i in 1:length(ndvi)){
writeRaster(ndvi[[i]],paste0(path_out,"Final_outputs/Correll_NAIP/Z",i,"_zeroed_NDVI.tif"))
}

}
#read files if they exist
file_list1<-unlist(map(paste0(path_out,"Final_outputs/Correll_NAIP/"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))
file_list2<-unlist(map(paste0(path_out,"Final_outputs/Correll_NAIP/"),~list.files(.,pattern = "PC1.tif$",full.names=T)))
#read as raster layers
ndvi<-map(file_list1,rast)
pca<-map(file_list2,rast)



## Environmental Predictor 4: Wind speed/direction (better flood indicator than precipitation?)


## Environmental Predictor 5: angle to horizon line? Measure of tall surrounding topography or objects? 





## 4. Calculate texture surfaces **redo homogeneity zones 2-7
# --------------------------------------------
# Read raster if it already exists. Otherwise process the texture layers.
if(file.exists(paste0(path_out,"Final_outputs/NAIP_texture/hom_txt_1.tif"))){
  #list raster files
  file_list1<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture/"),~list.files(.,pattern = "hom.*[0-9].tif$",full.names=T)))
  file_list2<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture/"),~list.files(.,pattern = "cor.*[0-9].tif$",full.names=T)))
  file_list3<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture/"),~list.files(.,pattern = "ent.*[0-9].tif$",full.names=T)))
  
  #read as raster layers
  txt_homo<-map(file_list1,rast)
  txt_corr<-map(file_list2,rast)
  txt_entro<-map(file_list3,rast)
  

}
if(!file.exists(paste0(path_out,"Final_outputs/NAIP_texture/hom_txt_1.tif"))){
  
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
  
  # write newly generated layers to file


  for(i in 1:length(txt_homo)){
    writeRaster(txt_homo[[i]],filename=paste0(path_out,"Final_outputs/NAIP_texture/hom_txt",i,".tif"),overwrite=T)
    writeRaster(txt_entro[[i]],filename=paste0(path_out,"Final_outputs/NAIP_texture/ent_txt",i,".tif"),overwrite=T)
    writeRaster(txt_corr[[i]],file=paste0(path_out,"Final_outputs/NAIP_texture/cor_txt",i,".tif"),overwrite=T)
  }
}




