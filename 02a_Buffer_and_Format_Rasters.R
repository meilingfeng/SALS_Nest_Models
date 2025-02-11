
library(tidyverse)
library(terra)


#######################################################################################
## extend the edges of each raster so that all predictor variables cover all (most) nest locations
#######################################################################################



## Set up file paths
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"





## Load coarse resolution raster data
#----------------------------------------------------------------------------

## Use focal zone of 9x9 cells for coarse data (30m)
coarse_buff<-function(x){terra::focal(rast(x), w=9, fun=mean,na.policy="only", na.rm=T)}

## 1. UVVR
######
  # transformed uvvr coordinate system from WGS to NAD83 in ArcPro
  # just want band 4 (UVVR). Band 1 is unveg, 2 is veg, and 3 is proportion of water.

  # first remove outliers from UVVR data. Values should not exceed 2 or be under 0.
if(!file.exists(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier.tif"))){ # check if the file already exists to prevent overwriting it
  
  # A. UVVR mean across 2014-2018
  uvvr<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_mean'},.cols = everything()) 
  # values above 2 are not accurate, set to NA
  uvvr[uvvr>2]<-NA
  
  # B. UVVR change (2014-2018)
  uvvr14<-rast(paste0(dat_path,"Environmental Predictors/UVVR/uvvr_14_utm18_mean_ext.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_14'},.cols = everything())
  # values above 2 are not accurate, set to NA
  uvvr14[uvvr14>2]<-NA
  uvvr18<-rast(paste0(dat_path,"Environmental Predictors/UVVR/uvvr_18_utm18_mean_ext.tif"))%>%
    dplyr::rename_with(function(x){x<-'uvvr_18'},.cols = everything())
  uvvr18[uvvr18>2]<-NA
  # subtract 2018 from 2014 to get change
  uvvr_diff<-uvvr18-uvvr14
  uvvr_diff<-uvvr_diff%>%
    dplyr::rename_with(function(x){x<-'uvvr_diff'},.cols = everything())
  
  
  writeRaster(uvvr,filename = paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier.tif"),overwrite=T)
  writeRaster(uvvr_diff,filename = paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier.tif"),overwrite=T)
  
}

  # now buffer the UVVR rasters with outliers removed

if(!file.exists(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier_buff.tif"))){
  # create buffer
uvvr_dif<-coarse_buff(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier.tif"))
  # write file
writeRaster(uvvr_dif,paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier_buff.tif"),overwrite=T)

  # create buffer
uvvr_mean<-coarse_buff(paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier.tif"))
  # write file
writeRaster(uvvr_mean,paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier_buff.tif"),overwrite=T)

}





## 2. Precipitation
#####
if(!file.exists(paste0(path_out,"Intermediate_outputs/Precip/precip_buff.tif"))){
  
  # create buffer
  precip<-coarse_buff(paste0(dat_path,"Environmental Predictors/Precip/PRISM_ppt_30yr_normal_800mM4_annual_bil.bil"))
  # write file
  writeRaster(precip,paste0(path_out,"Intermediate_outputs/Precip/precip_buff.tif"),overwrite=T)
}





## 3. Tidal restrictions
#####
if(!file.exists(paste0(path_out,"Intermediate_outputs/Tidal_restriction/tideres_buff.tif"))){
  
  # create buffer
  tideres<-coarse_buff(paste0(dat_path,"Environmental Predictors/DSL_tidal_restrictions/tideres_2020_v5.0.tif"))
  # write file
  writeRaster(tideres,paste0(path_out,"Intermediate_outputs/Tidal_restriction/tideres_buff.tif"),overwrite=T)
}






## Load fine resolution raster data
#---------------------------------------------------------------------
## Use focal zone of 17x17 cells for fine data (3m)
fine_buff<-function(x){terra::focal(rast(x), w=17, fun=mean,na.policy="only", na.rm=T)}


## 1. Texture
#####

if(!file.exists(paste0(path_out,"Intermediate_outputs/Texture/ent_txt_1.tif"))){

  # first calculate neighborhood texture metrics
  
    # Raster Quanitization 
      # use equal range method when comparing across several datasets (can set gloabl min and max range using max_val and min_val) - splits data into equal ranges
      # alternative option is "equal prob" which splits by quantiles, used in original paper (Haralick and Shanmugam 1973)
      # NDVI of live plants ranges 0:1, set range and breaks at 0.1 intervals
  ndvi<-map(unlist(map(paste0(dat_path,"Environmental Predictors/Correll_NAIP/NDVI"),~list.files(.,pattern = "NDVI.tif$",full.names=T))),rast)
  txt_entro<-txt_corr<-list()
  
  for (i in 1:length(ndvi)){
    ndvi_rq<- quantize_raster(r = ndvi[[i]], n_levels = 10, method = "equal range")
    
    # Calculate texture metrics
      # use window size of 3x3m 
      # calculate across multiple shifts (invariant)
    txt<- glcm_textures(ndvi_rq, w = c(3,3), n_levels =10, 
                        metrics = c("glcm_entropy", "glcm_correlation"),
                        quantization = "none", shift = list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)))
    
    txt_entro[[i]]<-txt[[1]]
    txt_corr[[i]]<-txt[[2]]
    
  }
  
  # write newly generated layers to file
  
  for(i in 1:length(txt_entro)){
    writeRaster(txt_entro[[i]],filename=paste0(path_out,"Intermediate_outputs/Texture/ent_txt",i,".tif"),overwrite=T)
    writeRaster(txt_corr[[i]],file=paste0(path_out,"Intermediate_outputs/Texture/cor_txt",i,".tif"),overwrite=T)
  }
}


  # then buffer the rasters
if(!file.exists(paste0(path_out,"Intermediate_outputs/Texture/ent_txt_1_buff.tif"))){
  
  # list files
cor_list<-unlist(map(paste0(path_out,"Intermediate_outputs/Texture"),~list.files(.,pattern = "^cor_txt_[0-9].tif$",full.names=T)))
  # apply buffer
cor_list<-map(cor_list,fine_buff)
  #write buffer files
for(i in 1:length(cor_list)){
  writeRaster(cor_list[[i]],paste0(path_out,"Intermediate_outputs/Texture/cor_txt_",i,"_buff.tif"))
}

  # list files
ent_list<-unlist(map(paste0(path_out,"Intermediate_outputs/Texture"),~list.files(.,pattern = "^ent_txt_[0-9].tif$",full.names=T)))
  # apply buffer
ent_list<-map(ent_list,fine_buff)
  # write buffer files
for(i in 1:length(ent_list)){
  writeRaster(ent_list[[i]],paste0(path_out,"Intermediate_outputs/Texture/ent_txt_",i,"_buff.tif"))
}

}



## 2. NDVI
#####

  # first buffer the original layer
if(!file.exists(paste0(path_out,"Intermediate_outputs/NDVI/1_NDVI_buff.tif"))){
  # list files
  ndvi_list<-unlist(map(paste0(dat_path,"Environmental Predictors/Correll_NAIP"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))
  # apply buffer
  ndvi_list<-map(ndvi_list,fine_buff)
  #write buffer files
  for(i in 1:length(ndvi_list)){
      writeRaster(ndvi_list[[i]],paste0(path_out,"Intermediate_outputs/NDVI/",i,"_NDVI_buff.tif"),overwrite=T)
    }
}

  # then also set all negative values in NDVI to zero. All negative values are non-plant matter.
if(!file.exists(paste0(path_out,"Intermediate_outputs/NDVI/Z1_zeroed_NDVI_buff.tif"))){
  ndvi_list<-map(unlist(map(paste0(path_out,"Intermediate_outputs"),~list.files(.,pattern = "NDVI_buff.tif$",full.names=T))),rast)
  #set values below 0 NDVI to 0 - barren land or water
  for(i in 1:length(ndvi_list)){
    dat<-ndvi_list[[i]]
    dat[dat<0]<-0
    ndvi_list[[i]]<-dat
  }
  for(i in 1:length(ndvi_list)){
    writeRaster(ndvi_list[[i]],paste0(path_out,"Intermediate_outputs/NDVI/Z",i,"_zeroed_NDVI_buff.tif"),overwrite=T)
  }
}


## 3. PCA *** continue here
#####

  # PCA raster has multiple layers, one for each principal component. We just want the first layer, PC1.
  # First extract PC1 to its own raster layer
if(!file.exists(paste0(path_out,"Intermediate_outputs/PCA/Z1_PC1.tif"))){
  # list files
  pca_list<-unlist(map(paste0(dat_path,"Environmental Predictors/Correll_NAIP"),~list.files(.,pattern = "PCA.tif$",full.names=T)))
  #extract just the first band for the pca data
  mapply(function(x,y){writeRaster(rast(x,lyrs=1),paste0(path_out,"Intermediate_outputs/PCA/Z",y,"_PC1.tif"),overwrite=T)},
         pca_list,c(1:length(pca_list)))
}

  # Now run buffer
if(!file.exists(paste0(path_out,"Intermediate_outputs/PCA/Z1_PC1_buff.tif"))){
  # list files
  pca_list<-unlist(map(paste0(path_out,"Intermediate_outputs/PCA"),~list.files(.,pattern = "PC1.tif$",full.names=T)))
  # apply buffer
  pca_list<-map(pca_list,fine_buff)
  #write buffer files
  for(i in 1:length(pca_list)){
    writeRaster(pca_list[[i]],paste0(path_out,"Intermediate_outputs/PCA/Z",i,"_PC1_buff.tif"),overwrite=T)
  }
}





## 4. High marsh
#####
if(!file.exists(paste0(path_out,"Intermediate_outputs/HIMARSH/Z1_vegcls_buff.tif"))){
#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

#read as raster layers
vg_cls<-map(file_list,rast)

mode_buff<-function(x){
  out<-terra::focal(x, w=9, fun=modal, na.policy="only",na.rm=T)
}

vg_cls<-map(vg_cls,mode_buff)

for (i in 1:length(vg_cls)){
  writeRaster(vg_cls[[i]],paste0(path_out,"Intermediate_outputs/HIMARSH/Z",i,"_vegcls_buff.tif"))
}

}
t<-mode_buff(vg_cls[[5]])

## 5. Elevation (USGS 1m DEM rescaled to 3m in Correll et al. 2018)
#####

if(!file.exists(paste0(path_out,"Intermediate_outputs/Elevation/Z1_DEM_buff.tif"))){
  # list files
  dem_list<-unlist(map(paste0(dat_path,"Environmental Predictors/Correll_DEM_RS"),~list.files(.,pattern = "DEMrs.tif$",full.names=T)))
  
  #read as raster layers
  dem<-map(dem_list,rast)
  
  # Now run buffer
  dem_list<-map(dem_list,fine_buff)
  #write buffer files
  for(i in 1:length(dem_list)){
    writeRaster(dem_list[[i]],paste0(path_out,"Intermediate_outputs/Elevation/Z",i,"_DEM_buff.tif"),overwrite=T)
  }
}

## 5. Elevation (USGS 1m DEM rescaled to 3m in Correll et al. 2018)
#####

if(!file.exists(paste0(path_out,"Intermediate_outputs/Elevation/Z1_DEM_buff.tif"))){
  # list files
  dem_list<-unlist(map(paste0(dat_path,"Environmental Predictors/Correll_DEM_RS"),~list.files(.,pattern = "DEMrs.tif$",full.names=T)))
  
  #read as raster layers
  dem<-map(dem_list,rast)

  # Now run buffer
  dem_list<-map(dem_list,fine_buff)
  #write buffer files
  for(i in 1:length(dem_list)){
    writeRaster(dem_list[[i]],paste0(path_out,"Intermediate_outputs/Elevation/Z",i,"_DEM_buff.tif"),overwrite=T)
  }
}
