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
dat_path<-"D:/Data/"
path_out<-"D:/Outputs/"



##  Load in data

#list raster files
file_list1<-unlist(map(paste0(dat_path,"Correll_NAIP/"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))
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


## Calculate texture surfaces **redo homogeneity zones 2-7
# --------------------------------------------

if(!file.exists(paste0(path_out,"hom_txt_1.tif"))){
  
  ## 3. a. Raster Quanitization 
  #----------------------------
  # use equal range method when comparing across several datasets (can set gloabl min and max range using max_val and min_val) - splits data into equal ranges
  # alternative option is "equal prob" which splits by quantiles, used in original paper (Haralick and Shanmugam 1973)
  # NDVI of live plants ranges 0:1, set range and breaks at 0.1 intervals
  txt_homo<-list()
  
  for (i in 1:length(ndvi)){
    ndvi_rq<- quantize_raster(r = ndvi[[i]], n_levels = 10, method = "equal range")
    
    ## 3. b. calc texture metric surfaces
    #----------------------------
    # use window size of 3x3m 
    # calculate across multiple shifts (invariant)
    txt<- glcm_textures(ndvi_rq, w = c(3,3), n_levels =10, 
                        metrics = c("glcm_homogeneity"),
                        quantization = "none", shift = list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)))
    
    txt_homo[[i]]<-txt[[1]]
    
  }
  
}



# If files don't already exist, write newly generated layers to file

if(!file.exists(paste0(path_out,"hom_txt_1.tif"))){
  for(i in 1:length(txt_homo)){
    writeRaster(txt_homo[[i]],filename=paste0(path_out,"hom_txt_",i,".tif"),overwrite=T)
  }
}


