## 6. Create prediction surface
# ---------------------------------------------------------------------------------------------
# align extent and resolution of all environmental rasters
# use a prediction surface resolution of 3 meters (matches the Correll data, plus middle ground between NAIP and UVVR)
# Set all coordinate systems to ESPG:4269 (NAD 1983 geographic coord system)
# Use Atlantic coast marsh extent (refer to Correll marsh layer) for each zone (n = 8)- predict by zone to reduce computation time

#file path names
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#functions
source("C:/Users/mefen/OneDrive/Documents/Github/UCONN/SHARP/Functions/randomPoints.R")

#nest data
nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_10_20_fail_thin_1m_dist_errors_removed.shp"))



## 1. Create prediction surface templates for each zone
#---------------------------------------------------------------------------

# Get zone extents from Correll Marsh Veg data

#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))
#read as raster layers
vg_cls<-map(file_list,rast)

# create templates from each zone
reso<-30 #resolution
temps<-c() #list to hold the zone templates
temps<-lapply(vg_cls,function(x)rast(extent=ext(x),crs="EPSG:26918",resolution=reso,vals=1))

vg_cls[[1]]

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
