library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(exactextractr)
#gdal, for gdalwarp (manupilating raster projections)
#first download gdal here https://gdal.org/download.html#current-release
#unzip gdal using this technique onto your C drive https://www.makeuseof.com/what-is-gz-file/
#library(devtools)
#install_github("JoshOBrien/gdalUtilities")
library(gdalUtilities)#raster manipulation


## Create prediction surface
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
#use UTM18 NAD83 coordinate system
temps<-lapply(vg_cls,function(x)rast(extent=ext(x),crs="EPSG:26918",resolution=reso,vals=1))




## 2. load all the remaining predictors
#---------------------------------------------------------------------------
#UVVR
uvvr_dif<-paste0(path_out,"Final_outputs/UVVR/uvvr_diff_noOutlier.tif")
uvvr_mean<-paste0(path_out,"Final_outputs/UVVR/uvvr_mean_noOutlier.tif")

#NAIP
#list raster files
cor_list<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture"),~list.files(.,pattern = "^cor_txt_[0-9].tif$",full.names=T)))
ent_list<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture"),~list.files(.,pattern = "^ent_txt_[0-9].tif$",full.names=T)))
hom_list<-unlist(map(paste0(path_out,"Final_outputs/NAIP_texture"),~list.files(.,pattern = "^hom_txt_[0-9].tif$",full.names=T)))

if(!file.exists(paste0(path_out,"Final_outputs/Correll_NAIP/Z1_zeroed_NDVI.tif"))){
ndvi_list<-map(unlist(map(paste0(dat_path,"Correll_NAIP"),~list.files(.,pattern = "NDVI.tif$",full.names=T))),rast)
#set values below 0 NDVI to 0 - barren land or water
for(i in 1:length(ndvi_list)){
  dat<-ndvi_list[[i]]
  dat[dat<0]<-0
  ndvi_list[[i]]<-dat
}
for(i in 1:length(ndvi_list)){
  writeRaster(ndvi_list[[i]],paste0(path_out,"Final_outputs/Correll_NAIP/Z",i,"_zeroed_NDVI.tif"))
}
}

ndvi_list<-unlist(map(paste0(path_out,"Final_outputs/Correll_NAIP"),~list.files(.,pattern = "NDVI.tif$",full.names=T)))

if(!file.exists(paste0(path_out,"Final_outputs/Correll_NAIP/Z1_PC1.tif"))){
pca_list<-unlist(map(paste0(dat_path,"Correll_NAIP"),~list.files(.,pattern = "PCA.tif$",full.names=T)))
  #extract just the first band for the pca data
mapply(function(x,y){writeRaster(rast(x,lyrs=1),paste0(path_out,"Final_outputs/Correll_NAIP/Z",y,"_PC1.tif"))},
       pca_list,c(1:length(pca_list)))
}

pca_list<-unlist(map(paste0(path_out,"Final_outputs/Correll_NAIP"),~list.files(.,pattern = "PC1.tif$",full.names=T)))


#Precipitation
precip<-paste0(dat_path,"Precip/PRISM_ppt_30yr_normal_800mM4_annual_bil.bil")



#combine data that are at the full extent
file_list1<-c(uvvr_dif,uvvr_mean,precip)

#combine data that are at zone extents
files2<-list(file_list,cor_list,ent_list,ndvi_list,pca_list)
file_list2<-unlist(lapply(files2,sort))



## 3. Align full extent predictors to templates
#---------------------------------------------------------------------------

#for each full extent dataset
for(j in 2:length(file_list1)){
  
  input<-file_list1[[j]]

  # First set the resampling method and output data type to match the input data
  #####
    #if the data is numeric (decimals), use average resampling
    #if the data is nominal (classes), use mode resampling
    # on spat raster data types https://rdrr.io/cran/raster/man/dataType.html 
  if(grepl("INT",rast(input)@ptr[["dataType"]],fixed=T)){
    rs_method<-"mode"
    dat_type<-"UInt16"
  }
  if(grepl("FLT",rast(input)@ptr[["dataType"]])){
    rs_method<-"average"
    dat_type<-"Float32"
  }
  
  
  # Then align the input data to each of the zone templates
  #####
  for(i in 1:length(temps)){
    
    #create output file name
    output<-paste(substring(file_list1[[j]], 
                            #take the spot just before .tif in the file name
                            c(1,nchar(file_list1[[j]])-3), 
                            c(nchar(file_list1[[j]])-4,nchar(file_list1[[j]]))), 
                  # and add in the aligned resolution and zone to the output name
                  collapse=paste0("_",reso,"align_Z",i))

    
    #if the file has not been already aligned, align it
    if (!file.exists(output)) {
    
    
  gdalwarp(srcfile=input, # original file name
           dstfile=output, #aligned output file name
           t_srs=as.character(crs(temps[[i]])), tr=res(temps[[i]]), # aligned coordinate system, aligned resolution
           te=c(xmin(temps[[i]]), ymin(temps[[i]]), xmax(temps[[i]]), ymax(temps[[i]])), #aligned extent
           r=rs_method, # resampling method ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")
           dstnodata="None", # no data value
           of='GTiff', ot=dat_type, overwrite=TRUE)# output data type, Byte, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
    }
  
  }
}





## 4. Align zone extent predictors to templates
#---------------------------------------------------------------------------

#define alignment function
align_z_rast<-function(in_file,temp_list,collapse_string){
  
  input<-in_file
  
  # 1. first set the resampling method and output data type to match the input data
  #####
  #if the data is numeric (decimals), use average resampling
  #if the data is nominal (classes), use mode resampling
  # on spat raster data types https://rdrr.io/cran/raster/man/dataType.html 
  if(grepl("INT",rast(input)@ptr[["dataType"]],fixed=T)){
    rs_method<-"mode"
    dat_type<-"UInt16"
  }
  if(grepl("FLT",rast(input)@ptr[["dataType"]])){
    rs_method<-"average"
    dat_type<-"Float32"
  }
  
  
  # 2. Then align the input data to each of the zone templates
  #####
    
    #create output file name
    output<-paste(substring(input, 
                            #take the spot just before .tif in the file name
                            c(1,nchar(input)-3), 
                            c(nchar(input)-4,nchar(input))), 
                  # and add in the aligned resolution and zone to the output name
                  collapse=collapse_string)
    
    
    #if the file has not been already aligned, align it
    if (!file.exists(output)) {
      
      
      gdalwarp(srcfile=input, # original file name
               dstfile=output, #aligned output file name
               t_srs=as.character(crs(temp_list)), tr=res(temp_list), # aligned coordinate system, aligned resolution
               te=c(xmin(temp_list), ymin(temp_list), xmax(temp_list), ymax(temp_list)), #aligned extent
               r=rs_method, # resampling method ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")
               dstnodata="None", # no data value
               of='GTiff', ot=dat_type, overwrite=TRUE)# output data type, Byte, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
    }
    
}


#create a repeating list of templates for each dataset
temp_order<-rep(temps,times = length(files2))
#apply the raster alignment to each dataset zone and template zone pair
mapply(align_z_rast,file_list2,temp_order,collapse_string=paste0("_",reso,"align"))




## 3. Align vegetation class proportions to templates
#---------------------------------------------------------------------------
#for the veg class data, also calculate the proportion of each class at the higher resolution

#class values in layers
#HIMARSH = 1
#LOMARSH=2
#MUDFLAT=4
#PHRAG=5
#TERRBRD=8
#STREAM=7
#POOL/PANNE=6

#(focus on HIMARSH, LOMARSH for now)
hi_files<-list()
lo_files<-list()

for(i in 1:length(vg_cls)){
#write the new raster with veg class proportions to file
hi_files[[i]]<-paste0(path_out,"Final_outputs/Correll_Marsh_Zones/Z",i,"_himarsh_",reso,".tif")
lo_files[[i]]<-paste0(path_out,"Final_outputs/Correll_Marsh_Zones/Z",i,"_lomarsh_",reso,".tif")

#if the file does not exist...
if (!file.exists(hi_files[[i]])) {
#break the veg layer into its classes by setting each veg class value to 1 separately
hi<-lo<-vg_cls[[i]]
hi[hi!=1]<-0
lo[lo!=2]<-0

#count the number of 3x3m cells within each larger raster cell 30x30, using aggregate
hi<-terra::aggregate(hi,fact=10,fun='sum',dissolve=F)
lo<-terra::aggregate(hi,fact=10,fun='sum',dissolve=F)

#divide count by 100 to get proportion
hi<-hi/100
lo<-lo/100

#give the variable a name
hi<-dplyr::rename_with(hi,function(x){x<-'HIMARSH'},.cols = everything())
lo<-dplyr::rename_with(hi,function(x){x<-'LOMARSH'},.cols = everything())

#write the new raster with veg class proportions to file
writeRaster(hi,filename=hi_files[[i]], overwrite=TRUE)
writeRaster(lo,filename=lo_files[[i]], overwrite=TRUE)
}

}

#then align these to the templates
#apply the raster alignment to each dataset zone and template zone pair
mapply(align_z_rast,hi_files,temps,"align")
mapply(align_z_rast,lo_files,temps,"align")
  







## 5. combine all prediction varaibles into a single dataframe by cell ID
#---------------------------------------------------------------------------
#Bring in all the aligned data
file_list1_align<-list()
file_list2_align<-list()
hi_align<-c()
lo_align<-c()
set_t<-c()

#full extent datasets
  for (j in 1:length(file_list1)){
  for(i in 1:length(temps)){
  
  #create output file name
  set_t[i]<-paste(substring(file_list1[[j]], 
                          #take the spot just before .tif in the file name
                          c(1,nchar(file_list1[[j]])-3), 
                          c(nchar(file_list1[[j]])-4,nchar(file_list1[[j]]))), 
                # and add in the aligned resolution and zone to the output name
                collapse=paste0("_",reso,"align_Z",i))

  }
    file_list1_align[[j]]<-set_t
  }
file_list1_align<-unlist(file_list1_align)


#zone extent datasets
file_list2_align<-unlist(map(file_list2,function(input){
  paste(substring(input, 
                #take the spot just before .tif in the file name
                c(1,nchar(input)-3), 
                c(nchar(input)-4,nchar(input))), 
      # and add in the aligned resolution and zone to the output name
      collapse=paste0("_",reso,"align"))
}
))


# veg proportion data
hi_align<-unlist(map(hi_files,function(input){
  paste(substring(input, 
                       #take the spot just before .tif in the file name
                       c(1,nchar(input)-3), 
                       c(nchar(input)-4,nchar(input))), 
             # and add in the aligned resolution and zone to the output name
             collapse="align")
}
))

lo_align<-unlist(map(lo_files,function(input){
  paste(substring(input, 
                       #take the spot just before .tif in the file name
                       c(1,nchar(input)-3), 
                       c(nchar(input)-4,nchar(input))), 
             # and add in the aligned resolution and zone to the output name
             collapse="align")
}
))

all_files<-c(file_list2_align,file_list1_align,hi_align,lo_align)



#give each template cell a unique id
ids<-map(all_files[c(1:8)],rast)
for(i in 1:length(ids)){
  id<-ids[[i]]
  id[id>0]<-c(1:ncell(id[id>0]))
  id[id==0]<-NA
  ids[[i]]<-id}




mat<-list()
file_list_all_zones<-list()

#for each zone
for(j in 2:length(temps)){
  file_list_zone<-list()
# and for each dataset in that zone
for(i in seq(1,length(all_files),by=8)){
  #go to each unique dataset i (length of zones) and take zone j
  file_list_zone[i]<-all_files[(i+j-1)]
  file_list_zone<-file_list_zone[file_list_zone!="NULL"]
}
  dat_zone<-rast(unlist(file_list_zone))
  dat_zone<-rast(list(dat_zone,ids[[j]]))
  dat_zone<-mask(dat_zone,ids[[j]])
  names(dat_zone)<-c("vg_clss","cor_txt","ent_txt","NDVI","PC1","uvvr_diff","uvvr_mean","precip","himarsh","lomarsh","id")
  mat[[j]]<-as.data.frame(terra::as.matrix(dat_zone,wide=F))%>%
    mutate(id=paste0(id,"z",j))%>%
    filter(vg_clss!="NaN")
  file_list_all_zones[[j]]<-file_list_zone
}


#zone 1 is too big, process it in 2 parts
n_div<-4
dat_zone_list<-list()
row_start<-c()
row_end<-c()
file_list_zone<-list()

  # for each dataset in that zone, make a list of environmental data files
  for(i in seq(1,length(all_files),by=8)){
    #go to each unique dataset i (length of zones) and take zone j
    file_list_zone[i]<-all_files[i]
    file_list_zone<-file_list_zone[file_list_zone!="NULL"]
  }
  #read the list of data files for that zone into a raster stack
  dat_zone<-rast(unlist(file_list_zone))
  #add the cell IDs to the stack
  dat_zone<-rast(list(dat_zone,ids[[1]]))
  #set all cells outside the marsh to NA
  dat_zone<-mask(dat_zone,ids[[1]])
  #label the datasets
  names(dat_zone)<-c("vg_clss","cor_txt","ent_txt","NDVI","PC1","uvvr_diff","uvvr_mean","precip","himarsh","lomarsh","id")
  #divide zone into 4
  for(i in 1:n_div){
   row_start[i]<-((round(nrow(dat_zone)/n_div))*(i-1))+1
   row_end[i]<-(round(nrow(dat_zone)/n_div))*i
  }
  for(i in 1:n_div){
  dat_zone_sub<-dat_zone[c(row_start[i]:row_end[i]),c(1:ncol(dat_zone)),drop=F]
  
  #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
  dat_zone_list[[i]]<-as.data.frame(terra::as.matrix(dat_zone_sub,wide=F))%>%
    mutate(id=paste0(id,"z1"))%>%
    #remove cells outside the marsh boundary
    filter(vg_clss!="NaN")
  }
  #add the list of data files for that zone to a list of all zones' datafiles
  file_list_all_zones[[1]]<-file_list_zone

  #bind all the zone sections into one dataframe
  mat[[1]]<-do.call("rbind",dat_zone_list)


  
#write each zone prediction dataset to file
for(i in 1:length(mat)){
  write.csv(mat[[i]],file=paste0(path_out,"/Final_outputs/prediction_surfaces/Z",i,"_prediction_30m.csv"),row.names = F)
}

#save the list of predictor files for each zone as an R object
save(file_list_all_zones,file=paste0(path_out,"/predictor_files_all_zones_",reso,"m.rds"))

