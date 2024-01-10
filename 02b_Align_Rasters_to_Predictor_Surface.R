library(tidyverse)
#spatial analysis
library(sf)
#raster analysis
library(terra)
library(tidyterra)
#library(exactextractr)# 
#gdal, for gdalwarp (manupilating/aligning raster projections)
#first download gdal here https://gdal.org/download.html#current-release
#unzip gdal using this technique onto your C drive https://www.makeuseof.com/what-is-gz-file/
#library(devtools)
#install_github("JoshOBrien/gdalUtilities")
library(gdalUtilities)#raster manipulation

##################################################################################################################
## Create prediction surface
# ---------------------------------------------------------------------------------------------
# align extent and resolution of all environmental rasters
# use a prediction surface resolution of 30 meters based on min distance between nests and available computing power
# Set all coordinate systems to ESPG:4269 (NAD 1983 geographic coord system)
# Use Atlantic coast marsh extent (Correll marsh layer) for each zone (n = 8)- predict by zone to reduce computation time
######################################################################################################################


## Set up
#--------------------------------------------------------------------------
#file path names
dat_path<-"C:/Users/10788/Desktop/SaltMarsh/Data/"
path_out<-"C:/Users/10788/Desktop/SaltMarsh/Outputs/"


#functions (to select random background points)
source("C:/Users/10788/Desktop/SaltMarsh/SHARP/Functions/randomPoints.R")

# cleaned and filtered focal species nest data
nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_dist_err_removed.shp"))



## 1. Create prediction surface templates for each zone (this is what we will align the other raster layers to)
#---------------------------------------------------------------------------

# Get the zone extents from Correll Marsh Veg data

#list raster files
file_list<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))

#read as raster layers
vg_cls<-map(file_list,rast)

# create templates from each zone
reso<-30 #resolution
temps<-c() #list to hold the zone templates
#use UTM18 NAD83 coordinate system
temps<-lapply(vg_cls,function(x)rast(extent=ext(x),crs="EPSG:26918",resolution=reso,vals=1))




## 2. load all the environmental predictors
#---------------------------------------------------------------------------
#UVVR 
uvvr_dif<-paste0(path_out,"Intermediate_outputs/UVVR/uvvr_diff_noOutlier_buff.tif")
uvvr_mean<-paste0(path_out,"Intermediate_outputs/UVVR/uvvr_mean_noOutlier_buff.tif")

#Precipitation
precip<-paste0(path_out,"Intermediate_outputs/Precip/precip_buff.tif")

# Tidal restrictions
tideres<-paste0(path_out,"Intermediate_outputs/Tidal_restriction/tideres_buff.tif")

# Texture
#list raster files
cor_list<-unlist(map(paste0(path_out,"Intermediate_outputs/Texture"),~list.files(.,pattern = "^cor_txt_[0-9]_buff.tif$",full.names=T)))
ent_list<-unlist(map(paste0(path_out,"Intermediate_outputs/Texture"),~list.files(.,pattern = "^ent_txt_[0-9]_buff.tif$",full.names=T)))

# NDVI
ndvi_list<-unlist(map(paste0(path_out,"Intermediate_outputs/NDVI"),~list.files(.,pattern = "_zeroed_NDVI_buff.tif$",full.names=T)))

# Principal Component of raw NAIP reflection bands (R, G, B)
pca_list<-unlist(map(paste0(path_out,"Intermediate_outputs/PCA"),~list.files(.,pattern = "PC1.tif$",full.names=T)))

# Elevation
dem_list<-unlist(map(paste0(path_out,"Intermediate_outputs/Elevation"),~list.files(.,pattern = "DEM_buff.tif$",full.names=T)))


#combine data that are at the full range extent
file_list1<-c(uvvr_dif,uvvr_mean,precip,tideres)

#combine data that are at zone extents
files2<-list(file_list,cor_list,ent_list,ndvi_list,pca_list,dem_list)
file_list2<-unlist(lapply(files2,sort))



## 3. Align full extent predictors to zone templates
#---------------------------------------------------------------------------

#for each full extent dataset...
for(j in 1:length(file_list1)){
  
  input<-file_list1[[j]] # take 1 individual dataset
  
  # Align the input data to each of the zone templates
  #####
    # for each zone template...
  for(i in 1:length(temps)){
    
    #create an output file name by referencing the input dataset
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
           r="average", # resampling method ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")
           dstnodata="None", # no data value
           of='GTiff', ot="Float32", overwrite=TRUE)# output data type, Byte, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
    }
  
  }
}





## 4. Align zone extent predictors to zone templates
#---------------------------------------------------------------------------

#define alignment function
align_z_rast_average<-function(in_file,temp_list,collapse_string){
  
  input<-in_file
  
  
  # Align the input data to each of the zone templates
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
               r="average", # resampling method ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")
               dstnodata="None", # no data value
               of='GTiff', ot="Float32", overwrite=TRUE)# output data type, Byte, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
    }
    
}

#mode version for the vegetation classes
align_z_rast_mode<-function(in_file,temp_list,collapse_string){
  
  input<-in_file
  
  
  # Align the input data to each of the zone templates
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
             r="mode", # resampling method ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")
             dstnodata="None", # no data value
             of='GTiff', ot="Int8", overwrite=TRUE)# output data type, Byte, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
  }
  
}


#create a repeating list of templates for each dataset
temp_order<-rep(temps,times = length(files2)-1)

#apply the raster alignment function to each dataset and template zone pair
mapply(align_z_rast_average,file_list2[-c(1:8)],temp_order,collapse_string=paste0("_",reso,"align"))
mapply(align_z_rast_mode,file_list2[c(1:8)],temps,collapse_string=paste0("_",reso,"align"))




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

#(focus on HIMARSH for now)

hi_files<-list()

if (!file.exists(paste0(path_out,"Intermediate_outputs/HIMARSH/Z1_himarsh_",reso,".tif"))) {
  
for(i in 1:length(vg_cls)){
# distinguish high marsh by setting high marsh values to 1 and all other values to 0
hi<-vg_cls[[i]]
hi[hi!=1]<-0

#count the number of 3x3m cells within each larger raster cell (30x30), using aggregate
hi<-terra::aggregate(hi,fact=10,fun='sum',dissolve=F)

#divide count by 100 to get proportion
hi<-hi/100

#give the variable a name
hi<-dplyr::rename_with(hi,function(x){x<-'HIMARSH'},.cols = everything())

#create output file names for the high marsh proportions 
hi_files[[i]]<-paste0(path_out,"Intermediate_outputs/HIMARSH/Z",i,"_himarsh_",reso,".tif")

#write the new raster with veg class proportions to file
writeRaster(hi,filename=hi_files[[i]], overwrite=TRUE)
}

}

for(i in 1:length(vg_cls)){
hi_files[[i]]<-paste0(path_out,"Intermediate_outputs/HIMARSH/Z",i,"_himarsh_",reso,".tif")
}

#then align these to the templates
#apply the raster alignment to each dataset zone and template zone pair
mapply(align_z_rast_average,hi_files,temps,"align")




## 5. combine all prediction variables into a single dataframe by cell ID
#---------------------------------------------------------------------------
#Bring in all the aligned data
file_list1_align<-list()
file_list2_align<-list()
hi_align<-c()
set_t<-c()

# get a list of file names for full extent datasets
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


# get a list of file names for zone extent datasets
file_list2_align<-unlist(map(file_list2,function(input){
  paste(substring(input, 
                #take the spot just before .tif in the file name
                c(1,nchar(input)-3), 
                c(nchar(input)-4,nchar(input))), 
      # and add in the aligned resolution and zone to the output name
      collapse=paste0("_",reso,"align"))
}
))


# get high marsh proportion file names
hi_align<-unlist(map(hi_files,function(input){
  paste(substring(input, 
                       #take the spot just before .tif in the file name
                       c(1,nchar(input)-3), 
                       c(nchar(input)-4,nchar(input))), 
             # and add in the aligned resolution and zone to the output name
             collapse="align")
}
))


# combine the file names of all aligned predictor variables into a list
all_files<-c(file_list2_align,file_list1_align,hi_align)



#give each predictor surface template cell a unique id
ids<-map(all_files[c(1:8)],rast)
for(i in 1:length(ids)){
  id<-ids[[i]]
  id[!(id%in%c(0,7,9))]<-c(1:ncell(id[!(id%in%c(0,7,9))])) #remove marsh buffer zones
  id[id%in%c(0,7,9)]<-NA
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
  #dat_zone<-rast(unlist(file_list_zone))
  #dat_zone<-rast(list(dat_zone,ids[[j]]))
  #dat_zone<-mask(dat_zone,ids[[j]])
  #names(dat_zone)<-c("vg_clss","cor_txt","ent_txt","NDVI","PC1","uvvr_diff","uvvr_mean","precip","tideres","himarsh","id")
  #mat[[j]]<-as.data.frame(terra::as.matrix(dat_zone[-c(3,6,8)],wide=F))%>%#remove predictors not being used (cor, uvvr_diff, precip)
  #  mutate(id=paste0(id,"z",j))%>%
  #  filter(vg_clss!="NaN")
  file_list_all_zones[[j]]<-file_list_zone
}


#zone 1 is too big, process it in parts
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
#  dat_zone<-rast(unlist(file_list_zone))
  #add the cell IDs to the stack
#  dat_zone<-rast(list(dat_zone,ids[[1]]))
  #set all cells outside the marsh to NA
#  dat_zone<-mask(dat_zone,ids[[1]])
  #label the datasets
#  names(dat_zone)<-c("vg_clss","cor_txt","ent_txt","NDVI","PC1","uvvr_diff","uvvr_mean","precip","tideres","himarsh","id")
  #divide zone into parts
#  for(i in 1:n_div){
#   row_start[i]<-((round(nrow(dat_zone)/n_div))*(i-1))+1
#   row_end[i]<-(round(nrow(dat_zone)/n_div))*i
#  }
#  for(i in 1:n_div){
#  dat_zone_sub<-dat_zone[c(row_start[i]:row_end[i]),c(1:ncol(dat_zone)),drop=F] # if it throws an error, might need to sub 1 from the last row end
  
  #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
#  dat_zone_list[[i]]<-as.data.frame(terra::as.matrix(dat_zone_sub[-c(3,6,8)],wide=F))%>% #remove variables we arent using
#    mutate(id=paste0(id,"z1"))%>%
    #remove cells outside the marsh boundary
#    filter(vg_clss!="NaN")
#  }
  #add the list of data files for that zone to a list of all zones' datafiles
  file_list_all_zones[[1]]<-file_list_zone

  #bind all the zone sections into one dataframe
#  mat[[1]]<-do.call("rbind",dat_zone_list)


  
#write each zone prediction dataset to file
#for(i in 1:length(mat)){
#  write.csv(mat[[i]],file=paste0(path_out,"/Final_outputs/prediction_surfaces/Z",i,"_prediction_30m.csv"),row.names = F)
#}

#save the list of predictor files for each zone as an R object
save(file_list_all_zones,file=paste0(path_out,"/predictor_files_all_zones_",reso,"m.rds"))

