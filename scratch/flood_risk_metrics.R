library(tidyverse)
library(terra)


#######################################################################################
## extend the edges of each raster so that all predictor variables cover all (most) nest locations
#######################################################################################



## Set up file paths
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/" #input data
path_out<-"D:/Nest_Models/Outputs/" #outputs

## Get NOAA Co-Ops tide gauge datums and exceedance levels
  # gauges that intersect or come close to the Correll 2018 vegetation layer
names<-c("Eastport","Cutler Farris Wharf","Bar Harbor","Portland","Seavey Island","Boston","Providence","Conimicut Light",
         "Chatham","Fall River","Quonset Point","Newport","Woods Hole","New London","New Haven","Nantucket Island",
         "Bridgeport","Montauk","Kings Point","The Battery","Bergen Point West Reach","Sandy Hook","Atlantic City",
         "Marcus Hook","Reedy Point","Delaware City","Cape May","Lewes","Cambridge","Bishops Head","Ocean City Inlet",
         "Wachapreague","Kiptopeke")
states<-c("ME","ME","ME","ME","ME","MA","RI","RI","MA","MA","RI","RI","MA","CT","CT","MA","CT","NY","NY","NY","NY","NJ",
         "NJ","PA","DE","DE","NJ","DE","MD","MD","MD","VA","VA")
codes<-c(8410140,8411060,8413320,8418150,8419870,8443970,8454000,8452944,8447435,8447386,8454049,8452660,8447930,8461490,
        8465705,8449130,8467150,8510560,8516945,8518750,8519483,8531680,8534720,8540433,8551910,8551762,8536110,8557380,
        8571892,8571421,8570283,8631044,8632200)
exceed<-c() # meters above MHHW expected 1% of the time
max<-c()
hat<-c()
mhhw<-c()
mhw<-c()
msl<-c()

gauge<-data.frame(name=names,state=states,code=codes)
coops_search(station_name = gauge$codes, begin_date = 20150927, end_date = 20150928,
             product = "daily_mean", datum = "stnd", time_zone = "lst")




## Load fine resolution raster data
#---------------------------------------------------------------------
## Use focal zone of 17x17 cells for fine data (3m)
fine_buff<-function(x){terra::focal(rast(x), w=17, fun=mean,na.policy="only", na.rm=T)}

## 5. DEM and tide data (tide level for habitat class and tide exceedance for flood risk)
#####
if(!file.exists(paste0(path_out,"Intermediate_outputs/Tides/Z1_tide_level_buff.tif"))){
  #list raster files (changed Zone 4 name to DEMrs but originally just DEM)
  file_list<-unlist(map(paste0(dat_path,"Environmental Predictors/Correll_DEM_RS"),~list.files(.,pattern = "DEMrs.tif$",full.names=T)))
  
  #read as raster layers
  dem<-map(file_list,rast)
  
  # apply buffer
  dem_list<-map(dem,fine_buff)
  
  for (i in 1:length(vg_cls)){
    writeRaster(vg_cls[[i]],paste0(path_out,"Intermediate_outputs/Tides/Z",i,"_tide_level_buff.tif"))
  }
  
}
