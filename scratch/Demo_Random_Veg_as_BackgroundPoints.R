

# load cleaned focal species nest observation shapefile
nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[j],"_valid_nest_locations_2010_2024.shp"))


## 1. Load random veg plots in each zone
#----------------------------------------------------------------------------
veg<-st_read(paste0(path_out,"Final_outputs/Veg_locations/veg_locations_01_29_25.shp"))%>%
  filter(type=="Random"&
           # also filter records missing coordinate information or that have coordinate errors
           (ot_bnds!=1&wrng_st!=1&iso_rec!=1))%>%
  st_transform("EPSG:4269")%>% #"EPSG:26918"
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         bp="v",
         Year=year(mdy(date)),
         fate=NA)%>% #mark as a random veg location
  #st_transform("EPSG:26918")%>%
  dplyr:: select(id=veg_id,latitude=Lat,longitude=Long,site=Site,bp,Year,fate)

# crop random vegetation points to breeding range
veg2<-st_crop(veg,nests)
