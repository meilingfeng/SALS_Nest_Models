library(tidyverse)
library(lubridate)
library(stringi)
#spatial analysis
library(sf)
library(sp)
library(terra)
library(tmap)
# spatial datasets
library(spData)
library(geodata)
library(rnaturalearth)


########################################################
# Fix any coordinate issues with the vegetation plots
########################################################



## Set file path to data
# -------------------------------------------
dat_path<-"C:/Users/10788/Desktop/SaltMarsh/Data/"
path_out<-"C:/Users/10788/Desktop/SaltMarsh/Outputs/"



## 1. Load data
# -------------------------------------------

# Cleaned nest fate data
fates<- read.csv(paste0(path_out,"Intermediate_Outputs/new_nest_coords_12_9_22.csv"))


# Original Nest location data
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2001-2020.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site", "Year", "Species",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))


# vegetation plots
veg<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2020.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("veg.id"="VegPointID","type"="PointType","site.code"="Site","date"="SurveyDate",
                "coord.system"="Coordinate.System","utm.zone"="UTM.Zone","Easting","Northing","Lat","Long")%>%
  mutate(type=case_when(
    grepl("R",type)~"Random",
    grepl("N",type)~"Nest"),
         Year = as.numeric(substr(date,nchar(date)-3,nchar(date))))%>%
  distinct(.keep_all = T)


## Compare similar records between vegetation data and nest data
similar<-veg[veg$veg.id%in%nests$id,] # what nest records in the nest location data have vegetation survey data?
nrow(similar)
#5061 shared plots between veg and nest data

nrow(veg[veg$type=="Nest",])
#5075 total nest plots in the veg data 
# (about 14 extra nest records in the veg data)
extra<-veg[!(veg$veg.id%in%nests$id)&veg$type%in%c("Nest","NEST"),]

nrow(veg[veg$type=="Random",])
#5930 random vegeation plots not associated with nests






## 2. Site information - select Site code, site name, and state
# -------------------------------------------------------------
sites <- fates%>%
  dplyr::select(site.code,Site,State,utm.zone)%>%
  distinct(.keep_all = T)%>%
  #For site "WB", assign it to woodland beach in DE over Waterford beach in CT, records are plotting in DE
  filter(!(site.code=="WB"&State=="CT"))
# fill in site info for sub sites found in veg data 
veg_sites<-unique(veg[!(veg$site.code%in%fates$site.code),"site.code"])
veg_sites_state<-c("CT","NY","NY","NY","NY","NY","NY","NY","NY","NY","NY",NA,"ME")
veg_sites_utm<-c("19T","18T","18T","18T","18T","18T","18T","18T","18T","18T","18T",NA,"19T")
veg_sites<-data.frame(site.code=veg_sites,
                         Site=rep(NA,length(veg_sites)),
                         State=veg_sites_state,
                         utm.zone=veg_sites_utm)
sites<-rbind(sites,veg_sites)


# Fix site name discrepancies
# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek - doesn't seem to be an issue in this dataset
unique(veg$site.code)
# WA site is actually BI change the one record
veg<-veg%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WA","BI",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="WA",gsub("WA","BI",veg.id),veg.id)
  )



## 3. Plot Data (mark any typos and leave correct data untouched)
# -------------------------------------------------------------------
# Mark records missing coordinates
veg2<-filter(veg,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))%>%
  # also remove coordinates that are 0 or small values (less than 10, likely typos)
  filter(if_any(c(Easting,Northing,Lat,Long), ~ .>10))

veg$missing.coords<-0
veg[!(veg$veg.id%in%veg2$veg.id),]$missing.coords<-1

nrow(veg[veg$missing.coords==1,])
#3780 missing coordates
ggplot(veg,aes(x=as.factor(missing.coords),color=type))+
  geom_histogram(stat = "count",aes(fill=type))
#majority of missing coordinates are in nest data
nrow(veg[veg$missing.coords==1&veg$type=="Random",])
# 316 random plots are missing coords
rand_missing<-veg[veg$missing.coords==1&veg$type=="Random",]


# Merge utm zones to veg data
veg<-left_join(dplyr::select(veg,-utm.zone),sites,by="site.code")


# label the coordinate system and unit for all these nests as Decimal Degrees and all other nests as UTM
veg<-veg%>%
  mutate(Coordinate.System=case_when(
    (if_all(c(Lat,Long),~!is.na(.)) & missing.coords!=1) ~"Lat/Long(DD)",
    (if_all(c(Easting,Northing),~!is.na(.)) & missing.coords!=1) ~"UTM(m)",
    (if_all(c(Lat,Long),~!is.na(.)) & if_all(c(Easting,Northing),~!is.na(.))) ~"Lat/Long(DD), UTM(m)"))


# get coordinate systems and convert coordinates to spatial data
utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#convert to spatial points
plots_utm18 <- st_as_sf(filter(veg,utm.zone=="18T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(veg,utm.zone=="19T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(veg,Coordinate.System%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")
plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")



#admin boundaries
data(us_states)
ne<-filter(us_states,REGION=="Norteast")

#Marsh boundaries
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))

#plot nest sites
tm_shape(ne) + tm_borders() +
  tm_shape(plots) + tm_dots()


plots<-terra::extract(marsh,vect(plots),bind=T)%>%
  sf::st_as_sf()%>%
  mutate(coord.typo=ifelse(
    uvvr_mean_utm18_2=="NaN",1,0
  ))%>%
  dplyr::select(-uvvr_mean_utm18_2)

# verify typos in arcpro
output_shp<-dplyr::select(plots,veg.id,coord.typo)%>%
  distinct(veg.id,.keep_all = T)
if(!(file.exists(paste0(path_out,"Intermediate_Outputs/Nest_locations/veg_locations_error_edits_12_26_22.shp")))){
st_write(output_shp,
         paste0(path_out,"Intermediate_Outputs/Nest_locations/veg_locations_error_edits_12_26_22.shp"), delete_layer =T)
}

# Mark records with typos and extract them for edits
veg_edits<-st_read(paste0(path_out,"Intermediate_Outputs/Nest_locations/veg_locations_error_edits_12_26_22.shp"))%>%
  st_drop_geometry()%>%
  right_join(veg,by=c("veg_id"="veg.id"))%>%
  
  #extract typos and missing coordinates to work with
  filter(coord_typo==1)%>%
  rename(veg.id=veg_id)%>%
  dplyr::select(-coord_typo)
veg_edits<-rbind(veg_edits,filter(veg,missing.coords==1))%>%
  distinct(.keep_all = T)





## 4. Fill in missing nest coordinate data from nest dataset
#---------------------------------------------------------

#first see if shared nests have the same coordinate info for those that do have info in the veg data
test<-veg_edits%>%
  filter(veg.id%in%nests$id)%>%
  dplyr::select(veg.id,site.code2=site.code,Easting2=Easting,Northing2=Northing,Lat2=Lat,Long2=Long, utm.zone)%>%
  left_join(nests[,c("id","site.code","Easting","Northing","Lat","Long")],by=c("veg.id"="id"))

#nests that have data in nest doc but not veg
test3<-test%>%
  filter((is.na(Easting2)&is.na(Northing2)&is.na(Lat2)&is.na(Long2))&((!is.na(Lat)&!is.na(Long))|(!is.na(Easting)&!is.na(Northing))))
nrow(test3)

#nests that have data in both but differ
test2<-test%>%
  filter((!is.na(Easting2)|!is.na(Long2)) & (Easting2!=Easting | Long2 != Long))
nrow(test2)

utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"


#convert to spatial points
plots_veg18 <- st_as_sf(filter(test2,!is.na(Easting2)&utm.zone=="18T"), coords = c("Easting2", "Northing2"), crs = utm18)%>%
  st_transform(nad)
plots_veg19 <- st_as_sf(filter(test2,!is.na(Easting2)&utm.zone=="19T"), coords = c("Easting2", "Northing2"), crs = utm19)%>%
  st_transform(nad)
plots_veg<-rbind(plots_veg18,plots_veg19)

plots_nest18 <- st_as_sf(filter(test2,!is.na(Easting2)&utm.zone=="18T"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)
plots_nest19 <- st_as_sf(filter(test2,!is.na(Easting2)&utm.zone=="19T"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)
plots_nest<-rbind(plots_nest18,plots_nest19)

#admin boundaries
data(us_states)
ne<-filter(us_states,REGION=="Norteast")



#plot nest sites
tm_shape(ne) + tm_borders() +
  tm_shape(plots_veg) + tm_dots(col = "red")+
  tm_shape(plots_nest) + tm_dots(col = "blue")
# None are nests included in the batch edits, most seem to be really close or typos either way
#write.csv(test2,paste0(dat_path,"veg_nest_coord_differences.csv"))


# look at all veg coordinates (including those that were not flagged as typos) that differ from nest data
test<-veg%>%
  filter(veg.id%in%nests$id)%>%
  dplyr::select(veg.id,site.code2=site.code,Easting2=Easting,Northing2=Northing,Lat2=Lat,Long2=Long, utm.zone)%>%
  left_join(nests[,c("id","site.code","Easting","Northing","Lat","Long")],by=c("veg.id"="id"))

#nests that have data in both but differ
test2<-test%>%
  filter((!is.na(Easting2)|!is.na(Long2)) & (Easting2!=Easting | Long2 != Long))%>%
  #add typo status in veg and nest data
  left_join(fates%>%dplyr::select(id,coord.typo,Notes),by=c("veg.id"="id"))%>%
  mutate(coord.typo2=ifelse(veg.id%in%veg_edits$veg.id,1,0))
nrow(test2)


test2<-test2%>%
  # for nests flagged as typos in the nest data, but not in the veg data, use the veg coordinates (likely a typo in the nest data)
  # for nests that are not typos in either dataset (coordinates are very close for these) just keep the veg coordinates
  mutate(use.veg=ifelse(coord.typo2==0 | site.code%in%c("SJ","WB"), 1,0),
  # for nests flagged as typos in the veg data, but not in the nest data, use the nest coordinates (likely a typo in the veg data - except SJ and WB sites, these are just in lat long format (missing a decimal) in the veg data as opposed to UTM in the nest data)
         use.nest=ifelse((coord.typo==0|is.na(coord.typo)) & coord.typo2==1 & !(site.code%in%c("SJ","WB")), 1,0))
sum(test2$use.veg)
# 84 (4 resolve typos in the nest data)
sum(test2$use.nest)
# 18


# Fill in all veg data with missing or typo coordinates with corresponding nest data coordinates (except for those with correct veg data over nest data)
# applying to veg edits automatically retains veg coordinates not flagged as typos
veg_edits2<-dplyr::select(veg_edits,-c("site.code","Easting","Northing","Lat","Long"))%>%
  # filter to nests in nest data and that were not flagged to use veg data (these are just the WB and SJ sites)
  filter(veg.id%in%nests$id & !(veg.id%in%test2[test2$use.veg==1,]$veg.id))%>%
  left_join(dplyr::select(nests,id,site.code,Easting,Northing,Lat,Long),by=c("veg.id"="id"))%>%
  mutate(nest.data.add=1)
veg_edits<-rbind(mutate(filter(veg_edits,!(veg.id%in%veg_edits2$veg.id)),nest.data.add=0),veg_edits2)

#fix the 4 typos in the nest data by replacing with correct coordinates in veg data
fates2<-dplyr::select(fates,-c("site.code","Easting","Northing","Lat","Long"))%>%
  filter(id%in%test2[test2$use.veg==1,]$veg.id & coord.typo==1)%>%
  left_join(dplyr::select(veg,veg.id,site.code,Easting,Northing,Lat,Long),by=c("id"="veg.id"))%>%
  mutate(veg.data.replace=1)
fates<-rbind(mutate(filter(fates,!(id%in%fates2$id)),veg.data.replace=0),fates2)%>%
# add new note of changes (mostly a repeated number deleted)
        mutate(Notes=ifelse(veg.data.replace==1,"Replaced with correct coordinates in vegetation dataset. 01/03/23 EF. ",Notes),
               coord.typo=ifelse(veg.data.replace==1,0,coord.typo))


# Re-mark records missing coordinate info
# select records with a full set of UTM or lat long coordinate pairs
veg_edits3<-filter(veg_edits,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))%>%
  # also remove coordinates that are 0 or small values (less than 10, likely typos)
  filter(if_any(c(Easting,Northing,Lat,Long), ~ .>10))

# Mark records missing coordinates
veg_edits$missing.coords<-0
veg_edits[!(veg_edits$veg.id%in%veg_edits3$veg.id),]$missing.coords<-1

nrow(veg_edits[veg_edits$missing.coords==1,])
#377 missing coordates
ggplot(veg_edits,aes(x=as.factor(missing.coords),color=type))+
  geom_histogram(stat = "count",aes(fill=type))
#majority of missing coordinates are in nest data
nrow(veg_edits[veg_edits$missing.coords==1&veg_edits$type=="Random",])
# 316 random plots are missing coords



## 5. Apply batch coordinate edits (Same batch edits for nests appear to mostly apply to the random veg sites)
#----------------------------------------------------------------------------------------------------------------

# 5-1) Degrees Decimal Minutes records (ONLY the case for NJ sites AT, OC, and MW in 2014 and 2015)

veg_edits4<-veg_edits%>%
  # For records with coordinate information...
  filter(missing.coords!=1 &
           # and sites AT, OC, or MW in 2014...
           (site.code %in% c("AT","OC","MW") & Year==2014))%>%
  #convert to DDs
  mutate(Long=(as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60)),
         Lat=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
         #then remove the original values from the UTM coordinate columns
         Easting=NA,
         Northing=NA)


veg_edits<-rbind(veg_edits4,veg_edits[!(veg_edits$veg.id%in%veg_edits4$veg.id),])

#See if rand veg plots are also plotting at an offset
# verify in arcpro
#output_shp<-dplyr::select(veg_edits4,veg.id,Lat,Long)%>%
#  distinct(veg.id,.keep_all = T)%>%
#  st_as_sf(coords=c("Long","Lat"),crs=nad)%>%
#  st_transform("EPSG:26918")

#st_write(output_shp,
#         paste0(dat_path,"Nest Locations/veg_NJ_14_edits.shp"), delete_layer =T)
#Yes they are
diff_lat<--0.001570855 #From nest data cleaning
diff_long<-0.001876374

# replace all the NJ 2014 SESP SALS data with fixed nest coordinates. 
NJ14<-fates%>%filter(Year==2014 & site.code %in% c("AT","OC","MW"), Species %in% c("SESP", "SALS"))%>%
  rename(Lat2=Lat,Long2=Long)
veg_edits<- veg_edits%>%
  left_join(NJ14[,c("id","Lat2","Long2")],by=c("veg.id"="id"))%>%
  mutate(Lat = ifelse(veg.id%in%NJ14$id, Lat2, Lat),
         Long = ifelse(veg.id%in%NJ14$id, Long2, Long),
         #mark which edit was made
         replace_dat= ifelse(veg.id%in%NJ14$id, 1, 0))%>%
  dplyr::select(-c("Lat2","Long2"))

#For NJ 2014 non-SESP and SALS nests, apply the average difference to those coordinates
veg_edits<-veg_edits%>%
  mutate(Lat2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(veg.id%in%NJ14$id), Lat, NA),#this maintains the original coordinates to add to the notes at the end
         Long2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(veg.id%in%NJ14$id), Long, NA),
         Lat = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(veg.id%in%NJ14$id), Lat+diff_lat, Lat),
         Long = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(veg.id%in%NJ14$id), Long+diff_long, Long),
         coord_shift= ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(veg.id%in%NJ14$id), 1, 0))



# repeat DD conversion for records with lat/long reversed: At, OC, and MW in 2015
veg_edits5<-veg_edits%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & 
           # and sites AT,OC,MW in 2015...
           (site.code %in% c("AT","OC","MW") & Year==2015))%>%
  mutate(Long=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
         Lat=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
         #then remove the unconverted values from the UTM coordinate columns
         Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_DecMin_DD_reversed=1)
veg_edits<-rbind(veg_edits5,mutate(veg_edits[!(veg_edits$veg.id%in%veg_edits5$veg.id),],batch_DecMin_DD_reversed=0))


# 5-2) Batch conversion of UTM coords that are actually DD missing decimals
## HM, ER, BI, JC, SP, 2014 and HM, BI, ER 2015 are DD missing decimals in nest data, is NOT the case for rand veg data.
# just replace the nest records with fixed nest coordinates. 

CT_RI<-fates%>%
  filter(missing.coords!=1 & missing.location.rec!=1 & 
           ((site.code %in% c("HM","ER","BI","JC","SP") & Year==2014) | (site.code %in% c("HM","BI","ER") & Year==2015)))%>%
  rename(Lat3=Lat,Long3=Long)

veg_edits<- veg_edits%>%
  left_join(CT_RI[,c("id","Lat3","Long3")],by=c("veg.id"="id"))%>%
  mutate(Lat = ifelse(veg.id%in%CT_RI$id, Lat3, Lat),
         Long = ifelse(veg.id%in%CT_RI$id, Long3, Long),
         #mark which edit was made
         batch_dec_addition= ifelse(veg.id%in%CT_RI$id, 1, 0))%>%
  dplyr::select(-c("Lat3","Long3"))


# MR 2019, SJ 2019 , WB 2019 for the random veg plots also appear to be missing decimals (some also have lat long reversed)
veg_edits6<-veg_edits%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & type=="Random"&
           ((site.code %in% c("MR","WB","SJ") & Year==2019)))%>%
  #add a decimal to the easting and northing values to turn into decimal degrees
  mutate(Easting=as.numeric(paste0(substr(Easting,1,2),".",
                                substr(Easting,3,nchar(Easting)))),
         Northing=as.numeric(paste0(substr(Northing,1,2),".",
                               substr(Northing,3,nchar(Northing)))),
         Lat=ifelse(abs(Easting)>65,Northing,Easting),
         Long=ifelse(abs(Easting)>65,-1*abs(Easting),-1*abs(Northing)),
         #then remove the original values from the UTM coordinate columns
         Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_dec_addition=1)

veg_edits<-rbind(veg_edits6,veg_edits[!(veg_edits$veg.id%in%veg_edits6$veg.id),])



# 5-3) Move DD data in easting/westing columns into long/lat columns (only applies to nest observations n =49)
dat6<-veg_edits%>%
  # if Lat is missing and Easting is using values in Latitude range
  filter(abs(Easting)<46 & abs(Easting)>36 & is.na(Lat))%>%
  # fill in Lat from the Easting column
  mutate(Lat=abs(Easting))%>%
  dplyr::select(veg.id,Lat)

dat7<-veg_edits%>%
  # if Lat is missing and and Northing is using values in Latitude range
  filter(abs(Northing)<46 & abs(Northing)>36 & is.na(Lat))%>%
  # fill in Lat from the NOrthing column
  mutate(Lat=abs(Northing))%>%
  dplyr::select(veg.id,Lat)

dat8<-veg_edits%>%
  # if Long is missing and Easting is using values in Longitude range
  filter(abs(Easting)<79&abs(Easting)>64&is.na(Long))%>%
  # fill in Long from the Easting column
  mutate(Long=Easting)%>%
  dplyr::select(-Lat)

dat9<-veg_edits%>%
  # if Long is missing and Northing is using values in Longitude range
  filter(abs(Northing)<79&abs(Northing)>64&is.na(Long))%>%
  # fill in Long from the Northing column
  mutate(Long=Northing)%>%
  dplyr::select(-Lat)

dat10<-rbind(dat8,dat9)%>%
  left_join(rbind(dat6,dat7),by="veg.id")%>%
  #remove lat long from easting northing columns
  mutate(Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_move_DD_to_LatLong = 1)

veg_edits<-rbind(dat10,mutate(veg_edits[!(veg_edits$veg.id%in%dat10$veg.id),],batch_move_DD_to_LatLong=0))

# 5-4) label the coordinate system and unit for all these nests as Decimal Degrees and all other nests as UTM
veg_edits<-veg_edits%>%
  mutate(Coordinate.System=case_when(
    (if_all(c(Lat,Long),~!is.na(.)) & missing.coords!=1) ~"Lat/Long(DD)",
    (if_all(c(Easting,Northing),~!is.na(.)) & missing.coords!=1) ~"UTM(m)",
    (if_all(c(Lat,Long),~!is.na(.)) & if_all(c(Easting,Northing),~!is.na(.))) ~"Lat/Long(DD), UTM(m)"))

veg<-rbind(filter(mutate(veg,
                         Lat2=0,Long2=0,batch_move_DD_to_LatLong=0,
                         batch_dec_addition=0,batch_DecMin_DD_reversed=0,
                         replace_dat=0,coord_shift=0,nest.data.add=0
                         ),
                  !(veg.id%in%veg_edits$veg.id)
                  ),
            veg_edits)%>%
  distinct(.keep_all = T)



#Make sure all longitude values are negative
veg<-veg%>%
  mutate(Long=ifelse(Long!="NOT REC", -abs(as.numeric(Long)), Long))%>%
  dplyr::select(-coord.system)



## 6. Plot Data
# -------------------------------------------

# get coordinate systems and convert coordinates to spatial data
utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#convert to spatial points
plots_utm18 <- st_as_sf(filter(veg,utm.zone=="18T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(veg,utm.zone=="19T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(veg,Coordinate.System%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")
plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")



#admin boundaries
data(us_states)
ne<-filter(us_states,REGION=="Norteast")

#Marsh boundaries
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))

#plot nest sites
tm_shape(ne) + tm_borders() +
  tm_shape(plots) + tm_dots()

# mark if plotting outside marsh bounds
# select only the records that were just edited since we already marked them previously
plots2<-filter(plots,nest.data.add==1|batch_move_DD_to_LatLong==1|batch_dec_addition==1|batch_DecMin_DD_reversed==1|replace_dat==1|coord_shift==1)
plots3<-terra::extract(marsh,vect(plots2),bind=T)%>%
  sf::st_as_sf()%>%
  mutate(coord.typo=ifelse(
    uvvr_mean_utm18_2=="NaN",1,0
  ))%>%
  dplyr::select(-uvvr_mean_utm18_2)


#nest points shp
output_shp<-dplyr::select(plots3,veg.id,coord.typo)%>%
  left_join(veg,by="veg.id")%>%
  distinct(veg.id,.keep_all = T)
if(!(file.exists(paste0(path_out,"Intermediate_Outputs/Nest Locations/veg_edit_locations_12_26_22.shp")))){
st_write(plots3,
         paste0(path_out,"Intermediate_Outputs/Nest Locations/veg_edit_locations_12_26_22.shp"), delete_layer =T)
}

#read back file with typo adjustments in arcpro
plots2<-st_read(paste0(dat_path,"Nest Locations/veg_edit_locations_12_26_22.shp"))%>%
  st_transform(crs(plots))%>%
  dplyr::select(coord.typo=crd_typ, veg.id=veg_id)%>%
  left_join(st_drop_geometry(plots2),by="veg.id")



plots<-plots%>%
  filter(!(veg.id%in%plots2$veg.id))%>%
  # all remaining records in veg edits have typos, mark with 1
  mutate(coord.typo=ifelse(veg.id%in%veg_edits$veg.id,1,0))%>%
  # then add records that were changed in veg edits and had their typo flags adjusted
  rbind(plots2)

sum(plots$coord.typo)
#213 typos in total remain
982-213
#769 resolved
nrow(veg)-nrow(plots)
#401 still missing coordinate data


veg_final<-veg%>%
  left_join(st_drop_geometry(dplyr::select(plots,veg.id,coord.typo)),by="veg.id")%>%
  distinct(.keep_all = T)
veg_final[is.na(veg_final$coord.typo),]$coord.typo<-0
# add notes
veg_final<- mutate(veg_final,missing.coords2=ifelse(missing.coords==1,"Missing coordinate information. ",""),
        batch_DecMin_DD_reversed2=ifelse(batch_DecMin_DD_reversed==1,"Batch conversion for NJ sites (OC, AT, MW) in 2015: Degree Decimal Minutes to Decimal Degrees and reversed Lat and Long column values. ",""),
        replace_dat2=ifelse(replace_dat==1,"Batch conversion for NJ sites (OC, AT, MW) in 2014: Replaced coordinates with Sam R's original data. ",""),
        batch_dec_addition2=ifelse(batch_dec_addition==1,"Batch conversion for CT and RI sites (HM, BI, ER, JC, SP) in 2014 and just CT (HM, BI, ER) in 2015: Added missing decimal to Decimal Degrees in Easting/Northin cols. ",""),
        coord_shift2=ifelse(coord_shift==1,paste0("Database coordinates are shifted N. East. Added average difference from Sam R's original SALS and SESP nest records to original database coordinates (", round(Lat2,5), " N, ",round(Long2,5)," W). "),""),
        batch_move_DD_to_LatLong2 = ifelse(batch_move_DD_to_LatLong == 1, "Moved decimal degrees in the Easting Northing columns to Lat Long Columns. ", ""),
        coord.typo2=ifelse(coord.typo==1,"Nest locations still plot outside site area. ",""),
        nest.data.add2=ifelse(nest.data.add==1,"Replaced with coordinates in nest database. ",""),
        Notes_Emily=paste0(missing.coords2,batch_DecMin_DD_reversed2,replace_dat2,batch_dec_addition2, coord_shift2, coord.typo2, batch_move_DD_to_LatLong2,nest.data.add2),
        Notes_Emily=ifelse(Notes_Emily=="",Notes_Emily,paste0(Notes_Emily,"12/9/22 EF."))
  )%>%
  dplyr::select(-c("missing.coords2","batch_DecMin_DD_reversed2","replace_dat2","batch_dec_addition2","coord.typo2", "coord_shift2","batch_move_DD_to_LatLong2","Lat2","Long2","nest.data.add2"))


plots_final<-dplyr::select(plots,veg.id,date)%>%
  distinct(.keep_all = T)




## 7.  Make final shp and csv outputs
#------------------------------------------------------------

# read in original veg data variables
veg_orig<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2020.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select(-c("PointType","Site","UTM.Zone","Easting","Northing","Lat","Long","Coordinate.System"))
#join to new coordinate info
output_csv<-dplyr::select(veg_final,VegPointID=veg.id,coord.typo,SurveyDate=date,PointType=type, Site_Name=Site, Site=site.code,UTM.Zone=utm.zone,
                          Easting,Northing,Lat,Long,Notes_Emily)%>%
  right_join(veg_orig,by=c("VegPointID","SurveyDate"))%>%
  distinct(.keep_all = T)


#join to shapefile
output_shp<-plots_final%>%
  left_join(output_csv,by=c("veg.id"="VegPointID","date"="SurveyDate"))
#write files
  #veg data files
st_write(output_shp,
         paste0(path_out,"Final_outputs/Veg_Locations/veg_locations_12_29_22.shp"), delete_layer =T)
st_write(output_shp,
         paste0(path_out,"Final_outputs/Veg_locations/veg_locations_12_29_22.kml"), delete_layer =T)

write.csv(output_csv,paste0(path_out,"Final_outputs/new_veg_coords_12_29_22.csv"),row.names = F)
write.csv(veg_final,paste0(path_out,"Final_outputs/new_veg_coords_12_29_22_wEditflags.csv"),row.names = F)

  # adjusted nest data files
write.csv(fates,paste0(path_out,"Final_outputs/new_nest_coords_01_3_23.csv"),row.names = F)

    #convert to spatial points
plots_utm18 <- st_as_sf(filter(fates,utm.zone=="18T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(fates,utm.zone=="19T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(fates,Coordinate.System%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")

plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")

st_write(plots,
         paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"), delete_layer =T)
st_write(plots,
         paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.kml"), delete_layer =T)
 
