#Data tidying
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
# Fix any coordinate issues with the nest locations
########################################################


## Set file path to data
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"



## 1. Load Nest data
# -------------------------------------------

## Nest fates information - select nest ID and nest fate
fates<-read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2020.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","fate"="UltimateNestFate")%>%
  
  # add flags for each type of coordinate edit we are applying to the original data
  mutate(missing.location.rec=0, # nest fate record is missing a nest location record
         missing.site.info=0, # record is missing site information (state, zone)
         missing.coords=0, # nest location record is missing coordinates
         coord.typo=0, # nest coordinates are not plotting in Mid Atlantic tidal marsh zone
         batch_DecMin_DD_reversed=0, # Coordinates are in Decimal Minute format and lat is in the longitude column (vice vera). Coverted to Decimal Degrees and reversed lat long.
         batch_DecMin_DD=0, # Coordinates are in Decimal Minute format. Coverted to Decimal Degrees.
         batch_dec_addition=0, # Decimal Degree coodinates are missing a decimal. Added decimal.
         batch_dec_addition_reversed=0, # Decimal Degree coodinates are missing a decimal and lat is in the longitude column (vice vera). Added decimal and reversed lat long columns.
         batch_move_DD_to_LatLong=0) # Decimal Degrees were in the UTM northing easting columns. Moved to lat long.
  # manually changed the nest id "id19051" to read "ID19051"



## Nest location information - select nest ID, Site code, Year, Species, Coordinate information
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2001-2020.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site", "Year", "Species",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))
  # adjusted "SAlS" to "SALS" for 1 nest id

nrow(filter(nests,Species=="SALS",Year>=2010))


## 2. Load Site information - compile Site code, site name, and state
#-----------------------------------------------------------------------------------
## some site information from SHARP banding SOP
sites_sop<-read.csv(paste0(dat_path,"Demographic Database/Sites.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("site.code"="Site_Code","Site", "State")

## more sites from Kate
sites_kate<-read.csv(paste0(dat_path,"Demographic Database/Nest_Site_Metadata_kate.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("site.code"="SiteID", "State","Region")

## additional sites' information from Sam A's code
  sites_sam<- data.frame(site.code = c("AT", "BI", "CF", "CL", "DI", "EL", "ER", "HM", "FB", "FS", 
                                    "ID", "JC", "JO", "LU", "MN", "MQ", "MW", "Morris Island", "NC", "NO", 
                                    "OC", "PA", "PB", "PR", "SA", "SC", "SG", "SY", "WI", "NC2", 
                                    "MM", "SP", "FS1", "WA", "MIRI", "STJN", "WOBE"),
                      State = c("NJ", "CT", "MA", "NH", "ME", "ME", "CT", "CT", "ME", "NY",
                                "NY", "RI", "ME", "NH", "NY", "ME", "NJ", "NJ", "NY", "ME", 
                                "NJ", "XX", "ME", "MA", "NY", "ME", "VA", "ME", "VA", "NY", 
                                "ME", "RI", "NY", "CT", "DE", "DE", "DE"),
                      utm.zone = c("18T", "19T", "19T", "19T", "19T", "19T", "18T", "18T", "19T", "18T", 
                                  "18T", "19T", "19T", "19T", "18T", "19T", "18T", "18T", "18T", "19T", 
                                  "18T", NA, "19T", "19T", "18T", "19T", "18T", "19T", "18T", "19T", 
                                  "19T", "19T", "18T", "19T", "18T", "18T", "18T"))

  
## consolidate all the sites into one table
sites<-full_join(sites_sop,sites_sam,by=c("site.code","State"))%>%
  full_join(sites_kate,by=c("site.code","State"))%>%
  # remove any duplicate sites across the different sources
  distinct(site.code,.keep_all = T)


## Assign UTM zones to each site
    # Greater than 72 degrees Longitude is zone 18, less than that is zone 19
    # CT is the only state split between zones 18 and 19, so we will assign CT sites once we fix the coordinates in following sections
    # for now assign zones to sites in states that are completely within 1 zone
sites <- sites%>%
  mutate(utm.zone=case_when(
    # RI, MA, NH, ME are all within UTM zone 19
    is.na(utm.zone) & State %in% c("RI","MA","NH","ME") ~ "19T",
    # NY, NJ,VA, MD, DE are in zone 18
    is.na(utm.zone) & State %in% c("NY","NJ","VA","MD","DE") ~ "18T",
    !(is.na(utm.zone)) ~ utm.zone
  ))





## 3. Address missing records between nest locations, fates, and sites data
# -----------------------------------------------------------------------

# Not all nests in the location records have nest fate records
missing_fates<-nests%>%
  filter(!(id%in%fates$id))
nrow(missing_fates)# number of nests missing fates
  # write this list to file for documentation
if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata.csv"))){
write.csv(missing_fates,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata.csv"), row.names = F)
}


# Not all nests in the fate records have location data 
missing_locations<-fates%>%
  filter(!(id%in%nests$id))
nrow(missing_locations) #(44 nest fates are missing location data)
fates[fates$id%in%missing_locations$id,]$missing.location.rec<-1
# this removes 38 fate records from site WB and 5 records with NA site codes
#fates<-fates%>%
#  filter(id%in%nests$id)



# Site codes missing site information (site name, state, and zone) - see if SHARP folks know where these are
missing_states<- nests%>%
  filter(!(site.code%in%sites$site.code))%>%
  distinct(site.code,.keep_all = T)
  # mark with our editing flag defined above
fates[fates$id%in%missing_states$id,]$missing.site.info<-1

# Can some site names be merged?
site_list<-data.frame(site.code=sort(unique(nests$site.code)))%>%left_join(sites,by="site.code")

if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nest_site_list.csv"))){
write.csv(site_list, paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nest_site_list.csv"),row.names = F)
}




## 4. merge nest fates and locations
# -------------------------------------------
dat1<-left_join(fates,nests,by="id")

# Fix any site name discrepancies (different codes for the same site)
# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek 
unique(dat1$site.code) #- doesn't seem to be an issue in this dataset
# site "WA" is actually "BI", change the one record
dat1<-dat1%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WA","BI",site.code),
    id = ifelse(
      substr(id, start=1,stop=2)=="WA",gsub("WA","BI",id),id)
  )
  




## 5. Spatial formatting
# -------------------------------------------

## A) Identify nest locations missing coordinates
######
# For Nests with location data
      # only keep records with a full set of UTM or lat long coordinate pairs
dat2<-filter(dat1,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))%>%
      # also remove coordinates that are 0 or small values (if a coordinate is less than 30, it's likely a typo)
      filter(if_any(c(Easting,Northing,Lat,Long), ~ .>30))



#Mark records as missing coordinates only if they have location records.
dat1[!(dat1$id%in%dat2$id),]$missing.coords<-1
dat1[dat1$missing.location.rec==1,]$missing.coords<-0

nrow(dat1[dat1$missing.coords==1,]) # 402 nests with location records are missing coordinates




## B) Batch edit nest coordinate information with typos
######
# range of longitude (E-W) for eastern coastline should be -67 (Maine) to -79 (Virginia) in Decimal Degrees
# range of latitude (N-S) should be 45 (Maine) to 36 (Virginia) Decimal Degrees 
# Decimal Degree values (latitude/longitude) outside this range are likely typos



#B-1) Convert Degrees Decimal Minutes records  to Decimal Degrees (DDs) (ONLY the case for NJ sites AT, OC, and MW in 2014 and 2015)

dat3<-dat1%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & missing.location.rec!=1 & 
             # and sites AT, OC, or MW in 2014...
             (site.code %in% c("AT","OC","MW") & Year==2014))%>%
         #convert to DDs - add a decimal after first 2 digits and divide the remaining digits by 60
  mutate(Long=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
         Lat=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
         #then remove the original values from the UTM coordinate columns
         Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_DecMin_DD=1)

    # replace the records with their new edits
dat1<-rbind(dat3,dat1[!(dat1$id%in%dat3$id),])




#B-2) Remove the remaining spatial shift in 2014 NJ sites

#2014 NJ records are still shifted Northeast after DDs conversion, slight systematic error remaining
#look at Sam R's original data files for NJ sites
NJ14<-read.csv(paste0(dat_path,"Demographic Database/SESP 2011-2015.csv"))%>%
  dplyr::select(id=ID,Lat2=LAT,Long2=LONG,Year=YEAR)
NJ14<-read.csv(paste0(dat_path,"Demographic Database/SALS 2011-2015.csv"))%>%
  dplyr::select(id=ident,Lat2=Latitude,Long2=Longitude,Year)%>%
  rbind(NJ14)%>%
  mutate(yr=substr(Year,3,4),
         nest=substr(id,3,(nchar(id)-4)),
         sp=substr(id,(nchar(id)-3),nchar(id)),
         site=substr(id,1,2),
         id=paste0(site,yr,sp,nest))

NJ_test<-left_join(NJ14[NJ14$yr==14,c("id","Lat2","Long2")],dat2[,c("id","Easting","Northing")],by="id")%>%
  left_join(dat1[,c("id","Lat","Long")],by="id")%>%
  mutate(dif_lat=Lat2-Lat,
         dif_long=abs(Long2)-Long)
#take the average difference between original data points and database points.
diff_lat<-mean(NJ_test$dif_lat)
diff_long<-mean(NJ_test$dif_long)

# Sam's data are about 0.001 degrees off from ours and are plotting correctly.

# replace all the NJ 2014 data with Sam's SALS and SESP coordinates. 
NJ14<-NJ14%>%filter(yr=='14')
dat1<- dat1%>%
  left_join(NJ14[,c("id","Lat2","Long2")],by="id")%>%
  mutate(Lat = ifelse(id%in%NJ14$id, Lat2, Lat),
         Long = ifelse(id%in%NJ14$id, Long2, Long),
         #mark which edit was made
         replace_dat= ifelse(id%in%NJ14$id, 1, 0))%>%
  dplyr::select(-c("Lat2","Long2","batch_DecMin_DD"))

#For NJ 2014 non-SESP and SALS nests, apply the average difference to those coordinates
dat1<-dat1%>%
  mutate(Lat2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Lat, NA),
         Long2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Long, NA),
         Lat = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Lat+diff_lat, Lat),
         Long = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Long+diff_long, Long),
         coord_shift= ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), 1, 0))
  



# B-3) repeat Decimal Minute to Decimal Degree conversion for records with lat/long reversed: sites AT, OC, and MW in 2015
dat4<-dat1%>%
    # For records with coordinate information...
    filter(missing.coords!=1 & missing.location.rec!=1 & 
             # and sites AT,OC,MW in 2015...
             (site.code %in% c("AT","OC","MW") & Year==2015))%>%
  mutate(Long=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
         Lat=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
         #then remove the unconverted values from the UTM coordinate columns
         Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_DecMin_DD_reversed=1)
dat1<-rbind(dat4,dat1[!(dat1$id%in%dat4$id),])




# B-4) Add a missing decimal to Decimal Degree coordinates (Checked plotting as UTMs and converting dMin to DD, neither worked) 
## HM, ER, BI, JC, SP, 2014 and HM, BI, ER 2015 are DD missing decimals.
dat5<-dat1%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & missing.location.rec!=1 & 
           # That are also either CT and RI sites HM, ER, BI, JC, or SP in 2014 OR
           # just CT sites HM, BI, ER in 2015...
           ((site.code %in% c("HM","ER","BI","JC","SP") & Year==2014) | (site.code %in% c("HM","BI","ER") & Year==2015)))%>%
           #add a decimal to the easting and northing values to turn into decimal degrees
  mutate(Long=as.numeric(paste0(substr(Easting,1,2),".",
                    substr(Easting,3,length(Easting)))),
         Lat=as.numeric(paste0(substr(Northing,1,2),".",
                                substr(Northing,3,length(Northing)))),
         #then remove the original values from the UTM coordinate columns
         Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_dec_addition=1)

dat1<-rbind(dat5,dat1[!(dat1$id%in%dat5$id),])



# B-5) Move Decimal Degree data in UTM easting/westing columns into long/lat columns
dat6<-dat1%>%
  # if Lat is missing and Easting is using values in Latitude range
  filter(abs(Easting)<46 & abs(Easting)>36 & is.na(Lat))%>%
      # fill in Lat from the Easting column
  mutate(Lat=abs(Easting))%>%
  dplyr::select(id,Lat)
  
dat7<-dat1%>%
  # if Lat is missing and and Northing is using values in Latitude range
  filter(abs(Northing)<46 & abs(Northing)>36 & is.na(Lat))%>%
     # fill in Lat from the NOrthing column
  mutate(Lat=abs(Northing))%>%
  dplyr::select(id,Lat)
  
dat8<-dat1%>%
  # if Long is missing and Easting is using values in Longitude range
  filter(abs(Easting)<79&abs(Easting)>64&is.na(Long))%>%
    # fill in Long from the Easting column
  mutate(Long=Easting)%>%
  dplyr::select(-Lat)
  
dat9<-dat1%>%
  # if Long is missing and Northing is using values in Longitude range
  filter(abs(Northing)<79&abs(Northing)>64&is.na(Long))%>%
    # fill in Long from the Northing column
  mutate(Long=Northing)%>%
  dplyr::select(-Lat)
  
dat10<-rbind(dat8,dat9)%>%
  left_join(rbind(dat6,dat7),by="id")%>%
  #remove lat long from easting northing columns
  mutate(Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_move_DD_to_LatLong = 1)

dat1<-rbind(dat10,dat1[!(dat1$id%in%dat10$id),])
      



# B-6) label the coordinate system and unit for nests as Decimal Degrees or as UTM
dat1<-dat1%>%
           mutate(Coordinate.System=case_when(
           (if_all(c(Lat,Long),~!is.na(.)) & missing.coords!=1) ~"Lat/Long(DD)",
           (if_all(c(Easting,Northing),~!is.na(.)) & missing.coords!=1) ~"UTM(m)",
           (if_all(c(Lat,Long),~!is.na(.)) & if_all(c(Easting,Northing),~!is.na(.))) ~"Lat/Long(DD), UTM(m)"))





## C) Fill in remaining missing UTM zones for CT. Longitude >= 72 DD is zone 18 <72 is zone 19
######
calc_zones <- dat1 %>% 
  filter(!is.na(Long)|Long=="NOT REC")%>%
  mutate(utm.zone=case_when(
    (is.na(utm.zone) & abs(as.numeric(Long))>=72) ~ "18T",
    (is.na(utm.zone) & abs(as.numeric(Long))<72) ~ "19T",
    !(is.na(utm.zone)) ~ utm.zone
  ))%>%
  distinct(site.code,.keep_all = T)%>%
  dplyr::select(site.code,utm.zone)

#Join site information to nest data and merge UTM zone info
dat1<-left_join(dat1,sites,by=c("site.code"))%>%
  left_join(calc_zones,by="site.code")%>%
  mutate(utm.zone=case_when(
    is.na(utm.zone) & !(is.na(utm.zone.y)) ~ utm.zone.y,
    is.na(utm.zone) & !(is.na(utm.zone.x)) ~ utm.zone.x,
    !(is.na(utm.zone)) ~ utm.zone))%>%
  #Waterford Beach CT site, WB, should be UTM zone 18T
  mutate(utm.zone=ifelse(site.code=="WB","18T",utm.zone))%>%
  dplyr::select(-c("utm.zone.x","utm.zone.y","coord.system"))
  



## D) Make sure all longitude values are negative
######
dat1<-dat1%>%
  mutate(Long=ifelse(Long!="NOT REC", -abs(as.numeric(Long)), Long))





## 4. Plot the nest data to look for any remaining coordinate errors
# -----------------------------------------------------------------------------

# get coordinate systems
utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#convert coordinate columns to spatial points
plots_utm18 <- st_as_sf(filter(dat1,utm.zone=="18T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(dat1,utm.zone=="19T"& Coordinate.System=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(dat1,Coordinate.System%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")
plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")



# administrative boundaries, for spatial reference
data(us_states)
ne<-filter(us_states,REGION=="Norteast")

#Marsh boundaries, for spatial reference
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))#uvvr overall mean dataset that I transformed coordinate systems in Arc

#plot nest sites
tm_shape(ne) + tm_borders() +
tm_shape(plots) + tm_dots()

#Mark as error if nest is plotting outside marsh layer (using UVVR since its a rasterized version of national wetland inventory tidal marsh polygons)
plots<-terra::extract(marsh,vect(plots),bind=T)%>%
  sf::st_as_sf()%>%
  mutate(coord.typo=ifelse(
    uvvr_mean_utm18_2=="NaN",1,0
  ))%>%
  dplyr::select(-uvvr_mean_utm18_2)


# QA/QC error flag edits in ArcPro ** find a reproducible way of doing this in R 
if(!file.exists(paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_12_9_22.shp"))){
st_write(plots,
         paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_12_9_22.shp"))
}
# (removed error flags around nest points on marsh border and added error flags to nests plotting at the wrong site) 
error_edits<-st_read(paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_12_9_22.shp"))%>%
  st_drop_geometry()%>%
  dplyr::select(id,crd_typ)
plots<-plots%>%
  left_join(error_edits,by="id")%>%
  mutate(coord.typo=ifelse(!(site.code%in%c("AT","OC","MW") & Year==2014),crd_typ,coord.typo))%>%
  dplyr::select(-crd_typ)




## 5. Write outputs to file
#--------------------------------------------

#nest points shapefile and kml file
output_shp<-dplyr::select(plots,id,coord.typo)%>%
  left_join(dplyr::select(dat1,-c("coord.typo","Lat2","Long2")),by="id")%>%
  distinct(id,.keep_all = T)

st_write(output_shp,
         paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_12_9_22.shp"), delete_layer =T)
st_write(output_shp,
         paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_12_9_22.kml"), delete_layer =T)


#nest coordinate csv
output_csv<-dplyr::select(plots,id,coord.typo)%>%
  st_drop_geometry()%>%
  right_join(dplyr::select(dat1,-c("coord.typo")),by="id")%>%
  distinct(id,.keep_all = T)%>%
  mutate(coord.typo=ifelse(is.na(coord.typo),0,coord.typo))%>%
  ## Add edit note information
  mutate(missing.location.rec2=ifelse(missing.location.rec==1,"Nest fate record is missing nest location record. ",""),
         missing.site.info2=ifelse(missing.site.info==1,"Missing site information. ",""),
         missing.coords2=ifelse(missing.coords==1,"Missing coordinates in nest location record. ",""),
         batch_DecMin_DD_reversed2=ifelse(batch_DecMin_DD_reversed==1,"Batch conversion for NJ sites (OC, AT, MW) in 2015: Degree Decimal Minutes to Decimal Degrees and reversed Lat and Long column values. ",""),
         replace_dat2=ifelse(replace_dat==1,"Batch conversion for NJ sites (OC, AT, MW) in 2014: Replaced coordinates with Sam R's original data. ",""),
         batch_dec_addition2=ifelse(batch_dec_addition==1,"Batch conversion for CT and RI sites (HM, BI, ER, JC, SP) in 2014 and just CT (HM, BI, ER) in 2015: Added missing decimal to Decimal Degrees in Easting/Northin cols. ",""),
         coord_shift2=ifelse(coord_shift==1,paste0("Database coordinates are shifted N. East. Added average difference from Sam R's original SALS and SESP nest records to original database coordinates (", round(Lat2,5), " N, ",round(Long2,5)," W). "),""),
         batch_move_DD_to_LatLong2 = ifelse(batch_move_DD_to_LatLong == 1, "Moved decimal degrees in the Easting Northing columns to Lat Long Columns. ", ""),
         coord.typo2=ifelse(coord.typo==1,"Nest locations still plot outside site area. ",""),
         Notes=paste0(missing.location.rec2,missing.site.info2,missing.coords2,batch_DecMin_DD_reversed2,replace_dat2,batch_dec_addition2, coord_shift2, coord.typo2, batch_move_DD_to_LatLong2),
         Notes=ifelse(Notes=="",Notes,paste0(Notes,"12/9/22 EF."))
  )%>%
  dplyr::select(-c("missing.location.rec2","missing.site.info2","missing.coords2","batch_DecMin_DD_reversed2","replace_dat2","batch_dec_addition2","coord.typo2", "coord_shift2","batch_move_DD_to_LatLong2","Lat2","Long2"))

if(!file.exists(paste0(path_out,"Intermediate_outputs/new_nest_coords_12_9_22.csv"))){
write.csv(output_csv,paste0(path_out,"Intermediate_outputs/new_nest_coords_12_9_22.csv"),row.names = F)
}


#Site table info
output_site<-output_csv%>%
  dplyr::select(site.code,Site,State,utm.zone)%>%
  distinct(Site,.keep_all = T)%>%
  right_join(sites,by=c("site.code","Site","State"))%>%
  mutate(utm.zone=case_when(
         is.na(utm.zone.x) ~ utm.zone.y,
         !(is.na(utm.zone.x))~utm.zone.x))%>%
  dplyr::select(-c("utm.zone.x","utm.zone.y"))

write.csv(output_site,paste0(path_out,"Final_outputs/compiled_site_table_12_9_22.csv"),row.names = F)


# error tally
t<-summarise(output_csv,
             missing.location.rec=sum(missing.location.rec),
             missing.site.info=sum(missing.site.info),
             missing.coords=sum(missing.coords),
             batch_DecMin_DD_reversed=sum(batch_DecMin_DD_reversed),
             replace_dat=sum(replace_dat),
             batch_dec_addition=sum(batch_dec_addition),
             no_replacement=sum(no_replacement),
             batch_move_DD_to_LatLong = sum(batch_move_DD_to_LatLong ),
             coord.typo=sum(coord.typo,na.rm=T))

