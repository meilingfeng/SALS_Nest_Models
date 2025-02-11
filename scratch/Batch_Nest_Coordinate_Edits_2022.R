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
fates<-read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","fate"="UltimateNestFate")%>%
  
  # add flags for each type of coordinate edit we are applying to the original data
  mutate(missing.location.rec=0, # nest fate record is missing a nest location record
         missing.site.info=0, # record is missing site information (state, zone)
         missing.coords=0, # nest location record is missing coordinates
         out.bounds=0, # nest coordinates are not plotting in Mid Atlantic tidal marsh zone
         wrong.site=0, # nest is plotting in a site different than the one in its record
         iso.rec=0, # nest is in tidal marsh but plotting further than 500m away from any other nest record
         batch_DecMin_DD_reversed=0, # Coordinates are in Decimal Minute format and lat is in the longitude column (vice vera). Coverted to Decimal Degrees and reversed lat long.
         batch_DecMin_DD=0, # Coordinates are in Decimal Minute format. Coverted to Decimal Degrees.
         batch_dec_addition=0, # Decimal Degree coodinates are missing a decimal. Added decimal.
         batch_dec_addition_reversed=0, # Decimal Degree coodinates are missing a decimal and lat is in the longitude column (vice vera). Added decimal and reversed lat long columns.
         batch_move_DD_to_LatLong=0) # Decimal Degrees were in the UTM northing easting columns. Moved to lat long.
# manually changed the nest id "id19051" to read "ID19051"



## Nest location information - select nest ID, Site code, Year, Species, Coordinate information
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2002-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site", "Year", "Species", "Site"="SiteName",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))
# adjusted "SAlS" to "SALS" for 1 nest id



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

#add different code for woodland beach
sites<-sites%>%
  rbind(data.frame(site.code="WO",Site="Woodland Beach",State="DE",utm.zone="18T",Region=NA))



## 3. Address missing records between nest locations, fates, and sites data
# -----------------------------------------------------------------------
# Not all nests in the location records have nest fate records
missing_fates<-nests%>%
  filter(!(id%in%fates$id))
nrow(missing_fates)# number of nests missing fates, 492
# write this list to file for documentation
if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata_2002_2024.csv"))){
  write.csv(missing_fates,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata_@002_2024.csv"), row.names = F)
}


# Not all nests in the fate records have location data 
missing_locations<-fates%>%
  filter(!(id%in%nests$id))
nrow(missing_locations) #(44 nest fates are missing location data) (now only 8 with the 2024 update!)
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
#site_list<-data.frame(site.code=sort(unique(nests$site.code)))%>%left_join(sites,by="site.code")

#if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nest_site_list.csv"))){
#write.csv(site_list, paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nest_site_list.csv"),row.names = F)
#}




## 4. merge nest fates and locations
# -------------------------------------------
dat1<-full_join(fates,nests,by="id")

# Fix any site name discrepancies (different codes for the same site)
# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek 
# WO WB and WOBE are all Woodland Beach, not for 2024 data
sort(unique(dat1$site.code)) #- doesn't seem to be an issue in this dataset
# site "WA" is actually "BI", change the one record
dat1<-dat1%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WA","BI",site.code),
    id = ifelse(
      substr(id, start=1,stop=2)=="WA",gsub("WA","BI",id),id)
  )%>%
  #  mutate(site.code = ifelse(
  #    substr(site.code, start=1,stop=2)=="WO","WB",site.code),
  #    id = ifelse(
  #      substr(id, start=1,stop=2)=="WO",gsub("WO","WB",id),id)
  #  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=3)=="LtR","LR2",site.code),
    id = ifelse(
      substr(id, start=1,stop=3)=="LtR",gsub("LtR","LR",id),id)
  )%>%
  #remove any sites listed as "RD" or with numbers in their site code- these are rapid demo sites and we want intensive demo
  filter(!(site.code=="RD" | grepl("^[[:digit:]]+",site.code) | substring(id,1,2)=="RD"))
#FR, BH, and BF are all clustered at one site in CT... I'm guessing a naming convention but not addressing here since its all pre SHARP data




## 5. Spatial formatting
# -------------------------------------------

## A) Identify nest locations missing coordinates
######
# filter for records with a full set of UTM or lat long coordinate pairs
dat2<-filter(dat1,((if_all(c(Easting,Northing), ~ !is.na(.)))|(if_all(c(Lat,Long), ~ !is.na(.)))))%>%
  # also remove coordinates that are 0 or small values (if a coordinate is less than 30, it's likely a typo)
  filter((if_all(c(Easting,Northing), ~ .>30))|(if_all(c(Lat,Long), ~ abs(.)>30)))



#Mark records missing coordinates.
dat1[!(dat1$id%in%dat2$id),]$missing.coords<-1 # ids not in records with valid coordinates
dat1[dat1$id%in%dat2$id,]$missing.coords<-0 # ids are in records with valid coordinates

nrow(dat1[dat1$missing.coords==1,]) # 463, now 481 nests with location records are missing coordinates

## B) Batch edit nest coordinate information with typos (Don't need to do these for the new nest data in 2024, it has these updates in it)
######
# range of longitude (E-W) for eastern coastline should be -67 (Maine) to -79 (Virginia) in Decimal Degrees
# range of latitude (N-S) should be 45 (Maine) to 36 (Virginia) Decimal Degrees 
# Decimal Degree values (latitude/longitude) outside this range are likely typos



#B-1) Convert Degrees Decimal Minutes records to Decimal Degrees (DDs) (ONLY the case for NJ sites AT, OC, and MW in 2014 and 2015)

#dat3<-dat1%>%
# For records with coordinate information...
#  filter(missing.coords!=1 & missing.location.rec!=1 & 
# and sites AT, OC, or MW in 2014...
#             (site.code %in% c("AT","OC","MW") & Year==2014))%>%
#convert to DDs - add a decimal after first 2 digits and divide the remaining digits by 60
#  mutate(Long=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
#         Lat=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
#then remove the original values from the UTM coordinate columns
#         Easting=NA,
#         Northing=NA,
#mark which edits were made
#         batch_DecMin_DD=1)

# replace the records with their new edits
#dat1<-rbind(dat3,dat1[!(dat1$id%in%dat3$id),])




#B-2) Remove the remaining spatial shift in 2014 NJ sites

#2014 NJ records are still shifted Northeast after DDs conversion, slight systematic error remaining
#look at Sam R's original data files for NJ sites
#NJ14<-read.csv(paste0(dat_path,"Demographic Database/SESP 2011-2015.csv"))%>%
#  dplyr::select(id=ID,Lat2=LAT,Long2=LONG,Year=YEAR)
#NJ14<-read.csv(paste0(dat_path,"Demographic Database/SALS 2011-2015.csv"))%>%
#  dplyr::select(id=ident,Lat2=Latitude,Long2=Longitude,Year)%>%
#  rbind(NJ14)%>%
#  mutate(yr=substr(Year,3,4),
#         nest=substr(id,3,(nchar(id)-4)),
#         sp=substr(id,(nchar(id)-3),nchar(id)),
#         site=substr(id,1,2),
#         id=paste0(site,yr,sp,nest))%>%
#  filter(Year==2014)%>%
#  dplyr::select(id,Lat2,Long2)

#NJ_test<-left_join(NJ14,dat2[,c("id","Easting","Northing")],by="id")%>%
#  left_join(dat1[,c("id","Lat","Long")],by="id")%>%
#  mutate(dif_lat=Lat2-Lat,
#         dif_long=abs(Long2)-Long)
#take the average difference between original data points and database points.
#diff_lat<-mean(NJ_test$dif_lat)
#diff_long<-mean(NJ_test$dif_long)

# Sam's data are about 0.001 degrees off from ours and are plotting correctly.

# replace all the NJ 2014 data with Sam's SALS and SESP coordinates. 
#dat1<- left_join(dat1,NJ14,by="id")%>%
#  mutate(Lat = ifelse(id%in%NJ14$id, Lat2, Lat),
#         Long = ifelse(id%in%NJ14$id, Long2, Long),
#mark which edit was made
#         replace_dat= ifelse(id%in%NJ14$id, 1, 0))%>%
#  dplyr::select(-c("Lat2","Long2","batch_DecMin_DD"))

#For NJ 2014 non-SESP and SALS nests, apply the average difference to those coordinates
#dat1<-dat1%>%
# lat2/long2 hold the original coordinates
#  mutate(Lat2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Lat, NA),
#         Long2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Long, NA),
# then replace lat and long with the adjusted coords
#         Lat = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Lat+diff_lat, Lat),
#         Long = ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), Long+diff_long, Long),
#         coord_shift= ifelse(site.code %in% c("AT","OC","MW") & Year==2014 & !(Species%in%c("SESP","SALS")), 1, 0))




# B-3) repeat Decimal Minute to Decimal Degree conversion for records with lat/long reversed: sites AT, OC, and MW in 2015
#dat4<-dat1%>%
# For records with coordinate information...
#    filter(missing.coords!=1 & missing.location.rec!=1 & 
# and sites AT,OC,MW in 2015...
#             (site.code %in% c("AT","OC","MW") & Year==2015))%>%
#  mutate(Long=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
#         Lat=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
#then remove the unconverted values from the UTM coordinate columns
#         Easting=NA,
#         Northing=NA,
#mark which edits were made
#        batch_DecMin_DD_reversed=1)
#dat1<-rbind(dat4,dat1[!(dat1$id%in%dat4$id),])




# B-4) Add a missing decimal to Decimal Degree coordinates (Checked plotting as UTMs and converting dMin to DD, neither worked) 
## HM, ER, BI, JC, SP, 2014 and HM, BI, ER 2015 are DD missing decimals.
#dat5<-dat1%>%
# For records with coordinate information...
#  filter(missing.coords!=1 & missing.location.rec!=1 & 
# That are also either CT and RI sites HM, ER, BI, JC, or SP in 2014 OR
# just CT sites HM, BI, ER in 2015...
#           ((site.code %in% c("HM","ER","BI","JC","SP") & Year==2014) | (site.code %in% c("HM","BI","ER") & Year==2015)))%>%
#add a decimal to the easting and northing values to turn into decimal degrees
#  mutate(Long=as.numeric(paste0(substr(Easting,1,2),".",
#                    substr(Easting,3,length(Easting)))),
#         Lat=as.numeric(paste0(substr(Northing,1,2),".",
#                                substr(Northing,3,length(Northing)))),
#then remove the original values from the UTM coordinate columns
#         Easting=NA,
#         Northing=NA,
#mark which edits were made
#         batch_dec_addition=1)

#dat1<-rbind(dat5,dat1[!(dat1$id%in%dat5$id),])



# B-5) Move Decimal Degree data in UTM easting/northing columns into long/lat columns
#dat6<-dat1%>%
# if Lat is missing and Easting is using values in Latitude range
#  filter(abs(Easting)<46 & abs(Easting)>36 & is.na(Lat))%>%
# fill in Lat from the Easting column
#  mutate(Lat=abs(Easting))%>%
#  dplyr::select(id,Lat)

#dat7<-dat1%>%
# if Lat is missing and and Northing is using values in Latitude range
#  filter(abs(Northing)<46 & abs(Northing)>36 & is.na(Lat))%>%
# fill in Lat from the NOrthing column
#  mutate(Lat=abs(Northing))%>%
#  dplyr::select(id,Lat)

#dat8<-dat1%>%
# if Long is missing and Easting is using values in Longitude range
#  filter(abs(Easting)<79&abs(Easting)>64&is.na(Long))%>%
#    # fill in Long from the Easting column
#  mutate(Long=Easting)%>%
#  dplyr::select(-Lat)

#dat9<-dat1%>%
# if Long is missing and Northing is using values in Longitude range
#  filter(abs(Northing)<79&abs(Northing)>64&is.na(Long))%>%
# fill in Long from the Northing column
#  mutate(Long=Northing)%>%
#  dplyr::select(-Lat)

#dat10<-rbind(dat8,dat9)%>%
#  left_join(rbind(dat6,dat7),by="id")%>%
#remove lat long from easting northing columns
#  mutate(Easting=NA,
#         Northing=NA,
#mark which edits were made
#         batch_move_DD_to_LatLong = 1)

#dat1<-rbind(dat10,dat1[!(dat1$id%in%dat10$id),])




# B-6) label the coordinate system and unit for nests as Decimal Degrees or as UTM, prioritizing DD
#dat1<-dat1%>%
#           mutate(Coordinate.System=case_when(
#           (if_all(c(Lat,Long),~!is.na(.)) & missing.coords!=1) ~"Lat/Long(DD)",
#           (if_all(c(Easting,Northing),~!is.na(.)) & if_all(c(Lat,Long),~is.na(.)) & missing.coords!=1) ~"UTM(m)"))





## C) Fill in remaining missing UTM zones for CT. Longitude >= 72 DD is zone 18 <72 is zone 19
######
#calc_zones <- dat1 %>% 
#  filter(!is.na(Long)|Long=="NOT REC")%>%
#  mutate(utm.zone=case_when(
#    (is.na(utm.zone) & abs(as.numeric(Long))>=72) ~ "18T",
#    (is.na(utm.zone) & abs(as.numeric(Long))<72) ~ "19T",
#    !(is.na(utm.zone)) ~ utm.zone
#  ))%>%
#  distinct(site.code,.keep_all = T)%>%
#  dplyr::select(site.code,utm.zone)

#Join site information to nest data and merge UTM zone info
#dat1<-left_join(dat1,sites,by=c("site.code"))%>%
#  left_join(calc_zones,by="site.code")%>%
#  mutate(utm.zone=case_when(
#    is.na(utm.zone) & !(is.na(utm.zone.y)) ~ utm.zone.y,
#    is.na(utm.zone) & !(is.na(utm.zone.x)) ~ utm.zone.x,
#    !(is.na(utm.zone)) ~ utm.zone))%>%
#Waterford Beach CT site, WB, should be UTM zone 18T
#  mutate(utm.zone=ifelse(site.code=="WB","18T",utm.zone))%>%
#  dplyr::select(-c("utm.zone.x","utm.zone.y","coord.system"))


#Make sure all longitude values are negative

dat1<-dat1%>%
  mutate(Long=ifelse(Long!="NOT REC"&!(is.na(Long)), -abs(as.numeric(Long)), Long))


#nest coordinate csv
output_csv<-dplyr::select(plots,id,out.bounds,wrong.site,iso.rec)%>%
  st_drop_geometry()%>%
  right_join(dplyr::select(dat1,-c("out.bounds","wrong.site","iso.rec")),by="id")%>%
  distinct(id,.keep_all = T)%>%
  mutate_at(vars(starts_with(c("iso","out","wrong","missing","batch"))), ~replace_na(., 0))%>%
  ## Add edit note information
  mutate(missing.location.rec2=ifelse(missing.location.rec==1,"Nest fate record is missing nest location record. ",""),
         missing.site.info2=ifelse(missing.site.info==1,"Missing site information. ",""),
         missing.coords2=ifelse(missing.coords==1,"Missing coordinates in nest location record. ",""),
         batch_DecMin_DD_reversed2=ifelse(batch_DecMin_DD_reversed==1,"Batch conversion for NJ sites (OC, AT, MW) in 2015: Degree Decimal Minutes to Decimal Degrees and reversed Lat and Long column values. ",""),
         replace_dat2=ifelse(replace_dat==1,"Batch conversion for NJ sites (OC, AT, MW) in 2014: Replaced coordinates with Sam R's original data. ",""),
         batch_dec_addition2=ifelse(batch_dec_addition==1,"Batch conversion for CT and RI sites (HM, BI, ER, JC, SP) in 2014 and just CT (HM, BI, ER) in 2015: Added missing decimal to Decimal Degrees in Easting/Northin cols. ",""),
         coord_shift2=ifelse(coord_shift==1,paste0("Database coordinates are shifted N. East. Added average difference from Sam R's original SALS and SESP nest records to original database coordinates (", round(Lat2,5), " N, ",round(Long2,5)," W). "),""),
         batch_move_DD_to_LatLong2 = ifelse(batch_move_DD_to_LatLong == 1, "Moved decimal degrees in the Easting Northing columns to Lat Long Columns. ", ""),
         out.bounds2=ifelse(out.bounds==1,"Nest location still plotting outside Atlantic tidal marsh. ",""),
         wrong.site2=ifelse(wrong.site==1,"Nest location plotting in a site different from the one in its records. ",""),
         iso.rec2=ifelse(iso.rec==1,"Nest location in Atlantic tidal marsh but plotting more than 500m away from any other nest record. ",""),
         Notes=paste0(missing.location.rec2,missing.site.info2,missing.coords2,batch_DecMin_DD_reversed2,replace_dat2,batch_dec_addition2, coord_shift2, out.bounds2, wrong.site2, iso.rec2, batch_move_DD_to_LatLong2),
         Notes=ifelse(Notes=="",Notes,paste0(Notes,"01/19/25 EF."))
  )%>%
  dplyr::select(-c("missing.location.rec2","missing.site.info2","missing.coords2","batch_DecMin_DD_reversed2","replace_dat2","batch_dec_addition2","out.bounds2", "wrong.site2", "iso.rec2", "coord_shift2","batch_move_DD_to_LatLong2","Lat2","Long2"))

#----------------------------------
#batch veg edits


## 4. Apply batch coordinate edits (Same batch edits for nests appear to mostly apply to the random veg sites)
#----------------------------------------------------------------------------------------------------------------

# 4-1) Degrees Decimal Minutes records (ONLY the case for NJ sites AT, OC, and MW in 2014 and 2015)

# Do 2014 first, 2015 coordinate info is in reverse columns
veg3<-veg_edits%>%
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


veg_edits<-rbind(veg3,veg_edits[!(veg_edits$veg.id%in%veg3$veg.id),])%>%
  #mark which edits were made
  mutate(batch_DecMin_DD=ifelse(veg.id%in%veg3$veg.id,1,0))

# repeat DD conversion for records with lat/long reversed: At, OC, and MW in 2015
veg4<-veg_edits%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & 
           # and sites AT,OC,MW in 2015...
           (site.code %in% c("AT","OC","MW") & Year==2015))%>%
  mutate(Long=as.numeric(substr(Northing,1,2))+(as.numeric(paste0(substr(Northing,3,4),".",substr(Northing,5,length(Northing))))/60),
         Lat=as.numeric(substr(Easting,1,2))+(as.numeric(paste0(substr(Easting,3,4),".",substr(Easting,5,length(Easting))))/60),
         #then remove the unconverted values from the UTM coordinate columns
         Easting=NA,
         Northing=NA)

veg_edits<-rbind(veg4,veg_edits[!(veg_edits$veg.id%in%veg4$veg.id),])%>%
  #mark which edits were made
  mutate(batch_DecMin_DD=ifelse(veg.id%in%veg4$veg.id,1,batch_DecMin_DD))

# add notes
veg_final<- mutate(veg_final,missing.coords2=ifelse(missing.coords==1,"Missing coordinate information. ",""),
                   batch_DecMin_DD2=ifelse(batch_DecMin_DD==1,"Batch conversion for NJ sites (OC, AT, MW) in 2014-15: Degree Decimal Minutes to Decimal Degrees and reversed Lat and Long column values in 2015. ",""),
                   batch_dec_addition2=ifelse(batch_dec_addition==1,"Batch conversion for CT and RI sites (HM, BI, ER, JC, SP) in 2014 and just CT (HM, BI, ER) in 2015: Added missing decimal to Decimal Degrees in Easting/Northin cols. ",""),
                   coord_shift2=ifelse(coord_shift==1,paste0("Database coordinates are shifted N. East. Added average difference from Sam R's original SALS and SESP nest records to original database coordinates (", round(Lat2,5), " N, ",round(Long2,5)," W). "),""),
                   batch_move_DD_to_LatLong2 = ifelse(batch_move_DD_to_LatLong == 1, "Moved decimal degrees in the Easting Northing columns to Lat Long Columns. ", ""),
                   out.bounds2=ifelse(out.bounds==1,"Location still plotting outside Atlantic tidal marsh. ",""),
                   wrong.site2=ifelse(wrong.site==1,"Location plotting in a site different from the one in its records. ",""),
                   iso.rec2=ifelse(iso.rec==1,"Location in Atlantic tidal marsh but plotting more than 500m away from any other nest record. ",""),
                   Notes_Emily=paste0(missing.coords2,batch_DecMin_DD2,batch_dec_addition2, coord_shift2, batch_move_DD_to_LatLong2,out.bounds2,wrong.site2,iso.rec2),
                   Notes_Emily=ifelse(Notes_Emily=="",Notes_Emily,paste0(Notes_Emily,"01/19/25 EF."))
)%>%
  dplyr::select(-c("missing.coords2","batch_DecMin_DD2","batch_dec_addition2", "coord_shift2","batch_move_DD_to_LatLong2","Lat2","Long2","out.bounds2","wrong.site2","iso.rec2"))


# 4-2) See if random veg plots are also plotting at an offset
# verify in arcpro
#output_shp<-dplyr::select(veg_edits,veg.id,Lat,Long)%>%
#  distinct(veg.id,.keep_all = T)%>%
#  st_as_sf(coords=c("Long","Lat"),crs=nad)%>%
#  st_transform("EPSG:26918")

#st_write(output_shp,
#         paste0(dat_path,"Nest Locations/veg_NJ_14_edits.shp"), delete_layer =T)

#Yes they are.

#offsets from nest data cleaning
diff_lat<--0.001570855 
diff_long<-0.001876374


#For NJ 2014 non-SESP and SALS nests, apply the average difference to those coordinates
veg_edits<-veg_edits%>%
  mutate(Lat2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014, Lat, NA),#this maintains the original coordinates to add to the notes at the end
         Long2 = ifelse(site.code %in% c("AT","OC","MW") & Year==2014, Long, NA),
         Lat = ifelse(site.code %in% c("AT","OC","MW") & Year==2014, Lat+diff_lat, Lat),
         Long = ifelse(site.code %in% c("AT","OC","MW") & Year==2014, Long+diff_long, Long),
         coord_shift= ifelse(site.code %in% c("AT","OC","MW") & Year==2014, 1, 0))



# 4-3) Replace all the NJ 2014-2015 nest-veg data with corrected nest coordinates to do 4-1 and 4-2 for the nest-veg records. 
NJ1415<-nests_clean%>%
  filter(Year%in%c(2014,2015) & site.code %in% c("AT","OC","MW"))%>%
  rename(Lat_nests=Lat,Long_nests=Long)
veg_edits<- veg_edits%>%
  left_join(NJ1415[,c("id","Lat_nests","Long_nests")],by=c("veg.id"="id"))%>%
  mutate(Lat = ifelse(veg.id%in%NJ1415$id, Lat_nests, Lat),
         Long = ifelse(veg.id%in%NJ1415$id, Long_nests, Long),
         #mark which edit was made
         replace_w_nest_dat= ifelse(veg.id%in%NJ1415$id, 1, 0),
         #update missing coordinate info
         missing.coords= ifelse(veg.id%in%NJ1415$id&!is.na(Lat), 0, missing.coords))%>%
  dplyr::select(-c("Lat_nests","Long_nests"))


# 4-4) Batch conversion of UTM coords that are actually DD missing decimals
# MR 2019, SJ 2019, WB 2019 the random veg plots appear to be missing decimals (some also have lat long reversed)
veg5<-veg_edits%>%
  # For records with coordinate information...
  filter(missing.coords!=1 & type=="Random"&
           ((site.code %in% c("MR","WO","SJ") & Year==2019)))%>%
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

veg_edits<-rbind(veg5,mutate(veg_edits[!(veg_edits$veg.id%in%veg5$veg.id),],batch_dec_addition=0))



# 4-5) Move DD data in easting/westing columns into long/lat columns (only applies to nest observations n =49)
veg6<-veg_edits%>%
  # if Lat is missing and Easting is using values in Latitude range
  filter(abs(Easting)<46 & abs(Easting)>36 & is.na(Lat))%>%
  # fill in Lat from the Easting column
  mutate(Lat=abs(Easting))%>%
  dplyr::select(veg.id,Lat)

veg7<-veg_edits%>%
  # if Lat is missing and and Northing is using values in Latitude range
  filter(abs(Northing)<46 & abs(Northing)>36 & is.na(Lat))%>%
  # fill in Lat from the NOrthing column
  mutate(Lat=abs(Northing))%>%
  dplyr::select(veg.id,Lat)

veg8<-veg_edits%>%
  # if Long is missing and Easting is using values in Longitude range
  filter(abs(Easting)<79&abs(Easting)>64&is.na(Long))%>%
  # fill in Long from the Easting column
  mutate(Long=Easting)%>%
  dplyr::select(-Lat)

veg9<-veg_edits%>%
  # if Long is missing and Northing is using values in Longitude range
  filter(abs(Northing)<79&abs(Northing)>64&is.na(Long))%>%
  # fill in Long from the Northing column
  mutate(Long=Northing)%>%
  dplyr::select(-Lat)

veg10<-rbind(veg8,veg9)%>%
  left_join(rbind(veg6,veg7),by="veg.id")%>%
  #remove lat long from easting northing columns
  mutate(Easting=NA,
         Northing=NA,
         #mark which edits were made
         batch_move_DD_to_LatLong = 1)

veg_edits<-rbind(veg10,mutate(veg_edits[!(veg_edits$veg.id%in%veg10$veg.id),],batch_move_DD_to_LatLong=0))



