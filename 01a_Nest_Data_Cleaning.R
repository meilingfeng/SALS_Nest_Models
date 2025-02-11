#Data tidying
library(tidyverse)
library(lubridate)
#spatial analysis
library(sf)
library(terra)



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
  dplyr::select("id"="SHARPNestID","fate"="UltimateNestFate")
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

  # site information in the 2024 demo data updates
dat_sites<-  nests%>%
    select(site.code,Site,utm.zone)%>%
  #just need to know the lat zone #, simplify the long zone letters to "T" only
  mutate(utm.zone=gsub("S","T",utm.zone))%>%
  distinct(.keep_all = T)%>%
  group_by(site.code)%>%
  slice(which.max(!is.na(utm.zone)))

  
## consolidate all the sites into one table
sites<-full_join(sites_sop,sites_sam,by=c("site.code","State"))%>%
  full_join(sites_kate,by=c("site.code","State"))%>%
  # remove any duplicate sites across the different sources
  distinct(site.code,.keep_all = T)

# add info from sites with utm or site name info that don't have info in the nest data
# first select sites that have this info
sites_utm<-sites%>%
  filter(!is.na(utm.zone))%>%
  select(site.code,utm.zone)
sites_sitename<-sites%>%
  filter(!is.na(Site))%>%
  select(site.code,Site)
#then join to the nest data records that don't have this info
dat_sites_utm<-dat_sites%>%
  filter(is.na(utm.zone))%>%
  select(-utm.zone)%>%
  left_join(sites_utm,by="site.code")
#nest data has no missing site names
dat_sites<-dat_sites%>%
  filter(!(site.code%in%dat_sites_utm$site.code))%>%
  rbind(dat_sites_utm)

#join state and regions to the nest data
dat_sites<-dat_sites%>%
  left_join(select(sites,site.code,State,Region),by="site.code")%>%
  #remove wallops island in CT, its in VA
  filter(!(site.code=="WI"&State=="CT"))%>%
  #remove WB from DE, should be woodland beach (WO), add the state info to WO
  filter(!(site.code=="WB"&State=="DE"))
dat_sites[dat_sites$site.code=="WO","State"]<-"DE"
dat_sites[dat_sites$site.code=="WO","utm.zone"]<-"18T"

#add any sites not currently in the nest data
sites_final<-dat_sites%>%
  rbind(filter(sites,!(site.code%in%dat_sites$site.code)))



## Assign UTM zones to each site
    # Greater than 72 degrees Longitude is zone 18, less than that is zone 19
    # CT is the only state split between zones 18 and 19, so we will assign CT sites once we fix the coordinates in following sections
    # for now assign zones to sites in states that are completely within 1 zone
sites_final <- sites_final%>%
  mutate(utm.zone=case_when(
    # RI, MA, NH, ME are all within UTM zone 19
    is.na(utm.zone) & State %in% c("RI","MA","NH","ME") ~ "19T",
    # NY, NJ,VA, MD, DE are in zone 18
    is.na(utm.zone) & State %in% c("NY","NJ","VA","MD","DE") ~ "18T",
    !(is.na(utm.zone)) ~ utm.zone
  ))%>%
  distinct(site.code,.keep_all = T)


## 3. Identify missing records between nest locations, fates, and sites data
# -----------------------------------------------------------------------

#first, join all the data tables together
dat1<-full_join(fates,nests,by="id")%>%
  select(-utm.zone,-Site)%>%
  left_join(sites_final,by="site.code")%>%
  # add flags for each type of coordinate error we are looking for
  mutate(missing.location.rec=0, # nest fate record is missing a nest location record
         missing.site.info=0, # record is missing site information (state, zone)
         missing.coords=0, # nest location record is missing coordinates
         out.bounds=0, # nest coordinates are not plotting in Mid Atlantic tidal marsh zone
         wrong.site=0, # nest is plotting in a site different than the one in its record
         iso.rec=0 # nest is in tidal marsh but plotting further than 500m away from any other nest record
  ) 

# Not all nests in the location records have nest fate records
missing_fates<-nests%>%
  filter(!(id%in%fates$id))
nrow(missing_fates)# number of nests missing fates, 492
  # write this list to file for documentation
if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata_2002_2024.csv"))){
write.csv(missing_fates,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/nests_missing_fatedata_2002_2024.csv"), row.names = F)
}


# Not all nests in the fate records have location data 
  #identify which ids are in the fates data but not in the location data
missing_locations<-fates%>%
  filter(!(id%in%nests$id))
nrow(missing_locations) #(44 nest fates are missing location data) (now only 8 with the 2024 update!)
  #mark which records these are
dat1[dat1$id%in%missing_locations$id,]$missing.location.rec<-1



# Site codes missing site information (site name, state) - see if SHARP folks know where these are
missing_states<- sites_final%>%
  filter(is.na(State)&is.na(utm.zone))%>%
  distinct(site.code,.keep_all = T)
  # mark with our editing flag defined above
dat1[dat1$site.code%in%missing_states$site.code,]$missing.site.info<-1




## 4. Fix any additional naming discrepancies in sites and remove rapid demo
# -------------------------------------------------------------------------------
# different codes for the same site
# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek 
# WO WB and WOBE are all Woodland Beach (fixed in 2024 data)
# SWMA=SG (saxis and guard WMA = saxis WMA)
# WI=WN (Wallops Island=Wallops Island NWR)
# HAM=HO (Hammock=Hammock River)
#JH/CHAFE=JC/John Chafee
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
    substr(site.code, start=1,stop=3)=="HAM","HO",site.code),
    id = ifelse(
      substr(id, start=1,stop=3)=="HAM",gsub("HAM","HO",id),id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=5)=="CHAFE","JC",site.code),
    id = ifelse(
      substr(id, start=1,stop=2)=="JH",gsub("JH","JC",id),id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WN","WI",site.code),
    id = ifelse(
    substr(id, start=1,stop=2)=="WN",gsub("WN","WI",id),id)
    )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=4)=="SWMA","SG",site.code),
    id = ifelse(
      substr(id, start=1,stop=2)=="SW",gsub("SW","SG",id),id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=3)=="LtR","LR2",site.code),#LR2 is in CT, LR is in ME
    id = ifelse(
      substr(id, start=1,stop=3)=="LtR",gsub("LtR","LR",id),id)
  )%>%
  #remove any sites listed as "RD" or with numbers in their site code- these are rapid demo sites and we want intensive demo
  # {2,} means at least 2 continuous numbers, some intensive demo sites have 1 digit after their code
  filter(!(site.code=="RD" | 
             grepl("^[[:digit:]]{2,}",site.code) | grepl("^[[:digit:]]{2,}",Site) |
             grepl("[[:digit:]]{2,}$",site.code) | grepl("[[:digit:]]{2,}$",Site) |
             substring(id,1,2)=="RD"))
#FR, BH, and BF are all clustered at one site in CT... I'm guessing a naming convention but not addressing here since its all pre SHARP data




## 5. Identify nest locations missing valid coordinates
# -------------------------------------------
# filter for records NOT missing coordinates- with a full set of UTM or lat long coordinate pairs
dat2<-filter(dat1,((if_all(c(Easting,Northing), ~ !is.na(.)))|(if_all(c(Lat,Long), ~ !is.na(.)))))%>%
# also coordinates that are NOT 0 or small values (if a coordinate is less than 30, it's likely a typo)
      filter((if_all(c(Easting,Northing), ~ .>30))|(if_all(c(Lat,Long), ~ abs(.)>30)))


#Mark records missing coordinates.
dat1[!(dat1$id%in%dat2$id),]$missing.coords<-1 # ids not in records with valid coordinates

sum(dat1$missing.coords) # 463, now 481 in 2024, nests with location records are missing valid coordinates


#Make sure all longitude values are negative
dat1<-dat1%>%
  mutate(Long=ifelse(Long!="NOT REC"&!(is.na(Long)), -abs(as.numeric(Long)), Long))



#remove duplicate records (from fates data)
dat1<-dat1%>%
  group_by(id)%>%
  distinct(fate,.keep_all = T)%>%
  ungroup()
  
  #for ids with conflicting fates, use the fate that is not UNKNOWN
dat1<-dat1%>%
  #order so that duplicate records with unknown fates come second. Duplicate only removes the records after the first occurrence.
  mutate(unknown=ifelse(grepl("UNKNOWN",fate),1,0))%>%
  group_by(id)%>%
  arrange(unknown)%>%
  filter(!duplicated(id))%>%
  select(-unknown)%>%
  ungroup()

## 6. Plot the nest data to look for any remaining coordinate errors
# -----------------------------------------------------------------------------
# get coordinate systems
utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#are any utm coordinates missing utm zones?
nrow(dat1%>%filter(coord.system=="UTM Zone"&is.na(utm.zone)))
#no

#convert coordinate columns to spatial points
plots_utm18 <- st_as_sf(filter(dat1,utm.zone=="18T"& grepl("UTM",coord.system)), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(dat1,utm.zone=="19T"& grepl("UTM",coord.system)), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(dat1,coord.system=="Lat/Long(DD)"), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")

#project all points to UTM 18/ NAD83
plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")%>%
  distinct(.keep_all = T)
plots100<-plots%>%
  #buffer to sample neighboring points
  st_buffer(150)
plots20<-plots%>%
  #buffer to account for spatial error in points that plot along marsh boundary
  st_buffer(20)
plots500<-plots%>%
  #buffer to remove isolated records
  st_buffer(500)




#Marsh boundaries:
# UVVR -rasterized version of national wetland inventory tidal marsh polygons
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))#uvvr overall mean dataset that I transformed coordinate systems in Arc

# Fine res marsh area from correll et al (less holes at marsh sites than UVVR layer)
# if the files don't already exist, set area outside marsh to NA
if(!file.exists(paste0(dat_path,"Correll_Marsh_Zones/Zone1_upland_stream_removed.tif"))){
  
  marsh_zones<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "DEM.tif$",full.names=T)))
  masks<-list()
  for(i in 1:length(marsh_zones)){
    masks[[i]]<-rast(marsh_zones[[i]][[1]])%>%
      terra::project("EPSG:26918")
  }
for(i in 1:length(masks)){
  mask<-masks[[i]]
  mask[mask==0|mask==9|mask==7]<-NA #(Upland and open water because these are outside the boundaries of the other layers and used as bordering cover types) 
  masks[[i]]<-mask
}

for(i in 1:length(masks)){
  writeRaster(masks[[i]],paste0(dat_path,"Correll_Marsh_Zones/Zone",i,"_upland_stream_removed.tif"))
}
#otherwise read in the marsh area files
}else{
  marsh_zones<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "_upland_stream_removed.tif$",full.names=T)))
  masks<-list()
  for(i in 1:length(marsh_zones)){
    masks[[i]]<-rast(marsh_zones[[i]][[1]])%>%
      terra::project("EPSG:26918")
}
}


# Check whether the nest is within 20m of tidal marsh using the correll layers and UVVR
# Mark as coordinate error if out of marsh bounds.
bound_check<-plots
for(i in 1:length(masks)){
  temp<-terra::extract(masks[[i]],vect(plots20),fun=min,na.rm=T,bind=T)%>%
    sf::st_as_sf()%>%
    st_drop_geometry()
  bound_check<-bound_check%>%
    left_join(temp)
}

bound_check2<-terra::extract(marsh,vect(plots20),fun=min,na.rm=T,bind=T)%>%
  sf::st_as_sf()%>%
  st_drop_geometry()%>%
  right_join(bound_check)%>%
  mutate(out.bounds2=ifelse(
    uvvr_mean_utm18_2=="NaN"&Z1_DEM=="NaN"&Zone2_DEM=="NaN"&Zone3_DEM=="NaN"&Zone4_DEM=="NaN"&Zone5_DEM=="NaN"&
      Zone6_DEM=="NaN"&Zone7_DEM=="NaN"&Zone8_DEM=="NaN",1,0
  ))%>%
  dplyr::select(id,out.bounds2)

# Mark as error if nest has a different site code than its neighbors 
# some records are in marsh sites but have the wrong site code. Could be recording error or plotting error for some. Chose to just remove all.
nest_neighb <-st_join(plots100,plots)%>%
  group_by(id.x,site.code.x,site.code.y)%>%
  count()%>%
  ungroup()%>%
  group_by(id.x,site.code.x)%>%
  #arrange sites within each nest's buffer by the frequency they occur
  arrange(id.x,n)%>%
  # get the most frequent site in nests neighbors and the number of different sites in its neighbors' attributes
  summarize(site.code.y=last(site.code.y),
            n_neighbors=last(n),
            n_sites=n())%>%
  ungroup()
  # mark record if it does not match the most frequent site in its neighbors and has more than 1 site in its neighbors
nest_neighb2<-nest_neighb%>%
  mutate(wrong.site2=ifelse((site.code.x!=site.code.y)&(n_sites>1),1,0))%>%
  select(id=id.x,wrong.site2)


#Mark as error if nest is an isolated record located more than 500m away from any other records
isolated_rec<-st_join(plots500,plots)%>%
  group_by(id.x)%>%
  count()%>%
  ungroup()%>%
  mutate(iso.rec2=ifelse(n==1,1,0))%>%
  select(id=id.x,iso.rec2)

#join all the coordinate error flags together
plots<-plots%>%
  left_join(st_drop_geometry(bound_check2),by="id")%>%
  mutate(out.bounds=out.bounds2)%>%
  left_join(st_drop_geometry(nest_neighb2),by="id")%>%
  mutate(wrong.site=wrong.site2)%>%
  left_join(st_drop_geometry(isolated_rec),by="id")%>%
  mutate_at(vars(starts_with(c("iso","out","wrong"))), ~replace_na(., 0))%>%
  #if the record is out of marsh boundaries, prioritize that flag over the isolation flag
  mutate(iso.rec=ifelse(iso.rec2==1&out.bounds==0,1,0))%>%
  dplyr::select(-out.bounds2, -wrong.site2, -iso.rec2)

 # check in ArcPro
if(!(file.exists(paste0(path_out,"Intermediate_outputs/Nest_locations/nest_edit_locations_01_29_25.shp")))){
  st_write(plots,
           paste0(path_out,"Intermediate_outputs/Nest_locations/nest_edit_locations_01_29_25.shp"), delete_layer =T)
}
## Notes based on plotting the data
# Need to class state==CT and year<2011 as not isolated. Exclude sites in CT pre-2010 from isolated record flag.
  # Pre SHARP CT sites tend to only have single records:
    #GM06SALS533,BP09SALS020,LI South08SALS011,RR09SESP006,M/H06SALS548,OR06RWBL514
  #PL19018,BB19011 also seem like valid single record sites in maine
  #LP13WILL002 seems valid single record in NY

#Need to class as out.bounds
  #HM15WILL100, HM15SALS001 are incorrectly plotting next to each other in a different marsh area. 
  #Same with JO15VIRA127 and JO13SALS110.
  #WR08SALS010 doesnt look right, ploting out of marsh boundary.
  #MN19021,NC18SALS102,FS14SALS183 are next to misclassified marsh pixels, should be out of bounds
  #FB21028, EL15SALS016

# Need to class as wrong.site
  #JC12SALS102,JC12SALS101 ,JC23SALS001
  #ER19088  
  #LU24023

#Site MS is right next to site DN, MR21033,MR21044,MR21056,MR20035,MR21036,MR21009,MR21040 are not in the wrong site
#north chapmans is also on the border of chapmans, 	CN23037,CN23012,CN23005,CN23050 not in the wrong site.

plots<-plots%>%
  mutate(iso.rec=ifelse(iso.rec==1&((State=="CT"&Year<2011)|id%in%c("PL19018","BB19011","LP13WILL002","SJ21005")),0,iso.rec),
         out.bounds=ifelse(id%in%c("HM15WILL100", "HM15SALS001","JO15VIRA127","JO13SALS110","WR08SALS010","MN19021",
                                   "NC18SALS102","FS14SALS183","FB21028", "EL15SALS016"),1,out.bounds),
         wrong.site=ifelse(id%in%c("JC12SALS102","JC12SALS101","ER19088","LU24023","JC23SALS001"),1,wrong.site),
         wrong.site=ifelse(id%in%c("MR21033","MR21044","MR21056","MR20035","MR21036","MR21009","MR21040",
                                   "CN23037","CN23012","CN23005","CN23050"),0,wrong.site))


## 5. Write outputs to file
#--------------------------------------------

#nest points shapefile and kml file
output_shp<-dplyr::select(plots,id,out.bounds,wrong.site,iso.rec)%>%
  left_join(dplyr::select(dat1,-c("out.bounds","wrong.site","iso.rec")),by="id")%>%
  distinct(id,.keep_all = T)

st_write(output_shp,
         paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_01_29_25.shp"), delete_layer =T)#Jan 29, 2025 layer uses the 2024 data update
st_write(output_shp,
         paste0(path_out,"Intermediate_outputs/Nest_locations/nest_locations_01_29_25.kml"), delete_layer =T)


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
         out.bounds2=ifelse(out.bounds==1,"Nest location still plotting outside Atlantic tidal marsh. ",""),
         wrong.site2=ifelse(wrong.site==1,"Nest location plotting in a site different from the one in its records. ",""),
         iso.rec2=ifelse(iso.rec==1,"Nest location in Atlantic tidal marsh but plotting more than 500m away from any other nest record. ",""),
         Notes=paste0(missing.location.rec2,missing.site.info2,missing.coords2, out.bounds2, wrong.site2, iso.rec2),
         Notes=ifelse(Notes=="",Notes,paste0(Notes,"01/19/25 EF."))
  )%>%
  dplyr::select(-c("missing.location.rec2","missing.site.info2","missing.coords2"))

if(!file.exists(paste0(path_out,"Intermediate_outputs/new_nest_coords_01_29_25.csv"))){
write.csv(output_csv,paste0(path_out,"Intermediate_outputs/new_nest_coords_01_29_25.csv"),row.names = F)
}


#Site table info
output_site<-output_csv%>%
  dplyr::select(site.code,Site,State,utm.zone)%>%
  distinct(site.code,.keep_all = T)

write.csv(output_site,paste0(path_out,"Final_outputs/compiled_site_table_01_29_25.csv"),row.names = F)
