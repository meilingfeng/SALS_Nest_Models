library(tidyverse)
library(lubridate)
#spatial analysis
library(sf)
library(terra)




########################################################
# Fix any coordinate issues with the vegetation plots
########################################################



## Set file path to data
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"



## 1. Load data
# -------------------------------------------

# Cleaned nest data including fates without coordinates
nests_clean<- read.csv(paste0(path_out,"Intermediate_Outputs/new_nest_coords_01_29_25.csv"))

# Original Nest location data
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2002-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site","Site"="SiteName", "Year", "Species",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))%>%
  #remove rapid demo records
  filter(!(site.code=="RD" | 
             grepl("^[[:digit:]]{2,}",site.code) | grepl("^[[:digit:]]{2,}",Site) |
             grepl("[[:digit:]]{2,}$",site.code) | grepl("[[:digit:]]{2,}$",Site) |
             substring(id,1,2)=="RD"))


# vegetation plots
veg<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("veg.id"="VegPointID","type"="PointType","site.code"="Site","Site"="SiteName","date"="SurveyDate",
                "coord.system"="Coordinate.System","utm.zone"="UTM.Zone","Easting","Northing","Lat","Long")%>%
  mutate(type=case_when(
    grepl("R",type)~"Random",
    grepl("N",type)~"Nest"),
         Year = as.numeric(substr(date,nchar(date)-3,nchar(date))))%>%
  distinct(.keep_all = T)%>%
  #remove rapid demo records
  filter(!(site.code=="RD" | 
             grepl("^[[:digit:]]{2,}",site.code) | grepl("^[[:digit:]]{2,}",Site) |
             grepl("[[:digit:]]{2,}$",site.code) | grepl("[[:digit:]]{2,}$",Site) |
             substring(veg.id,1,2)=="RD"))


## Compare similar records between vegetation data and nest data
similar<-veg[veg$veg.id%in%nests$id,] # what nest records in the nest location data have vegetation survey data?
nrow(similar)
#5061 shared plots between veg and nest data (7086 in 2024)

nrow(veg[veg$type=="Nest",])
#5075 total nest plots in the veg data (6829 in 2024) 
# (about 14 extra nest records in the veg data) (47 in 2024)
extra<-veg[!(veg$veg.id%in%nests$id)&veg$type%in%c("Nest","NEST"),]

nrow(veg[veg$type=="Random",])
#5930 random vegetation plots not at nests (7209 in 2024)




## 2. Site information - select Site code, site name, and state
# -------------------------------------------------------------
# Fix site name discrepancies in the veg data
# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek - doesn't seem to be an issue in this dataset
sort(unique(veg$site.code))
# WA site is actually BI change the one record
# SC15RAND049 has a conflict between site name and the site code used in the point ID- remove
# NC19RAND and MN19RAND and ID19RAND points all have a site name conflict- change from NC1 to NC and MN1 to MN and ID1/2/A/B to just ID
#	FS13SALS and SA13SALS listed as ID for site.code- change to FS and SA
veg<-veg%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WA","BI",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="WA",gsub("WA","BI",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WB","WO",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="WB",gsub("WB","WO",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=3)=="HAM","HO",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=3)=="HAM",gsub("HAM","HO",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=5)=="CHAFE","JC",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="JH",gsub("JH","JC",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=2)=="WN","WI",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="WN",gsub("WN","WI",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=4)=="SWMA","SG",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=2)=="SW",gsub("SW","SG",veg.id),veg.id)
  )%>%
  mutate(site.code = ifelse(
    substr(site.code, start=1,stop=3)=="LtR","LR2",site.code),
    veg.id = ifelse(
      substr(veg.id, start=1,stop=3)=="LtR",gsub("LtR","LR",veg.id),veg.id)
  )%>%
  filter(veg.id!="SC15RAND049")%>%
  mutate(site.code=case_when(
    site.code%in%c("NC1","NC2")~"NC",
    site.code%in%c("MN1","MN2")~"MN",
    site.code%in%c("ID3", "ID1", "ID2", "IDA", "IDB")~"ID",
    site.code%in%c("FS1")~"FS",
    substr(veg.id,1,8)=="FS13SALS"~"FS",
    substr(veg.id,1,4)=="SA13"~"SA",
    substr(veg.id,1,8)=="MN13SALS"~"MN",
    !(substr(veg.id,1,8)%in%c("SA13SALS","FS13SALS","MN13SALS")|site.code%in%c("NC1","NC2","MN1","MN2","ID3", "ID1", "ID2", "IDA", "IDB"))
    ~site.code))



#site info from nest data
sites <- nests_clean%>%
  dplyr::select(site.code,State,utm.zone)%>%
  distinct(site.code,.keep_all = T)
sites[sites$site.code=="SJ",]$utm.zone<-"18T"

# additional sites in veg data
veg_sites<-unique(veg[!(veg$site.code%in%nests_clean$site.code),"site.code"])
#none that aren't alternative names for ones already in nest data

# Join utm zones to veg data
veg<-left_join(dplyr::select(veg,-utm.zone),sites,by="site.code")




## 3. Mark missing coordinate info in the veg data
# -------------------------------------------------------------------
veg_edits<-veg
# Mark veg records missing valid coordinates
      # records that are NOT missing coordinates
veg2<-filter(veg_edits,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))%>%
      # Records that do NOT have coordinates that are 0 or small values (less than 30, likely typos)
  filter((if_all(c(Easting,Northing), ~ .>30))|(if_all(c(Lat,Long), ~ abs(.)>30)))

veg_edits$missing.coords<-0

      # Records not in filter are missing or have invalid coordinates
veg_edits[!(veg_edits$veg.id%in%veg2$veg.id),]$missing.coords<-1

      #3780 missing coordinates (442 in 2024 update)
nrow(veg_edits[veg_edits$missing.coords==1,])

ggplot(veg_edits,aes(x=as.factor(missing.coords),color=type))+
  geom_histogram(stat = "count",aes(fill=type))
#majority of missing coordinates are in random plot data
nrow(veg_edits[veg_edits$missing.coords==1&veg_edits$type=="Random",])
# 329 random plots are missing coords
rand_missing<-veg_edits[veg_edits$missing.coords==1&veg_edits$type=="Random",]




## 4. Misc Tidying
#------------------------------------------------------------------------------------------------
#set standard names for coordinate systems
veg_edits<-veg_edits%>%
  mutate(coord.system=case_when(
    #if lat long is filled and the record isn't missing valid coordinates, mark as DD
    (if_all(c(Lat,Long),~!is.na(.)) & missing.coords!=1) ~"Lat/Long(DD)",
    #if easting northing is filled and the record isn't missing valid coordinates, mark as UTM
    (if_all(c(Easting,Northing),~!is.na(.)) & missing.coords!=1) ~"UTM(m)",
    (if_all(c(Lat,Long),~!is.na(.)) & if_all(c(Easting,Northing),~!is.na(.))) ~"Lat/Long(DD), UTM(m)"))



#Make sure all longitude values are negative
veg_edits<-veg_edits%>%
  mutate(Long=ifelse(Long!="NOT REC", -abs(as.numeric(Long)), Long))



#remove duplicated record IDs (most are two different random veg plots with the same ID)
veg_edits2<-veg_edits%>%
  distinct(.keep_all = T)
#some IDs are shared between nests and their random veg point, only remove duplicates that are the same type
# remaining dups appear to be resurveys, remove those
veg_edits<-veg_edits[!duplicated(veg_edits[c('veg.id', 'type')]),]

#add a way to distinguish random and nest surveys for those that share the same ID
veg_edits<-veg_edits%>%
  arrange(type)%>%
  mutate(veg.id2=ifelse(duplicated(veg.id)&type=="Random",paste0(veg.id,"_rand"),veg.id))


# 5) Look for additional coordinate errors by plotting the data
##---------------------------------------------------------------------------------------
# get coordinate systems and convert coordinates to spatial data
utm18<- "+proj=utm +zone=18 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
nad<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

#convert to spatial points
plots_utm18 <- st_as_sf(filter(veg_edits,utm.zone=="18T"& coord.system=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)

plots_utm19 <- st_as_sf(filter(veg_edits,utm.zone=="19T"& coord.system=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)

plots_latlong <- st_as_sf(filter(veg_edits,coord.system%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")


plots<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")

plots150<-plots%>%
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
if(!exists("marsh")){
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))#uvvr overall mean dataset that I transformed coordinate systems in Arc

# Fine res marsh area from correll et al (less holes at marsh sites than UVVR layer)
  marsh_zones<-unlist(map(paste0(dat_path,"Correll_Marsh_Zones"),~list.files(.,pattern = "_upland_stream_removed.tif$",full.names=T)))
  masks<-list()
  for(i in 1:length(marsh_zones)){
    masks[[i]]<-rast(marsh_zones[[i]][[1]])%>%
      terra::project("EPSG:26918")
  }

}
  
# Check whether the record is within 20m of tidal marsh using the correll layers and UVVR
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
    mutate(out.bounds=ifelse(
      uvvr_mean_utm18_2=="NaN"&Z1_DEM=="NaN"&Zone2_DEM=="NaN"&Zone3_DEM=="NaN"&Zone4_DEM=="NaN"&Zone5_DEM=="NaN"&
        Zone6_DEM=="NaN"&Zone7_DEM=="NaN"&Zone8_DEM=="NaN",1,0
    ))%>%
    dplyr::select(veg.id2,out.bounds)
  
  
# Mark as error if record has a different site code than its neighbors 
  # some records are in marsh sites but have the wrong site code. Could be recording error or plotting error for some. Chose to just remove all.
  nest_neighb <-st_join(plots150,plots)%>%
    group_by(veg.id2.x,site.code.x,site.code.y)%>%
    count()%>%
    ungroup()%>%
    group_by(veg.id2.x,site.code.x)%>%
    #arrange sites within each record's buffer by the frequency they occur
    arrange(veg.id2.x,n)%>%
    # get the most frequent site in record's neighbors and the number of different sites in its neighbors' attributes
    summarize(site.code.y=last(site.code.y),
              n_neighbors=last(n),
              n_sites=n())%>%
    ungroup()
  # mark record if it does not match the most frequent site in its neighbors and has more than 1 site in its neighbors
  nest_neighb2<-nest_neighb%>%
    mutate(wrong.site=ifelse((site.code.x!=site.code.y)&(n_sites>1),1,0))%>%
    select(veg.id2=veg.id2.x,wrong.site)
  
  
#Mark as error if record is isolated- located more than 500m away from any other records
  isolated_rec<-st_join(plots500,plots)%>%
    group_by(veg.id2.x)%>%
    count()%>%
    ungroup()%>%
    mutate(iso.rec=ifelse(n==1,1,0))%>%
    select(veg.id2=veg.id2.x,iso.rec)
  
  plots2<-plots%>%
    left_join(st_drop_geometry(bound_check2),by="veg.id2")%>%
    left_join(st_drop_geometry(nest_neighb2),by="veg.id2")%>%
    left_join(st_drop_geometry(isolated_rec),by="veg.id2")%>%
    mutate_at(vars(starts_with(c("iso","out","wrong"))), ~replace_na(., 0))%>%
    mutate(iso.rec=ifelse(iso.rec==1&out.bounds==0,1,0))
  


# Verify typos in arcpro. Compare with visual detections Intermediate_Outputs/Nest_locations/veg_locations_error_edits_12_26_22.shp
output_shp<-dplyr::select(plots2,veg.id2,starts_with(c("iso","out","wrong")))%>%
  distinct(veg.id2,.keep_all = T)
  
st_write(output_shp,
           paste0(path_out,"Intermediate_Outputs/Nest_locations/veg_locations_error_edits_01_29_25.shp"), delete_layer =T)

#Notes:
# change to out.bounds HM13CLRA003, HM15RAND010, HM12RAND728, HM14RAND287, ER12SALS046, JO14SALS110, 	EL15SALS039
  # SC13RAND037,SC13RAND044,JO15RAND355, JO13RAND052, NO15RAND225, NC19209, NC19213,NC18RAND123,NC18RAND122,MQ17SALS001,EL11SALS030
  # MN19RAND155,MN19RAND155,ID15RAND126,FS14RAND102,FS14SALS183,FS14RAND123,FS14RAND103,FS15RAND102,FS15RAND101,SA13RAND003,MW13RAND206
# JO12RAND000,NO14RAND204,LU24023 don't appear to be in their sites, change to wrong.site
# FH24008 is in the right site, remove wrong.site flag
# OC11RAND002 change to isolated

plots2<-plots2%>%
  mutate(out.bounds=ifelse(veg.id2%in%c("HM13CLRA003", "HM15RAND010", "HM12RAND728", "HM14RAND287", "ER12SALS046", "JO14SALS110", 	"SC13RAND037",
                                   "SC13RAND044", "JO15RAND355", "JO13RAND052", "NO15RAND225", "NC19209", "NC19213","NC18RAND123","NC18RAND122",
                                   "MN19RAND155","MN19RAND155", 	"ID15RAND126", 	"FS14RAND102",	"FS14SALS183","FS14RAND123","FS14RAND103",
                                   "FS15RAND102", "FS15RAND101","SA13RAND003",	"MW13RAND206","	MN14RAND283","MQ17SALS001",
                                   "EL15SALS039","EL11SALS030"),1,out.bounds),
         wrong.site=ifelse(veg.id2%in%c("JO12RAND000", 	"NO14RAND204","LU24023"),1,wrong.site),
         wrong.site=ifelse(veg.id2%in%c("FH24008"),0,wrong.site),
         iso.rec=ifelse(veg.id2%in%c("OC11RAND002"),1,iso.rec))



veg_edits <-dplyr::select(plots2,veg.id2,out.bounds,wrong.site,iso.rec)%>%
  st_drop_geometry()%>%
  right_join(veg_edits,by="veg.id2")%>%
  mutate_at(vars(starts_with(c("iso","out","wrong"))), ~replace_na(., 0))


## 6. Fill in any missing/typo coordinate data from the cleaned nest dataset
#---------------------------------------------------------------------------------------------
# get all shared record ids between veg and nest data
veg_nests<-veg_edits%>%
  filter(veg.id2%in%nests_clean$id)%>%
  dplyr::select(veg.id2, site.code.veg=site.code, Easting.veg=Easting, Northing.veg=Northing, Lat.veg=Lat, Long.veg=Long, utm.zone,
                out.bounds.veg=out.bounds, wrong.site.veg=wrong.site, iso.rec.veg=iso.rec, missing.coords.veg=missing.coords)%>%
  left_join(nests_clean[,c("id","site.code","Easting","Northing","Lat","Long","out.bounds","wrong.site","iso.rec","missing.coords")],by=c("veg.id2"="id"))


#records that have different typo information between the datasets 
#use coordinate info from the data source without typos.
# don't include isolated record flags here because they might change as more data is added.
veg_nests<-veg_nests%>%
  mutate(iso.rec.veg=ifelse(iso.rec==1,1,iso.rec.veg))
veg12<-veg_nests%>%
  mutate(nest.error=ifelse(if_any(c(out.bounds,wrong.site,iso.rec,missing.coords),~.==1),1,0),
         veg.error=ifelse(if_any(c(out.bounds.veg,wrong.site.veg, iso.rec.veg, missing.coords.veg),~.==1),1,0))%>%
  #filter for records with coordinate errors in veg OR nest data. By including missing records and errors together, veg data missing coordinates won't use coordinates from nest data if there looks like something is wrong with the coordinate info.
  filter(veg.error!=nest.error)%>%
        #if there is no veg error and an error in the nest data, use veg location info. Vice versa.
        #mark the record ID with which data to use.
  mutate(coord.dat=ifelse(veg.error==0,"veg","nest"))
  #should fill in missing coordinates in veg data and fix typos
sum(veg12$nest.error)#n nest typos resolved
sum(veg12[veg12$missing.coords.veg!=1,]$veg.error)# n veg typos resolved
sum(veg12[veg12$missing.coords.veg==1,]$veg.error) # n missing veg coords filled




# Fill in veg data with missing or typo coordinates with corresponding nest data coordinates
# applying to veg edits automatically retains veg coordinates not flagged as typos
veg13<-veg_edits%>%
  # filter to nests in nest data and that were not flagged to use veg data (these are just the WB and SJ sites)
  filter(veg.id2%in%veg12[veg12$coord.dat=="nest",]$veg.id2)%>%
  select(-c("site.code","Easting","Northing","Lat","Long","Site","missing.coords","out.bounds","wrong.site","iso.rec","coord.system"))%>%
  left_join(select(nests_clean,id,site.code,Easting,Northing,Lat,Long,out.bounds,wrong.site,iso.rec,Site,missing.coords,coord.system),by=c("veg.id2"="id"))%>%
  mutate(data.source="Demo_Nests")
veg_edits<-rbind(mutate(filter(veg_edits,!(veg.id2%in%veg13$veg.id2)),data.source="Demo_Veg"),veg13)

#fix typos in the nest data by replacing with correct coordinates in veg data
nests_clean2<-nests_clean%>%
  filter(id%in%veg12[veg12$coord.dat=="veg",]$veg.id2)%>%
  select(-c("site.code","Easting","Northing","Lat","Long","Site","missing.coords","out.bounds","wrong.site","iso.rec","coord.system"))%>%
  left_join(dplyr::select(veg_edits,veg.id2,site.code,Easting,Northing,Lat,Long,out.bounds,wrong.site,iso.rec,Site,missing.coords,coord.system),by=c("id"="veg.id2"))%>%
  mutate(data.source="Demo_Veg")
nests_clean3<-rbind(mutate(filter(nests_clean,!(id%in%nests_clean2$id)),data.source="Demo_Nests"),nests_clean2)%>%
# add new note of changes (mostly a repeated number deleted)
        mutate(Notes=ifelse(data.source=="Demo_Veg","Replaced with correct coordinates in vegetation dataset. 01/20/25 EF. ",Notes))




# 7. Check for any changes in error flags after combining/adding data points across data sources
#-------------------------------------------------------------------------------------------------
#convert veg back to spatial points and check if there are any more typos now that we added new coordinate info from the nest data
plots_utm18 <- st_as_sf(filter(veg_edits,utm.zone=="18T"& coord.system=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)%>%
  rename("id"="veg.id2")

plots_utm19 <- st_as_sf(filter(veg_edits,utm.zone=="19T"& coord.system=="UTM(m)"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)%>%
  rename("id"="veg.id2")

plots_latlong <- st_as_sf(filter(veg_edits,coord.system%in%c("Lat/Long(DD)","Lat/Long(DD), UTM(m)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("id"="veg.id2","Long"="Easting","Lat"="Northing")

#Combine records from the nest data too
plots_utm18_nests <- st_as_sf(filter(nests_clean3,utm.zone=="18T"& grepl("UTM",coord.system)), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(nad)%>%
  mutate(type="Nest")

plots_utm19_nests <- st_as_sf(filter(nests_clean3,utm.zone=="19T"& grepl("UTM",coord.system)), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(nad)%>%
  mutate(type="Nest")

plots_latlong_nests <- st_as_sf(filter(nests_clean3,coord.system%in%c("Lat/Long(DD)")), coords = c("Long", "Lat"), crs = nad)%>%
  #rename Lat Long to Northing Easting to combine data below
  rename("Long"="Easting","Lat"="Northing")%>%
  mutate(type="Nest")


plots_veg<-rbind(plots_utm18,plots_utm19,plots_latlong)%>%
  st_transform("EPSG:26918")%>%
  select(id,type,site.code,wrong.site,iso.rec,out.bounds,Site,missing.coords, coord.system,State,Year)
plots_nests<-rbind(plots_utm18_nests,plots_utm19_nests,plots_latlong_nests)%>%
  st_transform("EPSG:26918")%>%
  select(id,type,site.code,wrong.site,iso.rec,out.bounds,Site,missing.coords, coord.system,State,Year)
plots<-rbind(plots_veg,plots_nests)%>%
  distinct(id,.keep_all = T)

plots150<-plots%>%
  #buffer to sample neighboring points
  st_buffer(150)
plots20<-plots%>%
  #buffer to account for spatial error in points that plot along marsh boundary
  st_buffer(20)
plots500<-plots%>%
  #buffer to remove isolated records (larger than nest records because random points are generated within the entire study plot)
  st_buffer(500)

# Check whether the record is within 20m of tidal marsh using the correll layers and UVVR
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
  mutate(out.bounds=ifelse(
    uvvr_mean_utm18_2=="NaN"&Z1_DEM=="NaN"&Zone2_DEM=="NaN"&Zone3_DEM=="NaN"&Zone4_DEM=="NaN"&Zone5_DEM=="NaN"&
      Zone6_DEM=="NaN"&Zone7_DEM=="NaN"&Zone8_DEM=="NaN",1,0
  ))%>%
  dplyr::select(id,out.bounds)


# Mark as error if record has a different site code than its neighbors 
# some records are in marsh sites but have the wrong site code. Could be recording error or plotting error for some. Chose to just remove all.
nest_neighb <-st_join(plots150,plots)%>%
  group_by(id.x,site.code.x,site.code.y)%>%
  count()%>%
  ungroup()%>%
  group_by(id.x,site.code.x)%>%
  #arrange sites within each record's buffer by the frequency they occur
  arrange(id.x,n)%>%
  # get the most frequent site in record's neighbors and the number of different sites in its neighbors' attributes
  summarize(site.code.y=last(site.code.y),
            n_neighbors=last(n),
            n_sites=n())%>%
  ungroup()
# mark record if it does not match the most frequent site in its neighbors and has more than 1 site in its neighbors
nest_neighb2<-nest_neighb%>%
  mutate(wrong.site=ifelse((site.code.x!=site.code.y)&(n_sites>1),1,0))%>%
  select(id=id.x,wrong.site)


#Mark as error if record is isolated- located more than 500m away from any other records
isolated_rec<-st_join(plots500,plots)%>%
  group_by(id.x)%>%
  count()%>%
  ungroup()%>%
  mutate(iso.rec=ifelse(n==1,1,0))%>%
  select(id=id.x,iso.rec)

plots2<-plots%>%
  left_join(st_drop_geometry(rename(bound_check2,out.bounds2=out.bounds)),by="id")%>%
  left_join(st_drop_geometry(rename(nest_neighb2,wrong.site2=wrong.site)),by="id")%>%
  left_join(st_drop_geometry(rename(isolated_rec,iso.rec2=iso.rec)),by="id")%>%
  mutate_at(vars(c("iso.rec2","out.bounds2","wrong.site2")), ~replace_na(., 0))%>%
  mutate(iso.rec2=ifelse(iso.rec2==1&out.bounds2==0,1,0))




# Notes to add to previous ArcPro visual notes:
  #LU15SALS005,EL15SALS016,	SC13RAND042, 	MQ17SALS001, 	NC19RAND118, 	NC19205,JO14STSP056,NC18SALS102,MQ17SALS001,MQ16SALS013
  # MN14RAND283,	MQ16SALS013,OC15SESP312,OC14WILL220,AT15SESP117,NC19068,NC19204 look out of marsh bounds
  #MR19RAND293 looks isolated
  #OC11RAND002 is not isolated any more, remove from previous edits
  #NC19209, NC19213, HM13CLRA003, ER12SALS046 are in bounds now, remove from previous edits
  # actually in the right site: CN23RAND005,CN23RAND012,CN23RAND037,N23RAND050
  #wrong site JC23SALS003,JC23SALS001

plots2<-plots2%>%
  mutate(out.bounds2=ifelse(id%in%c("HM15RAND010","HM12RAND728","HM14RAND287","JO14SALS110","SC13RAND037","OC14WILL220",
                                       "SC13RAND044","JO15RAND355","JO13RAND052","NO15RAND225","NC18RAND123","NC18RAND122","AT15SESP117",
                                       "MN19RAND155","MN19RAND155", "ID15RAND126","FS14RAND102","FS14SALS183","FS14RAND123","FS14RAND103",
                                       "FS15RAND102","FS15RAND101","SA13RAND003","MW13RAND206","MN14RAND283","OC15SESP312",
                                       "LU15SALS005","EL11SALS030","EL15SALS016","SC13RAND042", "MQ17SALS001", "NC19RAND118","MN14RAND283",
                                       "MQ16SALS013","NC19205","NC19204","NC19068","JO14STSP056","NC18SALS102","MQ17SALS001","MQ16SALS013"),1,out.bounds2),
         wrong.site2=ifelse(id%in%c("JO12RAND000", 	"NO14RAND204","JC23SALS003","JC23SALS001"),1,wrong.site2),
         wrong.site2=ifelse(id%in%c("CN23RAND005","CN23RAND012","CN23RAND037","N23RAND050"),0,wrong.site2),
         iso.rec2=ifelse(id%in%c("MR19RAND293"),1,iso.rec2))

# Notes to add from nest data cleaning
# Need to class state==CT and year<2011 as not isolated. Exclude sites in CT pre-2010 from isolated record flag.
# Pre SHARP CT sites tend to only have single records:
#GM06SALS533,BP09SALS020,LI South08SALS011,RR09SESP006,M/H06SALS548,OR06RWBL514
#PL19018,BB19011 also seem like valid single record sites in maine
#LP13WILL002 seems valid single record in NY

#Need to class as out.bounds
#HM15WILL100, HM15SALS001 are incorrectly plotting next to each other in a different marsh area. 
#Same with JO15VIRA127 and JO14SALS110.
#WR08SALS010 doesnt look right, ploting out of marsh boundary.
#NC18SALS102,FS14SALS183 are next to misclassified marsh pixels, should be out of bounds; veg edits fixed this for MN19021 and JO15VIRA127, remove
#FB21028, EL15SALS016
# new- MM20026 looks isolated


# Need to class as wrong.site
#JC12SALS102,JC12SALS101 ,JC23SALS001
#ER19088  
#LU24023

#Site MS is right next to site DN, MR21033,MR21044,MR21056,MR20035,MR21036,MR21009,MR21040 are not in the wrong site
#north chapmans is also on the border of chapmans, 	CN23037,CN23012,CN23005,CN23050 not in the wrong site.

plots2<-plots2%>%
mutate(iso.rec2=ifelse(iso.rec2==1&((State=="CT"&Year<2011)|id%in%c("PL19018","BB19011","LP13WILL002","SJ21005","MM20026")),0,iso.rec2),
       out.bounds2=ifelse(id%in%c("HM15WILL100", "HM15SALS001","JO14SALS110","WR08SALS010",
                                 "NC18SALS102","FS14SALS183","FB21028", "EL15SALS016"),1,out.bounds2),
       wrong.site2=ifelse(id%in%c("JC12SALS102","JC12SALS101","ER19088","LU24023","JC23SALS001"),1,wrong.site2),
       wrong.site2=ifelse(id%in%c("MR21033","MR21044","MR21056","MR20035","MR21036","MR21009","MR21040",
                                 "CN23037","CN23012","CN23005","CN23050"),0,wrong.site2))


# Once again check in ArcPro
if(!(file.exists(paste0(path_out,"Intermediate_outputs/Nest_locations/veg_nest_edit_locations_01_29_25.shp")))){
  st_write(plots2,
           paste0(path_out,"Intermediate_outputs/Nest_locations/veg_nest_edit_locations_01_29_25.shp"), delete_layer =T)
}


sum(plots_veg$out.bounds)+sum(plots_veg$iso.rec)+sum(plots_veg$wrong.site)
#176 typos in veg total remain (269 in 2024)
nrow(veg_edits)-nrow(plots_veg)
#401 (now 461) veg records still missing valid coordinate data


#combine records with coordinate info back with records missing coordinate info
#veg data
veg14 <-plots2%>%
  #replace old edit flags with new flags
  select(veg.id2=id,out.bounds=out.bounds2,wrong.site=wrong.site2,iso.rec=iso.rec2)%>%
  st_drop_geometry()
veg_final<-veg_edits%>%
  select(-c("out.bounds","wrong.site","iso.rec"))%>%
  left_join(veg14,by="veg.id2")%>%
  mutate_at(vars(starts_with(c("iso","out","wrong"))), ~replace_na(., 0))
# add notes
veg_final<- mutate(veg_final,missing.coords2=ifelse(missing.coords==1,"Missing coordinate information. ",""),
        out.bounds2=ifelse(out.bounds==1,"Location still plotting outside Atlantic tidal marsh. ",""),
        wrong.site2=ifelse(wrong.site==1,"Location plotting in a site different from the one in its records. ",""),
        iso.rec2=ifelse(iso.rec==1,"Location in Atlantic tidal marsh but plotting more than 500m away from any other nest record. ",""),
        Notes_Emily=paste0(missing.coords2,out.bounds2,wrong.site2,iso.rec2),
        Notes_Emily=ifelse(Notes_Emily=="",Notes_Emily,paste0(Notes_Emily,"01/19/25 EF."))
  )%>%
  dplyr::select(-c("missing.coords2","out.bounds2","wrong.site2","iso.rec2"))


#3780 missing coordinates; now reduced to 413
nrow(veg_final[veg_final$missing.coords==1,])
ggplot(veg_final,aes(x=as.factor(missing.coords),color=type))+
  geom_histogram(stat = "count",aes(fill=type))




# nest data
nest_updates <-plots2%>%
  #replace old edit flags with new flags
  select(id,out.bounds=out.bounds2,wrong.site=wrong.site2,iso.rec=iso.rec2)%>%
  # add any nest locations in the veg database to the nest data
  filter((id%in%nests_clean3$id))%>%
  st_drop_geometry()

new_nests_from_veg <-plots2%>%
  #replace old edit flags with new flags
  select(id,out.bounds=out.bounds2,wrong.site=wrong.site2,iso.rec=iso.rec2,type,Site,State,Year)%>%
  # add any nest locations in the veg database to the nest data
  filter(type=="Nest")%>%
  select(-type)%>%
  #remove random plot labeled as nest and nests already in nest data
  filter((id!="ID19RAND027")&!(id%in%nests_clean3$id))%>%
  #add coordinate info
  mutate(data.source="Demo_Veg",
         Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         coord.system="Lat/Long(DD)",
         site.code=substr(id,1,2),
         Species=substr(id,5,8),
         Notes="Nest location record added from demo vegetation database. 01/19/25 EF.")%>%
  st_drop_geometry()
  
nests_clean4<-nests_clean3%>%
  select(-starts_with(c("out.bounds","wrong.site","iso.rec")))%>%
  #merge error flag updates
  left_join(nest_updates,by=c("id"))%>%
  #add new nest location records from veg data (5)
  full_join(new_nests_from_veg,by=c("id","data.source","site.code","Species","Notes","Long","Lat","Site","Year","State","coord.system","out.bounds","wrong.site","iso.rec"))%>%
  mutate_at(vars(starts_with(c("iso","out","wrong","missing","batch"))), ~replace_na(., 0))




# 7. Write final files
#-------------------------------------------------------------------------------------------------

# read in original veg data variables
veg_orig<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select(-c("Site","UTM.Zone","Easting","Northing","Lat","Long","Coordinate.System"))%>%
  rename(veg.id=VegPointID,type=PointType)%>%
  mutate(type=case_when(
    grepl("R",type)~"Random",
    grepl("N",type)~"Nest"))%>%
  distinct(.keep_all = T)
veg_orig<-veg_orig[!duplicated(veg_orig[c('veg.id', 'type')]),]
  #add a way to distinguish random and nest surveys for those that share the same ID
veg_orig<-veg_orig%>%
  arrange(type)%>%
  mutate(veg.id2=ifelse(duplicated(veg.id)&type=="Random",paste0(veg.id,"_rand"),veg.id))%>%
  select(-type)
#join new coordinates to veg info
output_csv<-dplyr::select(veg_final,veg.id2,veg.id,date,type,Site,site.code,utm.zone,missing.coords,out.bounds,wrong.site,iso.rec,
                          Easting,Northing,Lat,Long,data.source,State,Year,Notes_Emily)%>%
  left_join(veg_orig,by=c("veg.id2","veg.id","date"="SurveyDate"))%>%
  distinct(veg.id2,.keep_all = T)


#join veg info to shapefile
output_shp<-plots2%>%
  select(veg.id2=id)%>%
  left_join(output_csv,by=c("veg.id2"))%>%
  filter(veg.id2%in%plots_veg$id)


#write veg files
st_write(output_shp,
         paste0(path_out,"Final_outputs/Veg_locations/veg_locations_01_29_25.shp"), delete_layer =T)
st_write(output_shp,
         paste0(path_out,"Final_outputs/Veg_locations/veg_locations_01_29_25.kml"), delete_layer =T)

write.csv(output_csv,paste0(path_out,"Final_outputs/Veg_locations/corrected_veg_coords_01_29_25.csv"),row.names = F)




  # adjusted nest data csv
write.csv(nests_clean4,paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"),row.names = F)

  #convert to updated nest data shapefile
nests_clean_plots <- plots2%>%
  filter(id%in%nests_clean4$id)%>%
  select(id)%>%
  left_join(nests_clean4,by="id")

st_write(nests_clean_plots,
         paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_29_25.shp"), delete_layer =T)
st_write(nests_clean_plots,
         paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_29_25.kml"), delete_layer =T)
 
