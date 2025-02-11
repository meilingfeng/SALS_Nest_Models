library(tidyverse)
library(sf)


# Original Nest location data with site codes in ids that differ from their recorded site code
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2002-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site","Site"="SiteName", "Year", "Species",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))%>%
  #remove rapid demo records
  filter(!(site.code=="RD" | 
             grepl("^[[:digit:]]{2,}",site.code) | grepl("^[[:digit:]]{2,}",Site) |
             grepl("[[:digit:]]{2,}$",site.code) | grepl("[[:digit:]]{2,}$",Site) |
             substring(id,1,2)=="RD"))%>%
  filter(substr(site.code,1,2)!=substr(id,1,2))%>%
  mutate(type="Nest",
         data.source="demo_nests")%>%
  dplyr::select(id,type,data.source,site.code,Site,Year,Lat,Long,Easting,Northing,utm.zone)


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
             substring(veg.id,1,2)=="RD"))%>%
  filter(substr(site.code,1,2)!=substr(veg.id,1,2))%>%
  filter(!(veg.id%in%nests$id))%>%
  dplyr::select(id=veg.id,type,site.code,Site,Year,Lat,Long,Easting,Northing,utm.zone)%>%
  mutate(data.source="demo_veg")

sitenameconflict<-rbind(nests,veg)


#How many are plotting in sites different from what's in their recorded site code
veg_clean<-read.csv(paste0(path_out,"Final_outputs/Veg_locations/corrected_veg_coords_01_29_25.csv"))

nest_clean<-read.csv(paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"))

wrongsite_veg<-veg_clean[!(veg_clean$veg.id2%in%nest_clean$id)&veg_clean$wrong.site==1,]%>%
  dplyr::select(id=veg.id,type,site.code,Site,Year,Lat,Long,Easting,Northing,utm.zone)%>%
  mutate(data.source="demo_veg")

wrongsite_nest<-nest_clean[nest_clean$wrong.site==1,]%>%
  mutate(type="Nest",
         data.source="demo_nests")%>%
  dplyr::select(id,type,data.source,site.code,Site,Year,Lat,Long,Easting,Northing,utm.zone)

wrongsites<-rbind(wrongsite_veg,wrongsite_nest)

write.csv(sitenameconflict,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/demodata_site_code_recordID_conflict.csv"),row.names = F)

write.csv(wrongsites,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/demodata_site_code_plotting_conflict.csv"),row.names = F)
