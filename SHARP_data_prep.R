library(tidyverse)
library(lubridate)
library(stringi)
library(sf)
library(raster)

fates<-read.csv("Data/NestFates_2001-2020.csv",na.strings=c("","NOT REC","NA"))
nests<-read.csv("Data/Nests_2001-2020.csv",na.strings=c("","NOT REC","NA"))
veg<-read.csv("Data/Veg_2011-2020.csv",na.strings=c("","NOT REC","NA"))


# 1. A. Format interesting variables
#--------------------------------------------------------------------------------------------------------------------------------------------
dat<-fates%>%
  dplyr::select("SHARPNestID",
                #Binary species response (fail or successful nest)
                "UltimateNestFate",
                
                #continuous species response (or surival rate)
                'MaxNumEggs',"NumFledged","MaxNumChicks",
                
                #Nesting date (temporal variable)
                "EstFirstEggDate","EstHatchDate")%>%
  
  left_join(nests[,c(1:13)],by="SHARPNestID")%>% #right_join to include nest observations in nests that are not in fates
  #could also look at egg vs chick mortality causes
  mutate(incu_start=day(mdy(EstFirstEggDate)),
         incu_end=day(mdy(EstHatchDate)),
         
         #incubation time
         incu_time=incu_end-incu_start,
         NumFledged=as.numeric(NumFledged),
         MaxNumEggs=as.numeric(MaxNumEggs),
         
         #survival rate
         surival_rate=round(NumFledged/MaxNumEggs,2),
         
         #fill in missing site, year, nest IDs based on SHARPNestID (this is the order the nest was found at a site in a given year)
         NestID=case_when(
           is.na(NestID)~as.numeric(stri_sub(SHARPNestID, -3)),
           !(is.na(NestID))~as.numeric(NestID)),
         Site=case_when(
           is.na(Site)~as.character(stri_sub(SHARPNestID, 3,-4)),
           !(is.na(Site))~as.character(Site)),
         Year=case_when(
           is.na(Year)~as.numeric(stri_sub(SHARPNestID, 1,2)),
           !(is.na(Year))~as.numeric(Year)),
         Species=case_when(
           is.na(Species)&length(SHARPNestID)==11~as.character(stri_sub(SHARPNestID, 5,8)),
           !(is.na(Species))~as.character(Species))
         )
#warnings come from missing date, egg, fledgling data, filled in as NA

# set survival rate to 0 for nests with eggs but no fledgling data
dat$surival_rate[is.na(dat$surival_rate)|dat$surival_rate=="NaN"]<-0
dat$surival_rate[dat$MaxNumEggs==0]<-NA


# 1. B. clean up species and spatial data
#--------------------------------------------------------------------------------------------------------------------
  # Complete species surveys (surveys for all spp) were not conducted
  # focus on SALS, SESP, WILL, CLRA, and maybe NESP (complicated due to hybrids and not at most sites).  
  # STSP in northern sites – those will be likely SALSxNESP hybrids and you may want to just ignore them for now (or ignore sites that have them)
  focal_spp<-c("SALS", "SESP", "WILL", "CLRA")
  #filter for focal species (removes 1,408 observations,nrow(dat)-nrow(dat2))
dat2<-filter(dat,Species%in%focal_spp)

  # only use nests with location data - at least 2 UTM coordinates or 2 latlong coordinates (removes 426 original nest observations of all species, 378 observations of focal spp, nrow(dat2)-nrow(dat3))
dat3<-filter(dat,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))
  # remove coordinates that are 0 or small values (less than 10, likely typos)
dat4<-filter(dat3,if_any(c(Easting,Northing,Lat,Long), ~ .>10))






# 2. Some data summaries
#--------------------------------------------------------------------------------------------------------------------------------------------
#number of sites and species across years
length(unique(dat4$Site)) #83 sites
length(unique(dat4$Species)) #4 species 

# max, min, mean number of nests found at each site per year. Number of years each site was surveyed. 
nests_per_site<-dat4%>%
  group_by(Year,Site)%>%
  summarise(n_nests=n())%>%
  ungroup()%>%
    group_by(Site)%>%
    summarise(max_nest_peryr=max(n_nests),
              min_nest_peryr=min(n_nests),
              mean_nest_peryr=round(mean(n_nests)),
              n_yrs_w_nests=n())

species_obs<-dat4%>%
  group_by(Species)%>%
  # number of sites each species was found nesting
  summarise(n_sites=length(unique(Site)),
  # number of years each species was surveyed and found
            n_years = length(unique(Year)))

#number of nests for each species each site and year
species_site_yr<-dat4%>%
  ungroup()%>%
  group_by(Species,Site,Year)%>%
  summarise(n_nests=n())




# 3. format spatial data for nest locations 
#--------------------------------------------------------------------------------------------------------------------------------------------

# get UTM zones for sites with zones reported for at least some data
utm <- nests %>% 
  filter(!is.na(UTM.Zone)) %>% 
  group_by(Site) %>% 
  summarize(UTM.Zone = as.factor(tail(UTM.Zone, 1)))%>%
  ungroup()

#range of longitude (E-W) for eastern seaboard should be ~66-99 deg, latitude should be 32-48. degrees outside this range are likely typos
dat_spatial<-dat4%>%
         #if Lat is missing and Easting is using degrees (2 digit integers) in Latitude (between 20 and 50), fill in Lat from the Easting column
  mutate(Lat=ifelse(abs(Easting)<50&abs(Easting)>20&is.na(Lat),abs(Easting),Lat),
         #if Lat is missing and Northing is using degrees (2 digit integers) in Latitude (between 20 and 50), fill in Lat from the Northing column
         Lat=ifelse(abs(Northing)<50&abs(Northing)>20&is.na(Lat),abs(Northing),Lat),
         #if Long is missing and Easting is using degrees (2 digit integers) in Longitude (between -90 and -60), fill in Long from the Easting column
         Long=ifelse(abs(Easting)<90&abs(Easting)>60&is.na(Long),Easting,Long),
         #if Long is missing and Northing is using degrees (2 digit integers) in Longitude (between -90 and -60), fill in Long from the Northing column
         Long=ifelse(abs(Northing)<90&abs(Northing)>60&is.na(Long),Northing,Long),
         Long=ifelse(Long>0,-Long,Long),
         #if Lat Long is available, remove UTM coords and use the Lat Long instead
         Easting=ifelse(is.na(Lat),Easting,NA),
         Northing=ifelse(is.na(Long),Northing,NA),
         #label whether observation uses UTM or lat long
         Coordinate.System=ifelse(is.na(Lat), "UTM",'Lat/Long(DD)'),
         UTM.Zone=as.factor(UTM.Zone))%>%
  #fill in missing UTM zones for observations at sites with recorded zones
  left_join(utm,by=c("Site"))%>%
  mutate(UTM.Zone=coalesce(UTM.Zone.x,UTM.Zone.y))%>%
  dplyr::select(-UTM.Zone.x,-UTM.Zone.y)%>%
  #some coordinates are really small and seem out of range of Eastern North America, probably more typos, remove here
  filter(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.)))



# get coordinate systems and convert coordinates to spatial data
unique(utm$UTM.Zone) #17, 18, 19

utm17<- "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
utm18<- "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
utm19<- "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

template<-raster("Data/prediction grids/raster_test_land.tif")

#convert to spatial points
plots_utm17 <- st_as_sf(filter(dat_spatial,UTM.Zone=="17T"&Coordinate.System=="UTM"), coords = c("Easting", "Northing"), crs = utm17)%>%
  st_transform(template@crs)
plots_utm18 <- st_as_sf(filter(dat_spatial,UTM.Zone=="18T"&Coordinate.System=="UTM"), coords = c("Easting", "Northing"), crs = utm18)%>%
  st_transform(template@crs)
plots_utm19 <- st_as_sf(filter(dat_spatial,UTM.Zone=="19T"&Coordinate.System=="UTM"), coords = c("Easting", "Northing"), crs = utm19)%>%
  st_transform(template@crs)
plots_latlong <- st_as_sf(filter(dat_spatial,Coordinate.System=="Lat/Long(DD)"), coords = c("Lat", "Long"), crs = utm18)%>%
  st_transform(template@crs)%>%
  rename(Lat=Easting,Long=Northing)

plots<-rbind(plots_utm17,plots_utm18,plots_utm19,plots_latlong)

#create 50m buffer around nest sites
#plots_buff<-st_buffer(plots,50)

#admin boundaries
states<-st_read("Data/geographies/tl_2018_us_state.shp")%>%
  #sharp states
  filter(POSTAL%in%c("ME","NH","MA","RI","CT","NY","NJ","DE","MD","VA","VT","PA"))%>%
  st_transform(template@crs)
bbox<-st_bbox(states)

#plot nest sites
ggplot() +
  geom_sf(data=states,color="black") +
  geom_sf(data=plots,color="red",size=3)+
  theme_bw()+
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax), expand = FALSE)



#outputs
#nest points shp
st_write(plots,"Outputs/nests_with_reasonable_coords.shp")

#observations with coordinate issues
fates_errors<-dat[!(fates$SHARPNestID%in%plots$SHARPNestID),]%>%mutate(error_coords=1,reasonable_coords=0)
#observations with reasonable or edited coordinates
fates_edited<-dat_spatial[(dat_spatial$SHARPNestID%in%plots$SHARPNestID),]%>%mutate(error_coords=0,reasonable_coords=1)

nrow(fates_errors)+nrow(fates_edited)
nrow(fates)

write.csv(rbind(fates_errors,fates_edited),"Outputs/fates_with_spatial_edits_Feng.csv",row.names = F)



# 4. Set up species response data
#--------------------------------------------------------------------------------------------------------------------------------------------

# two types of responses: nest success and nesting presence

# Cannot treat missing sites/species observations per year as absences. These missing data were not surveyed. 
# Use random veg plots for these absences



# nest success
  # Look at possible nest fates
unique(dat_spatial$UltimateNestFate)

success<-dat_spatial%>%
  #remove uncertain successes or fails
  filter(UltimateNestFate%in%c("UNKNOWN IF FLEDGED OR FAILED","UNKNOWN"))%>%
  mutate(success=ifelse(UltimateNestFate%in%c("INACTIVE","DEPREDATED","FLOODED","FAIL UNKNOWN")))
  

#download a phylogeny for your species from birdtree.org
