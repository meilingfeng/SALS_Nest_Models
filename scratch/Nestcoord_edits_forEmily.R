
# only use nests with location data - at least 2 UTM coordinates or 2 latlong coordinates 
all_5 <-filter(all_4,(if_all(c(Easting,Northing), ~ !is.na(.))|if_all(c(Lat,Long), ~ !is.na(.))))

# remove coordinates that are 0 or small values (less than 10, likely typos)
all_6 <-filter(all_5,if_any(c(Easting,Northing,Lat,Long), ~ .>10))

## EDIT TO TAKE OUT DELAWARE LISA NESTS TOO
no_VA_DE <- filter(all_6, Site.x != "SG" & Site.x != "WI" & Site.x != "216138" & Site.x != "WOBE" & Site.x != "MIRI" 
                   & Site.x != "STJN")

VA_DE <- filter(all_6, Site.x == "SG" | Site.x == "WI" | Site.x == "WOBE" | Site.x == "MIRI" 
                | Site.x == "STJN")

# trying again with fixing the NJ sites in 2014 and 2015 
#starts with 393 and 394

dat_spatial <- no_VA_DE %>%
  mutate(
    #if Lat is missing and Northing is using degrees (2 digit integers) in Latitude (between 20 and 50), fill in Lat from the Northing column
    Lat=ifelse(abs(Northing)<50&abs(Northing)>20&is.na(Lat),abs(Northing),Lat),
    #if Long is missing and Easting is using degrees (2 digit integers) in Longitude (between -90 and -60), fill in Long from the Easting column
    Long=ifelse(abs(Easting)<90&abs(Easting)>60&is.na(Long),Easting,Long),
    #if Long is missing and Northing is using degrees (2 digit integers) in Longitude (between -90 and -60), fill in Long from the Northing column
    Long=ifelse(abs(Northing)<90&abs(Northing)>60&is.na(Long),Northing,Long),
    Long=ifelse(Long>0,-Long,Long),
    Long = abs(ifelse(startsWith(as.character(Easting), "72"), abs(Easting), Long)),
    Long = abs(ifelse(startsWith(as.character(Easting), "71"), abs(Easting), Long)),
    Long = abs(ifelse(startsWith(as.character(Easting), "74"), abs(Easting), Long)),
    Long = abs(ifelse(startsWith(as.character(Northing), "72"), abs(Northing), Long)),
    Long = abs(ifelse(startsWith(as.character(Northing), "71"), abs(Northing), Long)),
    Long = abs(ifelse(startsWith(as.character(Northing), "74"), abs(Northing), Long)),
    Lat = abs(ifelse(startsWith(as.character(Easting), "41"), abs(Easting), Lat)),
    Lat = abs(ifelse(startsWith(as.character(Easting), "393"), abs(Easting), Lat)),
    Lat = abs(ifelse(startsWith(as.character(Easting), "394"), abs(Easting), Lat)),
    Lat = abs(ifelse(startsWith(as.character(Northing), "41"), abs(Northing), Lat)),
    Lat = abs(ifelse(startsWith(as.character(Northing), "393"), abs(Northing), Lat)),
    Lat = abs(ifelse(startsWith(as.character(Northing), "394"), abs(Northing), Lat)),
    ## removing from easting/northing
    Easting = ifelse(startsWith(as.character(Easting), "72"), NA, Easting),
    Easting = ifelse(startsWith(as.character(Easting), "71"), NA, Easting),
    Easting = ifelse(startsWith(as.character(Easting), "74"), NA, Easting),
    Northing = ifelse(startsWith(as.character(Easting), "41"), NA, Northing),
    Northing = ifelse(startsWith(as.character(Easting), "393"), NA, Northing),
    Northing = ifelse(startsWith(as.character(Easting), "394"), NA, Northing),
    #if Lat Long is available, remove UTM coords and use the Lat Long instead
    Easting=ifelse(is.na(Lat),Easting,NA),
    Northing=ifelse(is.na(Long),Northing,NA),
    #label whether observation uses UTM or lat long
    Coordinate.System=ifelse(is.na(Lat), "UTM",'Lat/Long(DD)')) 

# problems from above

## HM15SALS001
# all SG's and WI's

#SG21039

# remove all of the below
partial_gpscoord <- dat_spatial %>%
  filter(is.na(Lat) & is.na(Long))

# example BI14SALS211
test <- dat_spatial %>%
  filter(startsWith(as.character(Long), "72.") |
           startsWith(as.character(Long), "71.") | 
           startsWith(as.character(Long), "70.") |
           startsWith(as.character(Long), "74.")) %>%
  mutate(Long = Long*(-1))

test2 <- dat_spatial %>%
  mutate(Long = abs(Long)) %>%
  filter(startsWith(as.character(Long), "72") & !startsWith(as.character(Long), "72.")
         | (startsWith(as.character(Long), "71") & !startsWith(as.character(Long), "71."))
         | (startsWith(as.character(Long), "74") & !startsWith(as.character(Long), "74."))) %>%
  mutate(Lat = Lat*(0.00001)) %>%
  filter(Long > 7000000)  %>%
  mutate(Long = Long*(0.00001)) %>%
  mutate(Long = Long*(-1))

test3 <- dat_spatial %>%
  mutate(Long = abs(Long)) %>%
  filter(startsWith(as.character(Long), "72") & !startsWith(as.character(Long), "72.")
         | (startsWith(as.character(Long), "71") & !startsWith(as.character(Long), "71.")
            | (startsWith(as.character(Long), "74") & !startsWith(as.character(Long), "74.")))) %>%
  mutate(Lat = Lat*(0.00001)) %>%
  filter(Long < 7000000 & Long > 700000)  %>%
  mutate(Long = Long*(0.0001)) %>%
  mutate(Long = Long*(-1))

test4 <- dat_spatial %>%
  mutate(Long = abs(Long)) %>%
  filter(startsWith(as.character(Long), "72") & !startsWith(as.character(Long), "72.")
         | (startsWith(as.character(Long), "71") & !startsWith(as.character(Long), "71.")
            | (startsWith(as.character(Long), "74") & !startsWith(as.character(Long), "74.")))) %>%
  mutate(Lat = Lat*(0.00001)) %>%
  filter(Long < 700000 & Long > 70000)  %>%
  mutate(Long = Long*(0.001)) %>%
  mutate(Long = Long*(-1))

combine_1 <- rbind(test, test2)
combine_2 <- rbind(combine_1, test3)
combine_3 <- rbind(combine_2, test4)

## recombine those and the va sites 

final_maybe_1 <- rbind(partial_gpscoord, VA_DE)
final_maybe_2 <- rbind(final_maybe_1, combine_2)
final_maybe_3 <- rbind(final_maybe_2, combine_3) # find problem with maxnumeggs, all, and 

NJ_test_AT <- final_maybe_3 %>%
  filter(Site.x == "AT") %>%
  filter(startsWith(as.character(Long), "-74.")) %>%
  mutate(Lat = Lat + mean(testing_NJ_AT$Lat_diff, na.rm = TRUE)) %>%
  mutate(Long = Long + mean(testing_NJ_AT$Long_diff, na.rm = TRUE)) 

nj_testing_next_AT <- st_as_sf(NJ_test_AT, coords = c("Long", "Lat"), crs = 4326)
mapview::mapview(nj_testing_next_AT, layer.name = "Nests", cex = 5, zcol = "Species")

NJ_test_OC <- final_maybe_3 %>%
  filter(Site.x == "OC") %>%
  filter(startsWith(as.character(Long), "-74.")) %>%
  mutate(Lat = Lat + mean(testing_NJ_OC$Lat_diff, na.rm = TRUE)) %>%
  mutate(Long = Long + mean(testing_NJ_OC$Long_diff, na.rm = TRUE)) 

NJ_testing_next_OC <- st_as_sf(NJ_test_OC, coords = c("Long", "Lat"), crs = 4326)
mapview::mapview(NJ_testing_next_OC, layer.name = "Nests", cex = 5, zcol = "Species")

NJ_test_MW <- final_maybe_3 %>%
  filter(Site.x == "MW") %>%
  filter(startsWith(as.character(Long), "-74.")) %>%
  mutate(Lat = Lat + mean(testing_NJ_MW$Lat_diff, na.rm = TRUE)) %>%
  mutate(Long = Long + mean(testing_NJ_MW$Long_diff, na.rm = TRUE)) %>%
  filter(!is.na(Lat))

nj_testing_next_MW <- st_as_sf(NJ_test_MW, coords = c("Long", "Lat"), crs = 4326)
mapview::mapview(nj_testing_next_MW, layer.name = "Nests", cex = 5, zcol = "Species")

NJ_test_AT_OC <- rbind(NJ_test_AT, NJ_test_OC)
NJ_test_AT_OC_MW <- rbind(NJ_test_AT_OC, NJ_test_MW)

NJ_test_AT_OC_MW <- NJ_test_AT_OC_MW %>%
  relocate(SHARPNestID, Site.x, Site.y, Year, Species, Easting, Northing, Lat, Long, UltimateNestFate, MaxNumEggs, NumFledged, 
           VisitDate, VisitTime, Eggs, Nestlings, ChickAgeA, ChickAgeB, ChickAgeC, ChickAgeD, ChickAgeE, EggsWarm, Status, 
           Notes.x, Notes.y, Tapping, Starring, Pipping, Chicksseen_live_notes, Dead_hatching_eggs, Dead_chicks, Site, 
           SurveyDate, DistLipToGround, DistBottomToGround, AvgCent, NestLoc, WovenNestCanopy, WovenNestCanopyLiveDead, 
           PercentVisible, SpartinaPatensComp, SpartinaAlternifloraComp, DistichlisSpicataComp, JuncusGerardiiComp, Final_unvegcomp, 
           Final_otherspcomp, Final_vegcovernest, All, Interval)


##### START HERE ######## 

rest_test <- final_maybe_3 %>%
  filter(is.na(Long) | !(startsWith(as.character(Long), "-74")))

rest_test <- rest_test %>%
  relocate(SHARPNestID, Site.x, Site.y, Year, Species, Easting, Northing, Lat, Long, UltimateNestFate, MaxNumEggs, NumFledged, 
           VisitDate, VisitTime, Eggs, Nestlings, ChickAgeA, ChickAgeB, ChickAgeC, ChickAgeD, ChickAgeE, EggsWarm, Status, 
           Notes.x, Notes.y, Tapping, Starring, Pipping, Chicksseen_live_notes, Dead_hatching_eggs, Dead_chicks, Site, 
           SurveyDate, DistLipToGround, DistBottomToGround, AvgCent, NestLoc, WovenNestCanopy, WovenNestCanopyLiveDead, 
           PercentVisible, SpartinaPatensComp, SpartinaAlternifloraComp, DistichlisSpicataComp, JuncusGerardiiComp, Final_unvegcomp, 
           Final_otherspcomp, Final_vegcovernest, All, Interval)

## need to add AvgCent to NJ_test_AT_OC_MW

final_maybe_4 <- rbind(NJ_test_AT_OC_MW, rest_test)

## filter and fix HM, ER, BI's JK don't need to 

#test_ct <- final_maybe_4 %>%
#  filter(Year < 2016) %>%
#  filter(Site.x == "HM" | Site.x == "ER" | Site.x == "BI") %>%
#  mutate(Site.check = ifelse(startsWith(as.character(Long), "-71"), "BI", "Not_BI")) %>%
#  relocate(SHARPNestID, Site.x, Site.check)

# merge AT and ATT, BI and Barn Island, HM and Hammo, OC and Oyster creek
site_use<- c("AT", "BI", "CF", "CL", "DI", "EL", "ER", "HM", "FB", "FS", 
             "ID", "JC", "JO", "LU", "MN", "MQ", "MW", "Morris Island", "NC", "NO", 
             "OC", "PA", "PB", "PR", "SA", "SC", "SG", "SY", "WI", "NC2", 
             "MM", "SP", "FS1", "WA", "MIRI", "STJN", "WOBE")

state <- c("NJ", "CT", "MA", "NH", "ME", "ME", "CT", "CT", "ME", "NY",
           "NY", "RI", "ME", "NH", "NY", "ME", "NJ", "NJ", "NY", "ME", 
           "NJ", "XX", "ME", "MA", "NY", "ME", "VA", "ME", "VA", "NY", 
           "ME", "RI", "NY", "CT", "DE", "DE", "DE")

utm_zone <- c(18, 19, 19, 19, 19, 19, 18, 18, 19, 18, 
              18, 19, 19, 19, 18, 19, 18, 18, 18, 19, 
              18, NA, 19, 19, 18, 19, 18, 19, 18, 19, 
              19, 19, 18, 19, 18, 18, 18)

#WA site is actually BI change the one record

sites_states <- cbind(site_use, state)

sites_states_utm <- cbind(sites_states, utm_zone)

sites_states_utm <- tibble::as_tibble(sites_states_utm)

#before this need to rename sites for CT pre 2016. If 

final_maybe_5 <- rename(final_maybe_4, site_use = Site.x)

final_maybe_6 <- left_join(final_maybe_5, sites_states_utm)

final_maybe_7 <- relocate(final_maybe_6, SHARPNestID, site_use, Year, state, Coordinate.System, utm_zone)

final_maybe_7$utm_zone <- as.factor(final_maybe_7$utm_zone)
