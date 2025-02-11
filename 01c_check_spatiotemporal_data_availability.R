
library(sf)
library(tidyverse)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(dismo)



##############################################################################################
#Start Species Specific Analysis Here. Summarize data availability after spatial data cleaning.
##############################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## Set Species
# -------------------------------------
#Record nesting data summary for each species of interest
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")

for (i in 1:length(speciesnames)){

## Table for data summaries
#----------------------------------------------------------------------
data_sum<-data.frame(Step=c(1:10), 
                     Filter=c("Nest locations from demo database",
                              "Remove rapid/intensive demo records",
                               "Remove records missing coordinate/location information",
                               "Add records with coordinate information in the demo veg data",
                               paste0("Select species (",speciesnames[i],")"),
                               "Select years with standard SHARP demo protocol (2010-present)",
                               "Remove records plotting out of Northeast tidal marsh (20m buffer of Correll and Ganju layers)",
                                "Remove records plotting more than 500m away from any other records",
                                "Remove records plotting in a site that's different than what is recorded in the record",
                                "Remove records without fate information (NA, unknown, or inactive)"),
                     N_records_removed=rep(NA,10),
                     N_records=rep(NA,10),
                     #add the total number of sites and years with data availability, post filtering
                     N_sites=rep(NA,10),
                     N_years=rep(NA,10))
                     
## Data filtering steps (summary of scripts 01a and 01b)
#-----------------------------------------------------------------------------
#Original nest data
nests<-read.csv(paste0(dat_path,"Demographic Database/Nests_2002-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site", "Year", "Species",Site=SiteName,
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))

data_sum[1,"N_records"]<-nrow(nests)
data_sum[1,"N_records_removed"]<-0

#Step 1. Remove any sites listed as "RD"- these are rapid demo sites and we want intensive demo only
nests2<-nests%>%
  filter(!(site.code=="RD" | 
                        grepl("^[[:digit:]]{2,}",site.code) | grepl("^[[:digit:]]{2,}",Site) |
                        grepl("[[:digit:]]{2,}$",site.code) | grepl("[[:digit:]]{2,}$",Site) |
                        substring(id,1,2)=="RD"))
  #total records remaining
data_sum[2,"N_records"]<-nrow(nests2)
  #total records removed
data_sum[2,"N_records_removed"]<-nrow(nests)-nrow(nests2)


#Step 2. Remove records missing coordinate/location information
nests3<-read.csv(paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"))%>%
  filter(!grepl("^Nest location record added from demo vegetation database.",Notes))%>%
  filter(missing.coords!=1)
#total records remaining
data_sum[3,"N_records"]<-nrow(nests3)
#total records removed
data_sum[3,"N_records_removed"]<-nrow(nests2)-nrow(nests3)


#Step 3. Add records with coordinate infromation in the demo veg data
nests4<-read.csv(paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"))%>%
  filter(missing.coords!=1)
#total records remaining
data_sum[4,"N_records"]<-nrow(nests4)
#total records removed
data_sum[4,"N_records_removed"]<-nrow(nests3)-nrow(nests4)


#Step 4. Select species
nests5<-nests4%>%
  filter(Species==speciesnames[i])
#total records remaining
data_sum[5,"N_records"]<-nrow(nests5)
#total records removed
data_sum[5,"N_records_removed"]<-nrow(nests4)-nrow(nests5)


#Step 5. Select years with standard SHARP demo protocol (2010-present)
nests6<-nests5%>%
  filter(Year>=2010)
#total records remaining
data_sum[6,"N_records"]<-nrow(nests6)
#total records removed
data_sum[6,"N_records_removed"]<-nrow(nests5)-nrow(nests6)


#Step 6. Remove records plotting out of Northeast tidal marsh (20m buffer of Correll and Ganju layers)
nests7<-nests6%>%
  filter(out.bounds!=1)
#total records remaining
data_sum[7,"N_records"]<-nrow(nests7)
#total records removed
data_sum[7,"N_records_removed"]<-nrow(nests6)-nrow(nests7)


#Step 7. Remove records plotting more than 500m away from any other records
nests8<-nests7%>%
  filter(iso.rec!=1)
#total records remaining
data_sum[8,"N_records"]<-nrow(nests8)
#total records removed
data_sum[8,"N_records_removed"]<-nrow(nests7)-nrow(nests8)


#Step 8. Remove records plotting in a site that's different than what is recorded in the record
nests9<-nests8%>%
  filter(wrong.site!=1)
#total records remaining
data_sum[9,"N_records"]<-nrow(nests9)
#total records removed
data_sum[9,"N_records_removed"]<-nrow(nests8)-nrow(nests9)


#Step 9. Remove records without fate information
nests10<-nests9%>%
  filter(!is.na(fate)&fate%in%c("DEPREDATED","FLEDGED","FAIL UNKNOWN","FLOODED"))#removing inactive nests, those without eggs, and unknown fate
#total records remaining
data_sum[10,"N_records"]<-nrow(nests10)
#total records removed
data_sum[10,"N_records_removed"]<-nrow(nests9)-nrow(nests10)
unique(nests3$fate)

## Summary of data availability
#---------------------------------------------------------
data_sum[10,"N_sites"]<-length(unique(nests9$site.code))
data_sum[10,"N_years"]<-length(unique(nests9$Year))
#  speciesnames nest_n fate_n site_n                                                           year_n
#  1         SALS   2476   2301     42 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  2         SESP   1078   1036     42             2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019
#  3         CLRA    394    325     42 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021
#  4         WILL    357    284     42       2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  5         NESP    273    271     42       2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#  6         HYBR    100     99     42                                     2011, 2012, 2013, 2014, 2015

#write data summary to file
write.csv(data_sum,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/",speciesnames[i],"_valid_nest_locations_2010_2024_filtering_steps.csv"),row.names = F)


## Write filtered, species specific nests to a shapefile
#-------------------------------------------------------------------------------------------------
nest_points<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_29_25.shp"))%>%
    filter(id%in%nests9$id)
  
st_write(nest_points,paste0(path_out,"Final_outputs/Nest_locations/",speciesnames[i],"_valid_nest_locations_2010_2024.shp"),append = FALSE)
}
