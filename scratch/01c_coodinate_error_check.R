
#Data tidying
library(tidyverse)

#spatial analysis
library(sf)
library(sp)
library(terra)

########################################################
# Review nest and vegetation coordinate fixes in previous scripts
########################################################


## Set file path to data
# -------------------------------------------
path_out<-"D:/Nest_Models/Outputs/"
dat_path<-"D:/Nest_Models/Data/"


## Load data
# -------------------------------------------
nests<- read.csv(paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"))%>%
  dplyr::select(id,site.code,Year,Easting2=Easting,Northing2=Northing,Lat2=Lat,Long2=Long,
                missing.location.rec,missing.site.info,
                missing.coords,State,batch.move.cols=batch_move_DD_to_LatLong,batch.add.dec=batch_dec_addition,
                batch.DecMin=batch_DecMin_DD_reversed,batch.replace=replace_dat,batch.offset=coord_shift,
                out.bounds,wrong.site,iso.rec
                )%>%
  mutate(add.nest.data = 0,
         dataset="Nests")

veg<-read.csv(paste0(path_out,"Final_outputs/Veg_locations/corrected_veg_coords_01_29_25.csv"))%>%
  dplyr::select(type,id=veg.id,site.code,Easting2=Easting,Northing2=Northing,Lat2=Lat,Long2=Long,
                missing.coords,out.bounds,wrong.site,iso.rec)%>%
  mutate(dataset=ifelse(type=="Nest","Nests","Vegetation"),
         missing.location.rec=0,
         missing.site.info=0)%>%
  filter(!(id%in%nests$id))%>%#gets rid of most add.nest.data flags
  dplyr::select(-type)
         
dat1<-rbind(nests,veg)

nrow(veg[veg$dataset=="Vegetation",])


#how many SALS nests? in CT?
nests_og<-read.csv(paste0(dat_path,"Demographic Database/Nests_2002-2024.csv"),na.strings=c("","NOT REC","NA"))%>%
  dplyr::select("id"="SHARPNestID","site.code"="Site", "Year", "Species",
                "coord.system"="Coordinate.System", "utm.zone"="UTM.Zone", "Easting", "Northing", "Lat", "Long")%>%
  #Remove records missing site and year info (these were added as filler data to merge with veg data)
  filter(!is.na(site.code)&!is.na(Year))%>%
  filter(Species=="SALS")
nrow(nests_og)
nests%>%left_join(nests_og%>%dplyr::select(id,Species),by="id")%>%
  filter(Species=="SALS")%>%
  summarise(n=n(),
            error=sum(coord.typo),
            missing.lat=sum(missing.coords))
nests%>%left_join(nests_og%>%dplyr::select(id,Species),by="id")%>%
  filter(Species=="SALS"&State=="CT")%>%
  summarise(n=n(),
            error=sum(coord.typo),
            missing.lat=sum(missing.coords))


## Remaining Error Tally and Summary Plots
#--------------------------------------------------
sum1<-group_by(dat1,site.code,Year,dataset)%>%
  summarise(n_typo=sum(coord.typo,na.rm=T),
            total_n=n(),
            ratio=round(n_typo/total_n,2))%>%
  ungroup()%>%
  group_by(n_typo,ratio)%>%
  mutate(count=n())

ggplot(sum1%>%distinct(n_typo,ratio,count)%>%filter(n_typo!=0),aes(x=n_typo,y=ratio,size=count))+
  geom_point()+
  labs(size="Number of\nSites*Year",title ="Remaining Observations Plotting Outside Marsh Area",x="Remaining Observations with Plotting Errors",y="Plotting Errors:Total Observations at each Site*Year")+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)

ggsave(filename="remaining_errors_tally.png",path=paste0(dat_path,"Coordinate Edit Figures"),
       width=5,height=5,units="in",dpi=400)
# potential batch edits:
# Nest data
# CI 2009 (6), GI 2009 (4) 100% typos - Note: GI 2003 and 07 plot correctly, CI was only sampled in 2009 and has no correct records for reference. Lat values seem too high for GI 2009... check sheets?
# Veg data 
# JC 2014 (Longitude is incorrect)

if(!file.exists(paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/Nest_Veg_data_status_01_3_23.csv"))){
write.csv(sum1,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/Nest_Veg_data_status_01_3_23.csv"),row.names = F)
}
## Filter for those records
#typos<-nests%>%
#  filter(site.code %in% c("MW","AT","OC") & Year %in% c(2014,2015))#%>%
  #dplyr::select(id, site.code,Year, Easting, Northing, Lat, Long, Coordinate.System)




## Edits tally and summary plots
# ---------------------------------
sum2<-summarise(group_by(dat1,dataset),
             missing.coords=sum(missing.coords,na.rm = T),
             batch.DecMin=sum(batch.DecMin,na.rm = T),
             batch.replace=sum(batch.replace,na.rm = T),
             batch.add.dec=sum(batch.add.dec,na.rm = T),
             batch.offset=sum(batch.offset,na.rm = T),
             batch.move.cols = sum(batch.move.cols,na.rm = T),
             coord.typo=sum(coord.typo,na.rm=T),
             add.nest.data=sum(add.nest.data,na.rm=T))

#Combined veg and nest edits
sum3<-summarise(group_by(dat1,site.code,Year),
             missing.coords=sum(missing.coords,na.rm = T),
             batch.DecMin=sum(batch.DecMin,na.rm = T),
             batch.replace=sum(batch.replace,na.rm = T),
             batch.add.dec=sum(batch.add.dec,na.rm = T),
             batch.offset=sum(batch.offset,na.rm = T),
             batch.move.cols = sum(batch.move.cols,na.rm = T),
             coord.typo=sum(coord.typo,na.rm=T),
             add.nest.data=sum(add.nest.data,na.rm=T))%>%
  ungroup()%>%
  mutate(total=rowSums(across(-c(1:2))))%>%
  filter(total>14)%>%
  dplyr::select(-total,-add.nest.data,-coord.typo)

length(unique(sum3$site.code))

ggplot(sum3%>%pivot_longer(-c(1:2),names_to = "edit_type", values_to = "total_records"),aes(x=site.code,y=total_records,fill=as.factor(Year)))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(legend.position = "bottom")+
  facet_wrap(~edit_type,scales = "fixed", ncol = 3)+
  labs(fill="Year",x="Site Code",y= "Total Records",title="Vegetation and Nest Plot Edits")

# Just nest edits
sum4<-summarise(group_by(filter(dat1,dataset=="Nests"),site.code,Year),
                Missing.Coordinates=sum(missing.coords,na.rm = T),
                Decimal.Minute.Conversion=sum(batch.DecMin,na.rm = T),
                Replace.With.Original.Data=sum(batch.replace,na.rm = T),
                Add.Decimal.toDecDegrees=sum(batch.add.dec,na.rm = T),
                Add.Offset.toDecDegrees=sum(batch.offset,na.rm = T),
                Switch.UTM.and.Degree.Columns = sum(batch.move.cols,na.rm = T))%>%
  ungroup()%>%
  mutate(total=rowSums(across(-c(1:2))))%>%
  filter(total>14)%>%
  dplyr::select(-total)

ggplot(sum4%>%pivot_longer(-c(1:2),names_to = "edit_type", values_to = "total_records"),aes(x=site.code,y=total_records,fill=as.factor(Year)))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom")+
  facet_wrap(~edit_type,scales = "fixed", ncol = 3)+
  labs(fill="Year",x="Site Code",y= "Total Records",title="Nest Plot Edits")

ggsave(filename="edit_tally_nests.png",path=paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/Coordinate Edit Figures"),
       width=9,height=6,units="in",dpi=400)

# Just random veg plot edits
sum5<-summarise(group_by(filter(dat1,dataset=="Vegetation"),site.code,Year),
                Missing.Coordinates=sum(missing.coords,na.rm = T),
                Decimal.Minute.Conversion=sum(batch.DecMin,na.rm = T),
                Replace.With.Original.Data=sum(batch.replace,na.rm = T),
                Add.Decimal.toDecDegrees=sum(batch.add.dec,na.rm = T),
                Add.Offset.toDecDegrees=sum(batch.offset,na.rm = T),
                Switch.UTM.and.Degree.Columns = sum(batch.move.cols,na.rm = T))%>%
  ungroup()%>%
  mutate(total=rowSums(across(-c(1:2))))%>%
  filter(total>14)%>%
  dplyr::select(-total)


ggplot(sum5%>%pivot_longer(-c(1:2),names_to = "edit_type", values_to = "total_records"),aes(x=site.code,y=total_records,fill=as.factor(Year)))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom")+
  facet_wrap(~edit_type,scales = "fixed", ncol = 3)+
  labs(fill="Year",x="Site Code",y= "Total Records",title="Random Vegetation Plot Edits")

ggsave(filename="edit_tally_veg.png",path=paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/Coordinate Edit Figures"),
       width=9,height=6,units="in",dpi=400)
#Which records are not plotting out of the east coast area but in the wrong site?
# SC11STSP012, JO13STSP056, 	JO13WILL182, JO12SALS049 plotting at NO in Maine
#	NO12SALS139 ,JO13HYBR103,JO13SALS102 plotting in SC in Maine
# JC12SALS101,JC12SALS102, JC13SALS001 plotting in PC in Rhode Island
# BA11ABDU, PA11SALS008, HM11SALS028, HM11SESP004, HM12SESP001,ER11SESP006, ER12SALS006,ER12SESP006,ER12WILL011,ER12WILL019,ER15SALS418 plotting in BI in CT
# BI12SALS032, BI15SESP110,BI15SALS400,ER19088,ER11SALS029,ER11WILL012,ER13SALS048,ER14CLRA308,ER14SALS133,ER15SALS500,ER15SALS114 are plotting HM in CT
# HM14CLRA315,HM14SALS132,HM15SALS512 are plotting in ER in CT
# SA12SALS002 plotting in ID (Idlewild Park) in LI
# OC14SESP332, OC14WILL220 plotting near MW in NJ





## Plot distribution of error distances from the marsh area
# ---------------------------------------------------------------
library(tmap)
# spatial datasets
library(spData)
library(geodata)
library(rnaturalearth)
# coordinate system -NAD83 / UTM zone 18N
coord<- "EPSG:26918"

# load nest and veg points
nest_pts<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  dplyr::select(id,site_cd,Year,crd_typ)

veg_pts<-st_read(paste0(path_out,"Final_outputs/Veg_locations/veg_locations_12_29_22.shp"))%>%
  #use just the random veg points since nests are already covered in other dataset
  filter(PontTyp=="Random")%>%
  #need to use id and date for veg data identification since some plots were surveyed multiple times in a year
  dplyr::select(veg_id,Site,date,crd_typ)

# load UVVR marsh raster layer
marsh<-rast(paste0(dat_path,"Environmental Predictors/UVVR/UVVR_annual_mean/uvvr_mean_utm18_2.tif"))
# set extent to cover all veg and nest points (take the smaller range of min vals and larger max vals)
xmin<-ifelse(ext(veg_pts)$xmin<ext(nest_pts)$xmin,ext(veg_pts)$xmin,ext(nest_pts)$xmin)
ymin<-ifelse(ext(veg_pts)$ymin<ext(nest_pts)$ymin,ext(veg_pts)$ymin,ext(nest_pts)$ymin)
xmax<-ifelse(ext(veg_pts)$xmax>ext(nest_pts)$xmax,ext(veg_pts)$xmax,ext(nest_pts)$xmax)
ymax<-ifelse(ext(veg_pts)$ymax>ext(nest_pts)$ymax,ext(veg_pts)$ymax,ext(nest_pts)$ymax)

ext(marsh)<-c(xmin,xmax,ymin,ymax)

# set marsh area to values of 0
marsh[marsh>=0]<-0
marsh[is.na(marsh)]<-1

# create a raster distance surface from marsh area boundary
marsh_dst<-terra::gridDistance(marsh, target=0)

# ID which marsh raster cells each plot point falls in
  # turn points into coordinates
vg_crds<-st_coordinates(veg_pts)
nest_crds<-st_coordinates(nest_pts)

  # find unique cell locations for each coordinate point
vgcell<- cellFromXY(marsh,vg_crds)
nestcell<- cellFromXY(marsh,nest_crds)

  # mark the cell locations 

# Or just load in distances for each point location generated in arc
# a distance grid using euclidean distance and res of 100m from the uvvr layer
# extract values to points tool
nest_dst<-st_read(paste0(path_out,"Intermediate_outputs/Nest_locations/marsh_dist_nest.shp"))%>%
  dplyr::select(id,site_cd,crd_typ,Distance=RASTERVALU)%>%
  filter(crd_typ==1)%>%
  #points outside the extent of the UVVR layer, set to max distance
  mutate(Distance=ifelse(Distance==-9999,max(Distance),Distance))%>%
  mutate(Dist_Bin=case_when(
    Distance>=0&Distance<=100~"0-100",
    Distance>100&Distance<=200~"101-200",
    Distance>200&Distance<=400~"201-400",
    Distance>400&Distance<=800~"401-800",
    Distance>800&Distance<=1600~"801-1600",
    Distance>1600&Distance<=3200~"1601-3200",
    Distance>3200~"3201+"
  ),
  Dataset="Nest Plot")

veg_dst<-st_read(paste0(path_out,"Intermediate_outputs/Nest_locations/marsh_dist_veg.shp"))%>%
  #use just the random veg points since nests are already covered in other dataset
  filter(PontTyp=="Random"&crd_typ==1)%>%
  #need to use id and date for veg data identification since some plots were surveyed multiple times in a year
  dplyr::select(id=veg_id,site_cd=Site,crd_typ,Distance=RASTERVALU)%>%
  mutate(Distance=ifelse(Distance==-9999,max(Distance),Distance))%>%
  mutate(Dist_Bin=case_when(
    Distance>=0&Distance<=100~"0-100",
    Distance>100&Distance<=200~"101-200",
    Distance>200&Distance<=400~"201-400",
    Distance>400&Distance<=800~"401-800",
    Distance>800&Distance<=1600~"801-1600",
    Distance>1600&Distance<=3200~"1601-3200",
    Distance>3200~"3201+"
  ),
  Dataset="Random Veg Plot")

#merge two datasets
dat2<-rbind(nest_dst,veg_dst)

order<-c("0-100","101-200","201-400","401-800","801-1600","1601-3200","3201+")
# Histogram of error point distances from marsh area
ggplot(dat2, aes(x = Dist_Bin,fill=Dataset)) + 
  geom_histogram(stat = "count") +
  theme_classic(base_size = 12)+
  scale_x_discrete(labels=order)+
  labs(x="Distance from Marsh Border (m)", y="Number of Observations",fill="Dataset")

ggsave(filename="error_distances.png",path=paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/Coordinate Edit Figures"),
       width=9,height=5,units="in",dpi=400)
