library(tidyverse)
library(terra)
library(sf)
library(stringr)
library(vegan)


########################################
# Format Vegetation Characteristics
########################################


### Set up
########################################

## 1. file paths
#------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. Read in vegetation data
#------------------------------------
# There is veg data from 3 different surveys:

# b) rapid vegetation data (measured within a 50m radius circle)
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/Rapid_Vegetation_20230303_updated.csv"))%>%
  rename(id=BirdPtID)%>%
  mutate(year=as.factor(year(mdy(Date))))%>%
  filter(year!="2023")

## 3. Add the location data for the transect and rapid veg points
#--------------------------------------------------------------------

rapid_dat<-read.csv(paste0(dat_path,"Survey Database/All_points_attribute_table_20230426.csv"))%>%
  rename(id=region_ID)%>%
  dplyr::select("id","Long"="POINT_X","Lat"="POINT_Y")%>%
  right_join(rapid_dat,by="id")%>%
  filter(!(is.na(Lat)&is.na(Long)))



## 4. define function to calculate mode
#-----------------------------------------------------
find_mode <- function(x) {
  u <- unique(na.omit(x))
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}





### Format variables for each dataset
#####################################################


## 2. Rapid vegetation data (Survey database)
#--------------------------------------------------
rapid_dat2<-rapid_dat%>%
  ## indicate dataset is rapid_veg
  mutate(data="rapid_veg")%>%
  ## remove rows that have no data (all NAs)
  filter(!(if_all(contains(c("Dom","Sp","Dead","Water","Wrack")), ~is.na(.x))))%>%
  ## remove duplicate records
  distinct(pick(c("id","year","Date","Time")),.keep_all = TRUE)%>%
  ## Select Important Variables:
  dplyr::select(id, Lat, year,Date,SHARPTide,Time,data,
                #Categorical variables:
                # a) Percent cover within 50 m categories for low and high marsh, saltmarsh border, brackish terrestrial border, invasive species, etc
                LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, PannesChannelsCC,UplandCC,WrackCC=Wrack,OpenWaterCC=OpenWater,
                #Continuous variables:
                # b) number of dead snags (dead trees) >1m (DBH or height?) within 50m radius
                DeadSnags,
                #Percent cover of individual species
                # species names
                DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10,
                # species percent cover
                Sp1Percent,Sp2Percent,Sp3Percent,Sp4Percent,Sp5Percent,Sp6Percent,Sp7Percent,Sp8Percent,Sp9Percent,Sp10Percent)%>%
  ## Fix Formatting: 
  # reclassify categorical percent cover variables into 25% intervals
  mutate(across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),
                ~case_when(  # 0 is not present
                  .x=="0"~0,
                  # 0.5 is less than a percent, mark as 1
                  .x=="0.5"~1,
                  # groups together 1,2,3 -> 1-25%, marks as 2
                  .x%in%c("1","2","3")~2,
                  # 4 -> 26-50%, mark as 3
                  # 5 -> 50-75%, mark as 4
                  # 6-> 76-100%, mark as 5
                  .x=="4"~3,
                  .x=="5"~4,
                  .x=="6"~5,
                  .x=="NULL"|is.na(.x)~NA)))%>%
  #factor(LowMarshCC, levels=c('0', '<1', '1-25', '26-50','51-75','76-100'))
  # Aggregate categories for species names
  mutate(across(c(DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10),
                ~case_when(grepl("stichlis",.x,)~"distichlis_pct",
                           grepl("atens",.x)~"patens_pct",
                           grepl("lterniflora-tall",.x)~"alt_tall_pct",
                           grepl("lterniflora-short",.x)~"alt_short_pct",
                           grepl("ragmites",.x)~"phrag_pct",
                           grepl("gerardii",.x)~"gerardii_pct",
                           grepl("lterniflora",.x)&!grepl("lterniflora-tall|lterniflora-short",.x)~"alt_pct",
                           grepl("Conifer|Juniperus|Prunus|Pinus|Quercus|Acer",.x)&.x!="Angiosperm/Conifer shrub" ~"trees_pct",
                           grepl("shrub|Bocconia|Baccharis|Iva|Morella|Myrica|Elaeagnus umbellata",.x)~"shrubs_pct", #Elaeagnus umbellata is invasive
                           grepl("Bare Ground",.x)&.x!="Bare Ground (road/pavement)"~"bare_pct",#includes organic matter, sand, mud
                           grepl("Water|water",.x)&!grepl("ditch/creek|pool/panne",.x)~"open_water_pct", #(impoundment, Lemna is duckweed, usually just in standing water)
                           grepl("ditch/creek|Channel|pool/panne",.x)~"pannes_channels_pct",
                           grepl("pavement|Mowed Grass|impoundment",.x)~"modified_habitat_pct", #also ditches and channels?
                           grepl("Upland",.x)&.x!="Upland (pavement)"~"upland_pct",
                           grepl("Non-saltmarsh herbaceous spp.|Angiosperm|Unknown",.x)&.x!="Angiosperm/Conifer shrub"~"other_sp_pct",#divide into terrestrial and aquatic species? Celastrus orbiculatus and Lythrum salicaria is invasive
                           grepl("Pycnanthemum|Lythrum|Hibiscus|Helenium amarum|Lespedeza|Symphyotrichum|Erechtites|Verbena|Glaux|Phytolacca|Pluchea|Argentina|Tripleurospermum|Suaeda|Apocynum cannabinum|Ptilimnium|Triglochin|Atriplex|Toxicodendron|Impatiens|Amaranthus|Sesuvium|Mikania|Celastrus orbiculatus",.x)~ "other_marsh_sp_pct",
                           grepl("Alisma|latifolia|Sagittaria|Peltandra|Ludwigia|Nuphar|Pontederia|Lemna",.x)~"other_aquatic_sp_pct",
                           grepl("Lythrum|Rubus|Thinopyrum|Elaeagnus umbellata",.x)~ "invasive_pct", #includes herbaceous, grasses, shrubs
                           grepl("Juncus",.x)&!grepl("gerardii",.x)~"rushes_pct", #(juncus roemerianus,arcticus littoralis) , exclude gerardii?
                           grepl("Sedge|Bolboschoenus|Schoenoplectus|Carex|Scirpus|Cyperus|Eleocharis",.x)~"sedges_pct", #(Schoenoplectus pungens,tabernaemontani)
                           grepl("Poaceae sp.|Leersia|Spartina||Poa sp.|Panicum|Elymus|Festuca|Puccinellia|Ammophila|Zizania|Setaria|Echinochloa|Leptochloa",.x)&!grepl("atens|lterniflora",.x)~"grasses_pct", #exclude alterniflora and patens
                           grepl("Spartina",.x)&!grepl("atens|lterniflora",.x)~"cord_grasses_pct", #(pectinata,cynosuroides) exclude alt and patens
                           grepl("Solidago",.x)~"goldenrods_pct", #(sempervirens)
                           grepl("Plantago",.x)~"plantains_pct", #(maritima)
                           grepl("Typha",.x)~"cat_tails_pct", #(angustifolia,latifolia)
                           grepl("Limonium",.x)~"limonium_pct", #(carolinianum,nashii)
                           grepl("Polygonum",.x)~"polygonum_pct", #buckwheat and knotgrasses (pensylvanicum,punctatum,perfoliatum)
                           grepl("Salicornia",.x)~"salicornia_pct", #(depressa,bigelovii)
                           grepl("Wrack",.x)~"wrack_pct",
                           grepl("Thatch",.x)~"thatch_pct",
                           grepl("Algae|S. distichum",.x)~NA)))%>%
  # adjust the missing data values for snag count variable
  mutate(DeadSnags=ifelse(DeadSnags%in%c("NULL","-1"),NA,as.numeric(DeadSnags)))

# Restructure table to have species names as column names and percentages as the values
# first take the columns with species names and species percentages and make them into 1 column for names and 1 column for percents
rapid_dat3<-  pivot_longer(rapid_dat2,cols=contains("Dom"),names_to = "transect", values_to = "species")%>%
  dplyr::select(-contains("Percent"))%>%
  mutate(transect=substr(transect, nchar(transect), nchar(transect)),
         transect=case_when(transect=="0"~paste0("1",transect),
                            transect!=0~transect))

rapid_dat4<- pivot_longer(rapid_dat2,cols=contains("Percent"),names_to = "transect", values_to = "percent")%>%
  dplyr::select(id,transect,percent,year,Date,SHARPTide,Time,data)%>%
  mutate(transect=substr(transect, 3, 4),
         transect=case_when(substr(transect,2,2)=="P"~substr(transect,1,1),
                            substr(transect,2,2)!="P"~transect),
         percent=ifelse(percent=="NULL",NA,as.numeric(as.character(percent))))

# join the species column and percentage columns into one table
rapid_dat5<- left_join(rapid_dat3,rapid_dat4,by=c("id","transect","year","Date","SHARPTide","data","Time"))%>%
  dplyr::select(-transect)%>%
  # we grouped species into bigger categories, so now we need to sum together the individual percentages within each category
  group_by(id,year,Date,SHARPTide,Lat, data, Time,
           LowMarshCC,SaltMarshTBorderCC,HighMarshCC,BrackishTBorderCC,InvasivesCC,
           PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC,DeadSnags,species)%>%
  summarise(percent=sum(percent,na.rm=T))%>%
  ungroup()%>%
  # Then pivot wider so each species category is a variable (column) with percent cover in the plot as the values
  pivot_wider(names_from = species,values_from = percent)%>%
  # remove the NA species column (make sure its column 16)
  dplyr::select(-18)
# fill in NA's with 0's if data was recorded for either the species section or cover class section
is.zero <- function(x) {
  x == 0
}

rapid_dat6<-rapid_dat5%>%
  mutate(across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                  PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),
                ~ifelse(is.na(.x),0,.x)),
         across(c(alt_pct,alt_short_pct,alt_tall_pct,distichlis_pct,gerardii_pct,patens_pct,other_sp_pct,other_marsh_sp_pct,other_aquatic_sp_pct,
                  invasive_pct,sedges_pct,grasses_pct,modified_habitat_pct,shrubs_pct,trees_pct,rushes_pct,upland_pct,
                  phrag_pct,pannes_channels_pct,open_water_pct,bare_pct),
                ~ifelse(is.na(.x),0,.x)),
         # mark if all data is missing for the cover class or species sections, this part of the survey was probably not conducted
         CC_available=ifelse(if_all(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                                      PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),is.zero),0,1),
         sp_pct_available=ifelse(if_all(c(alt_pct,alt_short_pct,alt_tall_pct,distichlis_pct,gerardii_pct,patens_pct,other_sp_pct,other_marsh_sp_pct,other_aquatic_sp_pct,
                                          invasive_pct,sedges_pct,grasses_pct,modified_habitat_pct,shrubs_pct,trees_pct,rushes_pct,upland_pct,
                                          phrag_pct,pannes_channels_pct,open_water_pct,bare_pct),is.zero),0,1)
  )%>%
  filter(!(CC_available==0&sp_pct_available==0))

## Add dominant species and species presence variables:
# if species has more than 0%, mark as present with a 1, otherwise 0
rapid_dat7<-rapid_dat6%>%
  mutate(across(c(alt_pct,alt_short_pct,alt_tall_pct,distichlis_pct,gerardii_pct,patens_pct,other_sp_pct,other_marsh_sp_pct,other_aquatic_sp_pct,
                  invasive_pct,sedges_pct,grasses_pct,modified_habitat_pct,shrubs_pct,trees_pct,rushes_pct,upland_pct,
                  phrag_pct,pannes_channels_pct,open_water_pct,bare_pct),
                ~ifelse(.x>0,1,0)),
         across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                  PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),
                ~ifelse(.x>0,1,0)))

# list the species with the highest percent as the dominant species
rapid_dat7$dom_species<-colnames(
  rapid_dat7[,c("alt_pct","alt_short_pct","alt_tall_pct","distichlis_pct","gerardii_pct","patens_pct","other_sp_pct","other_marsh_sp_pct","other_aquatic_sp_pct",
                "invasive_pct","sedges_pct","grasses_pct","modified_habitat_pct","shrubs_pct","trees_pct","rushes_pct","upland_pct",
                "phrag_pct","pannes_channels_pct","open_water_pct","bare_pct")]
)[
  apply(rapid_dat7[,c("alt_pct","alt_short_pct","alt_tall_pct","distichlis_pct","gerardii_pct","patens_pct","other_sp_pct","other_marsh_sp_pct","other_aquatic_sp_pct",
                      "invasive_pct","sedges_pct","grasses_pct","modified_habitat_pct","shrubs_pct","trees_pct","rushes_pct","upland_pct",
                      "phrag_pct","pannes_channels_pct","open_water_pct","bare_pct")],1,which.max)
]
#format the species names by removing the percent column marker
rapid_dat7$dom_species<-gsub("_pct","",rapid_dat7$dom_species)

# list the cover class with the highest percent as the dominant cover class
rapid_dat7$dom_cover<-colnames(
  rapid_dat7[,c("LowMarshCC","HighMarshCC","SaltMarshTBorderCC","BrackishTBorderCC","InvasivesCC", 
                "PannesChannelsCC","UplandCC","WrackCC","OpenWaterCC")]
)[
  apply(rapid_dat7[,c("LowMarshCC","HighMarshCC","SaltMarshTBorderCC","BrackishTBorderCC","InvasivesCC", 
                      "PannesChannelsCC","UplandCC","WrackCC","OpenWaterCC")],1,which.max)
]
# format the cover class names by removing the cover class column marker
rapid_dat7$dom_cover<-gsub("CC","",rapid_dat7$dom_cover)

# make the dominant categories NA when there is no species or cover class data
rapid_dat7<-rapid_dat7%>%
  mutate(dom_cover=ifelse(CC_available==0,NA,dom_cover),
         dom_species=ifelse(sp_pct_available==0,NA,dom_species))

# adjust the columns markers to presence absence instead of percent
names(rapid_dat7)[c(18:38)] <- gsub("_pct","_pres",names(rapid_dat7[,c(18:38)]))
names(rapid_dat7)[c(8:16)] <- gsub("CC","_pres",names(rapid_dat7[,c(8:16)]))
# and join the percent and presence variables into one table
rapid_dat8<-left_join(rapid_dat6,rapid_dat7,by=c("id","year","Date","SHARPTide","Lat","DeadSnags","data","Time","CC_available","sp_pct_available"))



## Summarize data availability and filter for data completeness:
# number of missing tide records
table(rapid_dat8$SHARPTide)
# number missing time records
table(rapid_dat6$Time)

# number of records with both cover class and species data (Complete data)
table(rapid_dat6$CC_available,rapid_dat6$sp_pct_available)

#look for surveys that may have collected data over multiple days (records with same site and year, missing 1 section of data)
part_missing<-rapid_dat8%>%
  filter(sp_pct_available==0|CC_available==0)
split_survey<-part_missing%>%
  filter(duplicated(id,year))

# majority of records have complete data, remove those missing data
rapid_dat9<-rapid_dat8%>%
  filter(!(sp_pct_available==0|CC_available==0|Time=="0:00:00"))
# How many observations were removed during filtering?
nrow(rapid_dat8)-nrow(rapid_dat9) # filtering for complete data removed 893 records,
nrow(rapid_dat)-nrow(rapid_dat9) # in total removed 1,430 records during cleaning (including duplicate records)

# Look at duplicate data
nup<-rapid_datup%>%
  group_by(id,Date)%>%
  summarise(count=n())%>%
  filter(count>1)
# how many and which are duplicated?
dupsup<-nup%>%left_join(rapid_datup,by=c("id","Date"))%>%
  group_by_all() %>%
  mutate(duplicated = n() > 1)%>%
  filter(duplicated==F)
t<-anti_join(dups,dupsup,by=c("id","Date"))
write.csv(dups,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/rapid_veg_duplicates.csv"),row.names = F)
# which observers have duplicate surveys?
dups2<-dups%>%group_by(ObserverInitials)%>%
  summarise(num_dups=n())
write.csv(dups2,paste0(path_out,"Intermediate_outputs/Data_cleaning_notes/rapid_veg_duplicates_obsInitials.csv"),row.names = F)


# number of survey years per point
nyear<-rapid_dat9%>%
  group_by(id)%>%
  summarize(nyear=length(unique(year)))%>%
  ungroup()%>%
  arrange(desc(nyear))
head(nyear,20)
hist(nyear$nyear)

# points that were surveyed in 2011 and 2022 (either end of the time period)
max(as.numeric(as.character(rapid_dat9$year)));min(as.numeric(as.character(rapid_dat9$year)))
t<-rapid_dat9%>%
  group_by(id)%>%
  filter(("2011"%in%unique(year)|"2012"%in%unique(year))&("2021"%in%unique(year)|"2022"%in%unique(year)))
table(t$SHARPTide)
length(unique(t$id))
length(unique(rapid_dat9$id))

# points that were surveyed in 2014-2018 (middle of the time period)
t2<-rapid_dat9%>%
  filter(year%in%c("2014","2015","2016","2017","2018"))
table(t2$SHARPTide)
length(unique(t2$id))
length(unique(rapid_dat9$id))


## look at vegetation trends through time
t3<-t%>%
  group_by(id)%>%
  filter(year%in%c(max(as.numeric(as.character(year))),min(as.numeric(as.character(year)))))
ggplot(rapid_dat9%>%filter(id=="192567_HHO00"),aes(y=distichlis_pct,x=year,group=SHARPTide))+
  geom_point(aes(color=as.factor(SHARPTide)))+
  geom_smooth(method="lm",aes(linetype=as.factor(SHARPTide)),color="black")+
  #scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  #labs(color = "Nest Site", linetype= "Nest Site",y = "Raw Reflectance (PCA)", x="Proportion High Marsh") + 
  theme_classic(base_size = 12)

ggplot(t%>%filter(year%in%c("2011","2022")&id%in%unique(t$id)[1:10]),aes(group=id,x=as.numeric(as.character(year)),y=as.numeric(as.character(HighMarshCC))))+
  geom_line(aes(color=id))

