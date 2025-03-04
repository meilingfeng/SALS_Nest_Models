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
# rapid vegetation data (measured within a 50m radius circle)
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/Rapid_Vegetation_20230303_updated.csv"),na.strings=c("",".","NULL"))%>%
  rename(id=BirdPtID)%>%
  mutate(Date=mdy(Date),
         year=as.factor(year(Date)),
         ObserverInitials=ifelse(is.na(ObserverInitials),as.factor(CooperatorCode),as.factor(ObserverInitials)))%>%
  #2023 had not been QA/QC-ed yet 
  filter(year!=2023)


## 3. Add location data for the veg points
#--------------------------------------------------------------------
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/All_points_attribute_table_20230426.csv"))%>%
  rename(id=region_ID)%>%
  dplyr::select("id","Long"="POINT_X","Lat"="POINT_Y","PatchID","local_ID")%>%
  mutate(PatchID=ifelse(PatchID=="-9999",sub("[[:digit:]].*","",local_ID),PatchID))%>%
  dplyr::select(-local_ID)%>%
  right_join(rapid_dat,by="id")%>%
  filter(!(is.na(Lat)&is.na(Long)))




## 5. General data tidying
#------------------------------------------------------------
rapid_dat2<-rapid_dat%>%
  ## remove rows that have no data (all NAs) or that are missing tide and observer information
  filter(!(if_all(contains(c("Dom","Sp","Dead")), ~is.na(.x))))%>%
  filter(!is.na(ObserverInitials)|!is.na(SHARPTide))%>%
  ## Select important variables (Dont include categorical cover classes here, need to model them to account for error):
  dplyr::select(id, PatchID,Lat, Long, year,Date,SHARPTide,Time,ObserverInitials,
                #a) number of dead snags (dead trees) >1m (DBH or height?) within 50m radius
                DeadSnags,
                #b) Percent cover of dominant species (don't need to sum to 100%)
                # species names
                DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10,
                # species percent cover
                Sp1Percent,Sp2Percent,Sp3Percent,Sp4Percent,Sp5Percent,Sp6Percent,Sp7Percent,Sp8Percent,Sp9Percent,Sp10Percent)%>%
  ## condense tide into high vs low tide, use to account for observer ability to see the entire marsh.
  mutate(SHARPTide=case_when(
    SHARPTide%in%c("High","High/Rising","High/Falling")~"High",
    SHARPTide%in%c("Low","Low/Rising","Low/Falling")~"Low"),
    SHARPTide=as.factor(SHARPTide))%>%
  ## remove duplicate records and pick the latest survey within a given year
  distinct(pick(c("id","year","Date","Time")),.keep_all = TRUE)%>%
    # order in descending dates. Duplicate only removes the records after the first occurrence.
  group_by(id,year)%>%
  arrange(desc(Date))%>%
  filter(!duplicated(id,year))%>%
  ungroup()


## 5. Format vegetation variables
#--------------------------------------------------

## a) Find common species not listed on the veg survey as major communities/cover classes
#----
# count how common/rare species are in the study region
rapid_sp_sum1<-rapid_dat%>%
  pivot_longer(starts_with("DomSp"),names_to = "DomSp_num",values_to = "Species")
rapid_sp_sum2<-rapid_sp_sum1%>%
  group_by(Species)%>%
  summarise(n_sites=length(unique(PatchID)),
            n_years=length(unique(year)))%>%
  ungroup()
#store species commonality
if(!file.exists(paste0(path_out,"Final_outputs/Model_Results/plant_species_freq",".csv"))){
  write.csv(rapid_sp_sum2,paste0(path_out,"Final_outputs/Model_Results/plant_species_freq",".csv"), row.names = F)
}

# Focus on species that can be used as a reference throughout the breeding range (occur at 50 or more marsh patches)
      # 1/4 number of patches for reference...
length(unique(rapid_dat$PatchID))/4

# get rare species (<50 patches) and common species
rare<-rapid_sp_sum2%>%
  filter(n_sites<50)
common<-rapid_sp_sum2%>%
  filter(n_sites>=50)


# get other common species that aren't listed within a specific community/cover class in the veg survey sheet
common2<-mutate(common,
              #remove species that will get lumped into broader categories
              Species=case_when(
                grepl("lterniflora|hragmites|atens|stichlis|gerardii|Spartina cynosuroides|Bolboschoenus robustus|Scirpus robustus|Juncus sp|Juncus roemerianus|Limonium sp|Limonium carolinianum|Schoenoplectus pungens|Scirpus pungens|Schoenoplectus sp|Symphyotrichum tenuifolium|Symphyotrichum sp|Schoenoplectus americanus|Scirpus americanus|Limonium nashii|Triglochin maritima|Typha|Spartina pectinata|Solidago sempervirens|Baccharis halimifolia|Iva|Conifer|Acer|Catalpa|Betula|C. canadensis|Fagus|Fraxinus|Juglans|Juniperus|Pinus|Platanus|Populus|Prunus|Quercus|Rhus|Rhamnus|Salix|Upland|upland|Mowed Grass|Bare Ground (road/pavement)|Ditch|ditch|Creek|creek|Pool|pool|Panne|panne|Channel|channelAlgae|S. distichum|Fucus edentatus|Non-saltmarsh|Unknown|water|Water|Lythrum salicaria|Lonicera japonica",Species)~NA,
                !grepl("lterniflora|hragmites|atens|stichlis|gerardii|Spartina cynosuroides|Bolboschoenus robustus|Scirpus robustus|Juncus sp|Juncus roemerianus|Limonium sp|Limonium carolinianum|Schoenoplectus pungens|Scirpus pungens|Schoenoplectus sp|Symphyotrichum tenuifolium|Symphyotrichum sp|Schoenoplectus americanus|Scirpus americanus|Limonium nashii|Triglochin maritima|Typha|Spartina pectinata|Solidago sempervirens|Baccharis halimifolia|Iva|Conifer|Acer|Catalpa|Betula|C. canadensis|Fagus|Fraxinus|Juglans|Juniperus|Pinus|Platanus|Populus|Prunus|Carya|Quercus|Rhus|Rhamnus|Ulmus|Salix|Upland|upland|Mowed Grass|Bare Ground (road/pavement)|Ditch|ditch|Creek|creek|Pool|pool|Panne|panne|Channel|channel|Algae|S. distichum|Fucus edentatus|Non-saltmarsh|Unknown|water|Water|Lythrum salicaria|Lonicera japonica",Species)~Species
              ))%>%
  filter(!is.na(Species))%>%
  mutate(genus=sub("[[:space:]].*","",Species))%>%
  dplyr::select(Species,genus)
  
#if a rare species can't be combined into a cover/community category, combine them by genus and see if they are still rare
rare2<-mutate(rare,
          #remove species that will get lumped into broader categories
       Species=case_when(
         grepl("lterniflora|hragmites|atens|stichlis|gerardii|Spartina cynosuroides|Bolboschoenus robustus|Scirpus robustus|Juncus sp|Juncus roemerianus|Limonium sp|Limonium carolinianum|Schoenoplectus pungens|Scirpus pungens|Schoenoplectus sp|Symphyotrichum tenuifolium|Symphyotrichum sp|Schoenoplectus americanus|Scirpus americanus|Limonium nashii|Triglochin maritima|Typha|Spartina pectinata|Solidago sempervirens|Baccharis halimifolia|Iva|Conifer|Acer|Catalpa|Betula|C. canadensis|Fagus|Fraxinus|Juglans|Juniperus|Pinus|Platanus|Populus|Prunus|Quercus|Rhus|Rhamnus|Salix|Upland|upland|Mowed Grass|Bare Ground (road/pavement)|Ditch|ditch|Creek|creek|Pool|pool|Panne|panne|Channel|channelAlgae|S. distichum|Fucus edentatus|Non-saltmarsh|Unknown|water|Water|Lythrum salicaria|Lonicera japonica",Species)~NA,
         !grepl("lterniflora|hragmites|atens|stichlis|gerardii|Spartina cynosuroides|Bolboschoenus robustus|Scirpus robustus|Juncus sp|Juncus roemerianus|Limonium sp|Limonium carolinianum|Schoenoplectus pungens|Scirpus pungens|Schoenoplectus sp|Symphyotrichum tenuifolium|Symphyotrichum sp|Schoenoplectus americanus|Scirpus americanus|Limonium nashii|Triglochin maritima|Typha|Spartina pectinata|Solidago sempervirens|Baccharis halimifolia|Iva|Conifer|Acer|Catalpa|Betula|C. canadensis|Fagus|Fraxinus|Juglans|Juniperus|Pinus|Platanus|Populus|Prunus|Carya|Quercus|Rhus|Rhamnus|Ulmus|Salix|Upland|upland|Mowed Grass|Bare Ground (road/pavement)|Ditch|ditch|Creek|creek|Pool|pool|Panne|panne|Channel|channel|Algae|S. distichum|Fucus edentatus|Non-saltmarsh|Unknown|water|Water|Lythrum salicaria|Lonicera japonica",Species)~Species
       ))%>%
  filter(!is.na(Species))%>%
  mutate(genus=sub("[[:space:]].*","",Species))%>%
  right_join(rapid_sp_sum1,by="Species")%>%
  filter(!is.na(genus))
rare3<-rare2%>%
  group_by(genus)%>%
  summarise(n_sites=length(unique(PatchID)))%>%
  filter(n_sites>=50)
rare4<-filter(rare2,genus%in%rare3$genus)%>%
  distinct(Species,genus)


#Join common species/genus that aren't included in major veg/cover classes in veg data sheet into genus groups
additional_spp<-rbind(common2,rare4)%>%
  filter(Species!="NA")%>%
  group_by(genus)%>%
  summarise(Species=list(unique(Species)))%>%
  #remove uninformative groups
  filter(genus!="Angiosperm")
# 2 additional veg groups...




## b) Aggregate species/cover names
#----
rapid_dat3<-rapid_dat2%>%
  mutate(across(c(DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10),
                #Group species based on broader habitat and plants that SALS select for nesting
                ~case_when(#alt_tall_pct- tall alterniflora (associated with low marsh)
                           grepl("lterniflora-tall",.x)~"alt_tall_pct",
                           #phrag_pct can replace nesting habitat
                           grepl("hragmites",.x)~"phrag_pct",
                           #alt_short_pct- short alterniflora (associated with high marsh)
                           grepl("lterniflora-short",.x)~"alt_short_pct",
                           #alt_pct -unspecified alterniflora (add to both high and low marsh after).
                           grepl("lterniflora",.x)&!grepl("lterniflora-tall|lterniflora-short",.x)~"alt_pct",
                           #patens_pct (associated with high marsh)
                           grepl("atens",.x)~"patens_pct",
                           #distichlis_pct (associated with high marsh)
                           grepl("stichlis",.x,)~"distichlis_pct",
                           #gerardii_pct (associated with high marsh)
                           grepl("gerardii",.x)~"gerardii_pct",
                           #low_marsh_pct- other common low marsh species per veg protocol. Combine with alt_tall_pct + alt_pct after.
                           grepl("Spartina cynosuroides",.x)~"low_marsh_pct",
                           #high_marsh_pct- other common high marsh species per veg protocol. Combine with alt_tall_pct, patens_pct, distichlis_pct, and gerardii_pct after. Bolboschoenus, schoenoplectus, and Scirpus are name changes.
                           grepl("Bolboschoenus robustus|Scirpus robustus|Juncus sp|Juncus roemerianus|Limonium sp|Limonium carolinianum|Schoenoplectus pungens|Scirpus pungens|Schoenoplectus sp|Symphyotrichum tenuifolium|Symphyotrichum sp|Schoenoplectus americanus|Scirpus americanus|Limonium nashii|Triglochin maritima",.x)~"high_marsh_pct", #(pectinata,cynosuroides) exclude alt and patens
                           #Shrubs
                           grepl("Baccharis halimifolia",.x)~"Baccharis_pct",
                           grepl("Iva",.x)~"Iva_pct",
                           #brackish_border_pct- common species
                           grepl("Typha|Spartina pectinata",.x)~"brackish_border_pct", #(angustifolia,latifolia)
                           #saltmarsh border
                           grepl("Solidago sempervirens",.x)~"saltmarsh_border_pct", #(sempervirens,graminifolia)
                           #trees
                           grepl("Conifer|Acer|Catalpa|Betula|C. canadensis|Fagus|Fraxinus|Juglans|Juniperus|Pinus|Platanus|Populus|Prunus|Quercus|Rhus|Rhamnus|Carya|Salix|Ulmus",.x)&.x!="Angiosperm/Conifer shrub" ~"trees_pct",
                           #Natural and developed
                           grepl("Upland|upland|Mowed Grass|Bare Ground (road/pavement)",.x)~"upland_pct",
                           #Pannes, Pools, Creeks(channels, creeeks, ditches, pools/pannes)
                           grepl("Ditch|ditch|Creek|creek|Pool|pool|Panne|panne|Channel|channel",.x)~"pool_panne_channel_pct",
                           #Open water
                           grepl("Water|water",.x)&!grepl("ditch|creek|pool|panne",.x)~"water_pct", #(impoundment, Lemna is duckweed, usually just in standing water)
                           #Invasives (add phrag. also  Lythrum salicaria)
                           grepl("Lythrum salicaria|Lonicera japonica",.x)~"invasives_pct",
                           #other less common species that occur in greater than 20 marsh patches. Grouped by genus.
                           .x%in%unlist(additional_spp[1,]$Species)~paste(additional_spp[1,]$genus,"pct",sep = "_"),
                           .x%in%unlist(additional_spp[2,]$Species)~paste(additional_spp[2,]$genus,"pct",sep = "_"),
                           #remove algae, non-descript names, and rare species 
                           (grepl("Algae|S. distichum|Fucus edentatus|Non-saltmarsh|Unknown",.x)|is.na(.x)|.x==".")~NA,
                           TRUE~NA)))%>%
    # adjust the missing data values for snag count variable
  mutate(DeadSnags=ifelse(DeadSnags%in%c("NULL","-1"),NA,as.numeric(DeadSnags)))
  

## c) Restructure the table to have species names as column names and percentages as the values
#------
      
# turn the columns with species names and percentages and make them into 1 column for names and 1 column for percents
rapid_dat4<-  pivot_longer(rapid_dat3,cols=contains("Dom"),names_to = "count", values_to = "species")%>%
   dplyr::select(-contains("Percent"))%>%
   mutate(count=sub("DomSp","",count))

rapid_dat5<- pivot_longer(rapid_dat3,cols=contains("Percent"),names_to = "count", values_to = "percent")%>%
   dplyr::select(id,count,percent,year,Date,SHARPTide,Time,Lat,Long,PatchID,ObserverInitials)%>%
   mutate(count=sub("Percent","",count),
          count=sub("Sp","",count))

    # join the species column and percentage columns into one table
rapid_dat6<- left_join(rapid_dat4,rapid_dat5,by=c("id","count","year","Date","SHARPTide","PatchID","Time","ObserverInitials","Lat","Long"))%>%
  dplyr::select(-count)%>%
    # we grouped species into bigger categories, so now we need to sum together the individual percentages within each category
  group_by(id,year,Date,SHARPTide,Lat, Long,PatchID,Time,ObserverInitials,DeadSnags,species)%>%
   summarise(percent=sum(percent,na.rm=T))%>%
  ungroup()%>%
    # Then pivot wider so each species category is a variable (column) with percent cover in the plot as the values
  pivot_wider(names_from = species,values_from = percent)%>%
    # remove the NA species column (make sure its column 11)
  dplyr::select(-11)


rapid_dat7<-rapid_dat6%>%
         # mark if no common species/communities/cover classes were recorded in a survey
  mutate(sp_pct_available=ifelse(if_all(ends_with("_pct"),is.na),0,1),
         # turn all NA's into 0 counts
         across(ends_with("_pct"),
                ~ifelse(is.na(.x),0,.x)))%>%
         # remove surveys without common groups
  filter(sp_pct_available!=0)%>%
  dplyr::select(-sp_pct_available)



## 6. Summarize data availability
#--------------------------------------------------------------------------------------
  # number of survey years per point
nyear<-rapid_dat7%>%
  group_by(id)%>%
  summarize(nyear=length(unique(year)))%>%
  ungroup()%>%
  arrange(desc(nyear))
head(nyear,20)
hist(nyear$nyear)

  # points that were surveyed in 2011 and re-surveyed 2022
max(as.numeric(as.character(rapid_dat7$year)));min(as.numeric(as.character(rapid_dat7$year)))
t<-rapid_dat7%>%
  group_by(id)%>%
  filter(("2011"%in%unique(year)|"2012"%in%unique(year))&("2021"%in%unique(year)|"2022"%in%unique(year)))
length(unique(t$id))
length(unique(rapid_dat7$id))
#about half of the sites were resurvey sites

  # points that were surveyed in 2014-2018 (middle of the time period)
t2<-rapid_dat7%>%
  filter(year%in%c("2014","2015","2016","2017","2018"))
length(unique(t2$id))
length(unique(rapid_dat7$id))
#most sites surveyed between this period


## look at vegetation trends through time
t3<-t%>%
  group_by(id)%>%
  filter(year%in%c(max(as.numeric(as.character(year))),min(as.numeric(as.character(year)))))
ggplot(rapid_dat7%>%filter(id=="192567_p1"),aes(y=distichlis_pct,x=year,group=SHARPTide))+
  geom_point(aes(color=as.factor(SHARPTide)))+
  geom_smooth(method="lm",aes(linetype=as.factor(SHARPTide)),color="black")+
  #scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  #labs(color = "Nest Site", linetype= "Nest Site",y = "Raw Reflectance (PCA)", x="Proportion High Marsh") + 
  theme_classic(base_size = 12)

ggplot(t%>%filter(year%in%c("2011","2022")&id%in%unique(t$id)[1:10]),aes(group=id,x=as.numeric(as.character(year)),y=as.numeric(as.character(patens_pct))))+
  geom_line(aes(color=id))


rapid_dat8<-rapid_dat7%>%
  #make the high and low marsh pct cover variables include all high and low marsh species. Also lump invasive category.
  rowwise()%>%
  mutate(high_marsh_pct= sum(high_marsh_pct,distichlis_pct,gerardii_pct,patens_pct,alt_short_pct,na.rm = T),
         low_marsh_pct=sum(low_marsh_pct, alt_tall_pct,na.rm=T),
         invasives_pct=sum(invasives_pct,phrag_pct,na.rm=T),
         saltmarsh_border_pct=sum(Baccharis_pct,Iva_pct,saltmarsh_border_pct,na.rm=T))%>%
  ungroup()

write.csv(rapid_dat8, paste0(path_out,"Intermediate_outputs/Survey_Vegetation/processed_rapid_veg.csv"),row.names = F)


