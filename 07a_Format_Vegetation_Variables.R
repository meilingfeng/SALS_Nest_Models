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
dat_path<-"C:/Users/10788/Desktop/SaltMarsh/Data/"
path_out<-"C:/Users/10788/Desktop/SaltMarsh/Outputs/"
## 2. Read in vegetation data
#------------------------------------
# There is veg data from 3 different surveys:
# a) transect survey points (only for conservation/restoration sites, measured in a transect line from upland to coast, centered on rapid veg points)
transect_dat<-read.csv(paste0(dat_path,"Survey Database/Point_Vegetation_20221209.csv"))%>%
  # add consistent columns names to each dataset (id and year)
    rename(id=Point_ID)%>%
    mutate(year=as.factor(Year))
 #          Time=substr(PointInterceptDate,nchar(PointInterceptDate)-4,nchar(PointInterceptDate)),
 #        Date=substr(PointInterceptDate,1,nchar(PointInterceptDate)-5),
 #        year=year(mdy(Date)),
 #        Month=month(mdy(Date)),
  #       Day=day(mdy(Date)))%>%
 # filter(year!=2023) #this year's data has not been QA/QC yet
# b) rapid vegetation data (measured within a 50m radius circle)
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/Rapid_Vegetation_20230303.csv"))%>%
  rename(id=BirdPtID)%>%
  mutate(year=as.factor(year(ymd(Date))))
 # filter(year!="2023")
# c) demographic random vegetation data (measured within 1m square)
demo_dat<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2020.csv"))%>%
  rename(id=VegPointID)%>%
  mutate(year=as.factor(year(mdy(SurveyDate))))%>%
  mutate(across(c(TallestMP4, TallestMP3, TallestMP2, TallestMP1,TallestCent,contains(c("Avg","Comp","T1","T2","T3","T4","Dist","Percent","EntOrient","Thatch"))),
                ~as.numeric(.x)))

## 3. Add the location data for the transect and rapid veg points
#--------------------------------------------------------------------
transect_dat<-read.csv(paste0(dat_path,"Survey Database/All_points_attribute_table_20230426.csv"))%>%
  rename(id=region_ID)%>%
  dplyr::select("id","Long"="POINT_X","Lat"="POINT_Y")%>%
  right_join(transect_dat,by="id")%>%
  filter(!(is.na(Lat)&is.na(Long)))

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

## 1. Demographic Vegetation Data
#--------------------------------

## Surrounding Vegetation data
demo_veg<-demo_dat%>%
  mutate(#add note of which dataset this is
         data="demo_veg")%>%
  # Format original variables
  # Make all the groups in the categorical variables (species lists) consistent 
  mutate(across(c(TallestMP1Sp, TallestMP2Sp, TallestMP3Sp, TallestMP4Sp, TallestCentSp),
                ~case_when(grepl("stichlis",.x)~"distichlis",
                           grepl("atens",.x)~"patens",
                           grepl("lterniflora",.x)~"alterniflora",
                           grepl("gerardii",.x)~"gerardii",
                           grepl("ragmites",.x)~"phrag",
                           grepl("Iva",.x)~"shrub",
                           grepl("Open water|Bare Ground|Wrack",.x)~"unvegetated",
                           grepl("Not rec|NOT REC",.x)~NA,
                           #group rare species into "Other" category
                           !grepl("stichlis|atens|lterniflora|ragmites|Not rec|NOT REC|gerardii|Open water|Bare Ground|Wrack|Iva",.x)~"Other")))%>%
  #Calculate Final Variables
  rowwise()%>%
  #take average of transect data (percent coverage along 4 transects) from 2011, then use calculate dominant species below
  mutate(rackT=mean(c(RackT1,RackT2,RackT3,RackT4),na.rm=T),
         bareT=mean(c(BareGroundT1,BareGroundT2,BareGroundT3,BareGroundT4),na.rm=T),
         waterT=mean(c(WaterT1,WaterT2,WaterT3,WaterT4),na.rm=T),
         patensT=mean(c(SpartinaPatensT1,SpartinaPatensT2,SpartinaPatensT3,SpartinaPatensT4),na.rm=T),
         altT=mean(c(SpartinaAlternifloraT1,SpartinaAlternifloraT2,SpartinaAlternifloraT3,SpartinaAlternifloraT4),na.rm=T),
         dstchT=mean(c(DistichlisSpicataT1,DistichlisSpicataT2,DistichlisSpicataT3,DistichlisSpicataT4),na.rm=T),
         juncusT=mean(c(JuncusGerardiiT1,JuncusGerardiiT2,JuncusGerardiiT3,JuncusGerardiiT4),na.rm=T),
         atripT=mean(c(AtriplexPatulaT1,AtriplexPatulaT2,AtriplexPatulaT3,AtriplexPatulaT4),na.rm=T),
         limT=mean(c(LimoniumCarolinianumT1,LimoniumCarolinianumT2,LimoniumCarolinianumT3,LimoniumCarolinianumT4),na.rm=T),
         salT=mean(c(SalicorniaEuropaeaT1,SalicorniaEuropaeaT2,SalicorniaEuropaeaT3,SalicorniaEuropaeaT4),na.rm=T),
         pucT=mean(c(PuccinelliaAmericanaT1,PuccinelliaAmericanaT2,PuccinelliaAmericanaT3,PuccinelliaAmericanaT4),na.rm=T),
         triT=mean(c(TriglochinMaritimaT1,TriglochinMaritimaT2,TriglochinMaritimaT3,TriglochinMaritimaT4),na.rm=T),
         patens.dstchT=mean(c(Patens.DistichlisT1,Patens.DistichlisT2,Patens.DistichlisT3,Patens.DistichlisT4),na.rm=T),
         ivaT=mean(c(IvaFrutescensT1,IvaFrutescensT2,IvaFrutescensT3,IvaFrutescensT4),na.rm=T),
         
    # Continuous variables:
    # a) thatch depth averaged across 4 sides and center of each plot quadrant (1x1m)
        thatch_depth=mean(c(ThatchMP1,ThatchMP2,ThatchMP3,ThatchMP4,ThatchCent), na.rm=T),
    # b) average veg height averaged across sides and center
        avg_height=mean(c(AvgMP1,AvgMP2,AvgMP3,AvgMP4,AvgCent),na.rm=T),
    # c) tallest veg height across sides and center
        tallest_height=max(c(TallestMP4, TallestMP3, TallestMP2, TallestMP1,TallestCent), na.rm=T),
    # d) percent cover of S. patens within the plot quadrant
        patens_pc=SpartinaPatensComp,
    # e) percent cover of S. alterniflora
        alt_pc=SpartinaAlternifloraComp,
    # f) percent cover of Distichlis
        dstch_pc=DistichlisSpicataComp,
    # g) percent cover of Juncus
        jncsg_pc=JuncusGerardiiComp,
    # h) percent cover of water
        water_pc=WaterComp,
    # i) percent cover of patens and distichlis mix
        patens_dstch_pc=SpartinaPatensAndDistichlisSpicataComp,
    # j) percent cover of wrack
        rack_pc=RackComp,
    # k) percent cover of bare ground
        bare_pc=sum(c(UnvegetatedComp,BareGroundComp),na.rm=T),
    
    # Categorical variables:
    # l) Dominant tallest veg species across sides and center
        dom_tallest=list(find_mode(c(TallestMP1Sp, TallestMP2Sp, TallestMP3Sp, TallestMP4Sp, TallestCentSp))),
    # m) convert percent vegetation to variables indicating the presence of each species within the plot
    # phragmites, wrack, lysimachia, pectinata, triglochin, limonium only have 1-5 records, lump into an "other" category for rare species
        alt=ifelse(alt_pc>0|altT>0,1,0),
        patens=ifelse(patens_pc>0|patensT>0|patens.dstchT>0|patens_dstch_pc>0,1,0),
        juncus=ifelse(jncsg_pc>0|juncusT>0,1,0),
        distichlis=ifelse(dstch_pc>0|dstchT>0|patens.dstchT>0|patens_dstch_pc>0,1,0),
        wrack=ifelse(rack_pc>0|rackT>0,1,0),
        phrag=ifelse(PhagmitesAustralisComp>0,1,0),
        shrub=ifelse(IvaFrutescensComp>0|ivaT>0,1,0),
        other=ifelse(SalicorniaComp>0|OtherSpComp>0|atripT>0|limT>0|salT>0|pucT>0|triT>0, 1,0),
        unvegetated=ifelse(bare_pc>0|water_pc>0|bareT>0|waterT>0,1,0))%>% 
  # n) Dominant vegetation cover type in the plot
  dplyr::select(id,year,Lat,data,
                # select the vegetation percentage variables
                rack_pc,bare_pc,water_pc,patens_pc,alt_pc,dstch_pc,jncsg_pc,patens_dstch_pc,
                rackT,bareT,waterT,patensT,altT,dstchT,juncusT,patens.dstchT,
                #and other continuous/categorical variables:
                #categorical
                alt,patens,juncus,distichlis,wrack,phrag,other,shrub,unvegetated,dom_tallest,
                #continuous
                thatch_depth,tallest_height,avg_height)
# create a variable that lists the column name with the highest percent cover
demo_veg$dom_cover<-colnames(
  demo_veg[,c("rack_pc","bare_pc","water_pc","patens_pc","alt_pc","dstch_pc","jncsg_pc","patens_dstch_pc",
              "rackT","bareT","waterT","patensT","altT","dstchT","juncusT","patens.dstchT")]
  )[
    apply(demo_veg[,c("rack_pc","bare_pc","water_pc","patens_pc","alt_pc","dstch_pc","jncsg_pc","patens_dstch_pc",
                      "rackT","bareT","waterT","patensT","altT","dstchT","juncusT","patens.dstchT")],1,which.max)
    ]
demo_veg<-demo_veg%>%mutate(dom_cover=case_when(grepl("dstch",dom_cover)~"distichlis",
                                 grepl("patens",dom_cover)~"patens",
                                 grepl("alt",dom_cover)~"alterniflora",
                                 grepl("jncsg|juncus",dom_cover)~"gerardii",
                                 grepl("patens_dstch",dom_cover)~"patens.distichlis",
                                 grepl("rack",dom_cover)~"wrack",
                                 grepl("bare",dom_cover)~"unvegetated",
                                 grepl("water",dom_cover)~"unvegetated",
                                 !grepl("dstch|patens|alt|jncsg|juncus|patens_dstch|rack|bare|water",dom_cover)~dom_cover))
  
demo_veg<-dplyr::select(demo_veg,-rackT,-bareT,-waterT,-patensT,-altT,-dstchT,-juncusT,-patens.dstchT)



## Nest Characteristics data
demo_nest<-demo_dat%>%
  # mark what data it is
  mutate(data="demo_nest")%>%
  # use only nest points, not random veg points
  filter(PointType=="Nest")%>%
  # Select Important Variables
  dplyr::select(id, Site, Lat, 
                #Categorical variables for Nests
                NestLoc, # nest is under thatch (dead vegetation), live vegetation, both, or none (exposed)
                WovenNestCanopy, # Does the nest have a woven canopy over it? (None, Partial, Complete)
                WovenNestCanopyLiveDead, # If it has a canopy, it is made with live or thatch/dead material or both?
                WovenNestCanopyVegTypes, # If it has a canopy, what plant species is it made of?
                VertStructureLiveDead, # Is the vertical structure (stems) the nest is attached to live or dead material or both?
                VertStructureVegTypes, # plant species in vertical structure
                NumVertStructureAttachments, # the the vertical structure made of single stems or multiple stems?
                #Continuous variables for Nests
                EntOrientation,# degree bearing of the nest entrance
                DistLipToCent, #Distance from lip to nest center
                DistLipToGround, #Distance from lip of nest to ground
                DistBottomToGround, # Distance from bottom of nest to ground (maybe subtract these to get nest depth?)
                PercentVisible # Percent canopy cover over the nest
  )%>%
  # Aggregate the categories for categorical variables
  mutate(VrtStLD=as.factor(case_when(VertStructureLiveDead%in%c("both", "Both","BOTH")~ "both", 
                                     VertStructureLiveDead%in%c("LIVE","Live")~"live",
                                     VertStructureLiveDead=="Thatch/Dead"~"dead",
                                     VertStructureLiveDead=="NOT REC"~NA)),
         NestLoc=as.factor(case_when(NestLoc%in%c("Under Thatch","Under thatch")~"under thatch",
                                     NestLoc%in%c("Under Both","Under both")~"under both",
                                     NestLoc%in%c("Under Live Growth","Under live growth")~"under live growth",
                                     NestLoc=="Exposed"~ "exposed",
                                     is.na(NestLoc)|NestLoc=="NOT REC"~NA)),
         NmVrtSA= as.factor(case_when(NumVertStructureAttachments%in%c("Multiple stems", "Multiple Stems")~ "multiple", 
                                      NumVertStructureAttachments%in%c("Single stem", "Single Stem","One stem")~ "single",
                                      is.na(NumVertStructureAttachments)|NumVertStructureAttachments=="NOTREC"~NA)),
         WvnNstC= as.factor(case_when(WovenNestCanopy%in%c("Complete","COMPLETE")~"complete",
                                      WovenNestCanopy%in%c("partial","Partial")~"partial",
                                      WovenNestCanopy=="None"~ "none",
                                      is.na(WovenNestCanopy)|WovenNestCanopy=="NOT REC"~NA)),
         WvnNCLD=as.factor(case_when(WovenNestCanopyLiveDead%in%c("Thatch/dead", "Thatch/Dead")~"dead",
                                     WovenNestCanopyLiveDead=="Both"~"both",
                                     WovenNestCanopyLiveDead=="Live"~ "live",
                                     is.na(WovenNestCanopyLiveDead)|WovenNestCanopyLiveDead=="NOT REC"~NA)),
         
    #convert woven canopy vegetation to binary variables indicating the presence of each species within the canopy
        canopy_alt=ifelse(grepl("alterniflora|ALTERNIFLORA|alternifora",WovenNestCanopyVegTypes),1,0),
         canopy_patens=ifelse(grepl("patens|PATENS",WovenNestCanopyVegTypes),1,0),
         canopy_juncus=ifelse(grepl("Juncus|JUNCUS",WovenNestCanopyVegTypes),1,0),
         canopy_distichlis=ifelse(grepl("Distichlis",WovenNestCanopyVegTypes),1,0),
         canopy_phrag=ifelse(grepl("Phragmites",WovenNestCanopyVegTypes),1,0),
         canopy_wrack=ifelse(grepl("Wrack",WovenNestCanopyVegTypes),1,0),
         canopy_shrub=ifelse(grepl("Iva|Baccharis",WovenNestCanopyVegTypes),1,0),
         #lysimachia, pectinata, triglochin, limonium only have 1-5 records, so lump into an "other" category for rare species *maybe break up into more categories like the rapid/transect veg?
         canopy_other=ifelse(grepl("Salicornia|Triglochin|Limonium|Carex|Schoenoplectus|Solidago|pectinata|Eleocharis|baltica|Plantago|Potentilla|Lysimachia",WovenNestCanopyVegTypes),1,0),
    #convert vertical structure vegetation to binary variables indicating the presence of each species within the vertical structure     
         vert_alt=ifelse(grepl("alterniflora|ALTERNIFLORA|alternifora",VertStructureVegTypes),1,0),
         vert_patens=ifelse(grepl("patens|PATENS",VertStructureVegTypes),1,0),
         vert_juncus=ifelse(grepl("Juncus|JUNCUS",VertStructureVegTypes),1,0),
         vert_distichlis=ifelse(grepl("Distichlis",VertStructureVegTypes),1,0),
         vert_phrag=ifelse(grepl("Phragmites",VertStructureVegTypes),1,0),
         vert_wrack=ifelse(grepl("Wrack",VertStructureVegTypes),1,0),
         vert_shrub=ifelse(grepl("Iva|Baccharis",VertStructureVegTypes),1,0),
         vert_other=ifelse(grepl("Salicornia|Triglochin|Limonium|Carex|Schoenoplectus|Solidago|pectinata|Eleocharis|baltica|Plantago|Potentilla|Lysimachia",VertStructureVegTypes),1,0)
  )%>%
  dplyr::select(-WovenNestCanopy,-WovenNestCanopyLiveDead,-WovenNestCanopyVegTypes,-VertStructureVegTypes,-VertStructureLiveDead,-NumVertStructureAttachments)





## 2. Rapid vegetation data (Survey database)
#--------------------------------------------------
#Understand which and how many species was recorded in the datasheet
rapid_domsp<-cbind(rep(rapid_dat$id,10),unlist(rapid_dat[,c(27,29,31,33,35,37,39,41,43,45)]),rep(rapid_dat$year,10))%>%
  as.data.frame.matrix()
colnames(rapid_domsp)<-c("id","species","year")
count_unique<-function(x){
  specieslist=unique(x$species)[-1]
  countsite=c()
  countyear=c()
  for (i in 1:length(specieslist)){
    sitelist<-rapid_domsp$id[grepl(paste0('^',specieslist[i],'$'),rapid_domsp$species)]
    yearlist<-rapid_domsp$year[grepl(paste0('^',specieslist[i],'$'),rapid_domsp$species)]
    uniquesite<-unique(sitelist)
    countsite<-c(countsite,length(uniquesite))
    uniqueyear<-unique(yearlist)
    countyear<-c(countyear,length(uniqueyear))
  }
  output=data.frame(species=specieslist,countsite=countsite,countyear=countyear)
  return(output)
}
count<-as.data.frame(count_unique(rapid_domsp))
frequency<-as.data.frame(table(unlist(rapid_dat[,c(27,29,31,33,35,37,39,41,43,45)])))
colnames(frequency)<-c("species","freq")
rapid_sum<-merge(count,frequency,by="species")
rapid_sum
if(!file.exists(paste0(path_out,"Final_outputs/Model_Results/plant_species_freq",".csv"))){
  write.csv(rapid_sum,paste0(path_out,"Final_outputs/Model_Results/plant_species_freq",".csv"), row.names = F)
}
rapid_dat2<-rapid_dat%>%
  ## indicate dataset is rapid_veg
  mutate(data="rapid_veg")%>%
  ## remove rows that have no data (all NAs)
  filter(!(if_all(contains(c("Dom","Sp","Dead","Water","Wrack")), ~is.na(.x))))%>%
  ## remove duplicate records
  distinct(pick(c("id","year","Date","Time")),.keep_all = TRUE)%>%
  ## Select Important Variables:
  dplyr::select(id, Lat, Long, year,Date,SHARPTide,Time,data,
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
                           # 5 -> 51-75%, mark as 4
                           # 6-> 76-100%, mark as 5
                           .x=="4"~3,
                           .x=="5"~4,
                           .x=="6"~5,
                           .x=="NULL"|is.na(.x)~NA)))%>%
                            #factor(LowMarshCC, levels=c('0', '<1', '1-25', '26-50','51-75','76-100'))
    # Aggregate categories for species names
  mutate(across(c(DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10),
                #Group species based on marsh habitat
                ~case_when(
                           grepl("lterniflora-tall",.x)~"alt_tall_pct",
                           grepl("lterniflora-short",.x)~"alt_short_pct",
                           grepl("lterniflora",.x)&!grepl("lterniflora-tall|lterniflora-short",.x)~"alt_pct",
                           grepl("atens",.x)~"patens_pct",
                           grepl("Phragmites",.x)~"phrag_pct",
                           grepl("stichlis",.x,)~"distichlis_pct",
                           grepl("gerardii",.x)~"gerardii_pct",
                           grepl("Iva",.x)~"iva_pct",
                           grepl("Spartina cynosuroides",.x)~"low_marsh_pct",
                           grepl("robustus|Salicornia|americanus|Juncus roemerianus|Limonium|pungens|Glaux",.x)~"high_marsh_pct", #(pectinata,cynosuroides) exclude alt and patens
                           grepl("Typha augustifolia|Spartina pectinata",.x)~"brackish_border_pct", #(angustifolia,latifolia)
                           grepl("Solidago sempervirens|Baccharis halimifolia",.x)~"saltmarsh_border_pct", #(sempervirens,graminifolia)
                           grepl("Conifer",.x)&.x!="Angiosperm/Conifer shrub" ~"trees_pct",
                           grepl("Water|pool/panne",.x)~"water_pct", #(impoundment, Lemna is duckweed, usually just in standing water)
                           grepl("Upland",.x)~"upland_pct",
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
  group_by(id,year,Date,SHARPTide,Lat, Long,data, Time,
           LowMarshCC,SaltMarshTBorderCC,HighMarshCC,BrackishTBorderCC,InvasivesCC,
           PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC,DeadSnags,species)%>%
   summarise(percent=sum(percent,na.rm=T))%>%
  ungroup()%>%
    # Then pivot wider so each species category is a variable (column) with percent cover in the plot as the values
  pivot_wider(names_from = species,values_from = percent)%>%
    # remove the NA species column (make sure its column 16)
  dplyr::select(-19)
    # fill in NA's with 0's if data was recorded for either the species section or cover class section
is.zero <- function(x) {
  x == 0
}

rapid_dat6<-rapid_dat5%>%
  mutate(across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                  PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),
                           ~ifelse(is.na(.x),0,.x)),
         across(c(alt_pct,alt_short_pct,alt_tall_pct,patens_pct,phrag_pct,distichlis_pct,gerardii_pct,iva_pct,
                  low_marsh_pct,high_marsh_pct,brackish_border_pct,saltmarsh_border_pct,trees_pct,
                  water_pct,upland_pct),
                ~ifelse(is.na(.x),0,.x)),
      # mark if all data is missing for the cover class or species sections, this part of the survey was probably not conducted
         CC_available=ifelse(if_all(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                                 PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),is.zero),0,1),
         sp_pct_available=ifelse(if_all(c(alt_pct,alt_short_pct,alt_tall_pct,patens_pct,phrag_pct,distichlis_pct,gerardii_pct,iva_pct,
                                          low_marsh_pct,high_marsh_pct,brackish_border_pct,saltmarsh_border_pct,trees_pct,
                                          water_pct,upland_pct),is.zero),0,1)
                )%>%
  filter(!(CC_available==0&sp_pct_available==0))

## Add dominant species and species presence variables:
  # if species has more than 0%, mark as present with a 1, otherwise 0
rapid_dat7<-rapid_dat6%>%
  mutate(across(c(alt_pct,alt_short_pct,alt_tall_pct,patens_pct,phrag_pct,distichlis_pct,gerardii_pct,iva_pct,
                  low_marsh_pct,high_marsh_pct,brackish_border_pct,saltmarsh_border_pct,trees_pct,
                  water_pct,upland_pct),
                ~ifelse(.x>0,1,0)),
         across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, 
                  PannesChannelsCC,UplandCC,WrackCC,OpenWaterCC),
                ~ifelse(.x>0,1,0)))

  # list the species with the highest percent as the dominant species
rapid_dat7$dom_species<-colnames(
  rapid_dat7[,c("alt_pct","alt_short_pct","alt_tall_pct","patens_pct","phrag_pct","distichlis_pct","gerardii_pct","iva_pct",
                "low_marsh_pct","high_marsh_pct","brackish_border_pct","saltmarsh_border_pct","trees_pct",
                "water_pct","upland_pct")]
)[
  apply(rapid_dat7[,c("alt_pct","alt_short_pct","alt_tall_pct","patens_pct","phrag_pct","distichlis_pct","gerardii_pct","iva_pct",
                      "low_marsh_pct","high_marsh_pct","brackish_border_pct","saltmarsh_border_pct","trees_pct",
                      "water_pct","upland_pct")],1,which.max)
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
rapid_dat7
  # adjust the columns markers to presence absence instead of percent
names(rapid_dat7)[c(19:33)] <- gsub("_pct","_pres",names(rapid_dat7[,c(19:33)]))
names(rapid_dat7)[c(9:17)] <- gsub("CC","CC_pres",names(rapid_dat7[,c(9:17)]))
  # and join the percent and presence variables into one table
rapid_dat8<-left_join(rapid_dat6,rapid_dat7,by=c("id","year","Date","SHARPTide","Lat","Long","DeadSnags","data","Time","CC_available","sp_pct_available"))



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
ggplot(rapid_dat9%>%filter(id=="192567_p1"),aes(y=distichlis_pct,x=year,group=SHARPTide))+
  geom_point(aes(color=as.factor(SHARPTide)))+
  geom_smooth(method="lm",aes(linetype=as.factor(SHARPTide)),color="black")+
  #scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  #labs(color = "Nest Site", linetype= "Nest Site",y = "Raw Reflectance (PCA)", x="Proportion High Marsh") + 
  theme_classic(base_size = 12)

ggplot(t%>%filter(year%in%c("2011","2022")&id%in%unique(t$id)[1:10]),aes(group=id,x=as.numeric(as.character(year)),y=as.numeric(as.character(HighMarshCC))))+
  geom_line(aes(color=id))

## 3. Transect vegetation data (Survey database)
#---------------------------------------------------------
transect_domsp<-table(transect_dat$ScientificName)%>%
  sort()%>%
  as.data.frame()
colnames(transect_domsp)<-c("species","frequency")

#ScientificName, 1-10,Year, Point_ID
transect_dat2<-transect_dat%>%
  ## indicate dataset is transect_veg
  mutate(data="transect_veg")%>%
  ## remove rows that have no data (all NAs)
  filter(!(if_all(contains("X"), ~.x==FALSE)))
    # number of rows removed (16,197 removed of 45,734)
  nrow(transect_dat)-nrow(transect_dat2)

transect_dat3<-transect_dat2%>%
  ## remove duplicate records
  distinct(pick(c("id","Day","Month","Year","Time","ScientificName")),.keep_all = TRUE)
  # number of rows removed (12,333 removed)
  nrow(transect_dat2)-nrow(transect_dat3)
  

transect_dat4<-transect_dat3%>%
  ## Select Important Variables:
  dplyr::select(id,Long,Lat,SurveyDate.1,Month,Day,Year,Time,data,
                #Presence of veg species along 10 points on a transect from upland to coast
                ScientificName, contains("X"))%>%
  ## Fix Formatting: 
  # Aggregate categories for species names
  mutate(ScientificName=case_when(grepl("Ailanthus",ScientificName)~"ailanthus_pct",
                                  grepl("lterniflora-tall",ScientificName)~"alt_tall_pct",
                                  grepl("lterniflora-short",ScientificName)~"alt_short_pct",
                                  grepl("lterniflora",ScientificName)&!grepl("lterniflora-tall|lterniflora-short",ScientificName)~"alt_pct",
                                  grepl("atens",ScientificName)~"patens_pct",
                                  grepl("Phragmites",ScientificName)~"phrag_pct",
                                  grepl("stichlis",ScientificName,)~"distichlis_pct",
                                  grepl("gerardii",ScientificName)~"gerardii_pct",
                                  grepl("Iva",ScientificName)~"iva_pct",
                                  grepl("Spartina cynosuroides",ScientificName)~"low_marsh_pct",
                                  grepl("robustus|Salicornia|americanus|Juncus roemerianus|Limonium|pungens|Glaux|Tripleurospermum|Pluchea|Triglochin",ScientificName)~"high_marsh_pct", #(pectinata,cynosuroides) exclude alt and patens
                                  grepl("Typha|Spartina pectinata",ScientificName)~"brackish_border_pct", #(angustifolia,latifolia)
                                  grepl("Panicum virgatum|Atriplex|Solidago|Baccharis",ScientificName)~"saltmarsh_border_pct", #(sempervirens,graminifolia)
                                  grepl("Non-saltmarsh",ScientificName)~"non-saltmarsh_pct", #(juncus roemerianus,arcticus littoralis) , exclude gerardii?
                                  grepl("Conifer|Prunus|Pinus|Quercus|Acer|Rhus|Liquidambar|Robinia",ScientificName)&ScientificName!="Angiosperm/Conifer shrub" ~"trees_pct",
                                  grepl("Water|water|ditch/creek|Channel|pool/panne",ScientificName)~"water_pct", #(impoundment, Lemna is duckweed, usually just in standing water)
                                  grepl("Upland",ScientificName)~"upland_pct",
                                  grepl("Bare Ground",ScientificName)~"bare_pct",#includes organic matter, sand, mud
                                  grepl("Wrack",ScientificName)~"wrack_pct",
                                  grepl("Thatch",ScientificName)~"thatch_pct"))%>%
  group_by(ScientificName,id,Year,Month,Day,Time,data)%>%
  summarise(across(contains("X"),~ifelse(sum(.x)>0,1,0)))%>%
  ungroup()%>%
  # Restructure table to have species names as column names and count as the values
  pivot_longer(cols = contains("X"),names_to = "sub_point",values_to = "presence")%>%
  pivot_wider(names_from = ScientificName, values_from = presence)
  transect_dat4[is.na(transect_dat4)]<-0
  ## Calculate variables:
  # diversity across the transect (beta diversity) (low diversity means habitat is similar, high diversity means habitat transitions from upland to coast)
  transect_dat4$diversity<-NA
  dis<-distinct(transect_dat4,id,Year,Month,Day,Time)
  for(i in 1:nrow(dis)){
    x <- transect_dat4[(transect_dat4$id==dis[i,]$id)&(transect_dat4$Year==dis[i,]$Year)&
                         (transect_dat4$Month==dis[i,]$Month)&(transect_dat4$Day==dis[i,]$Day),-c(1:7)]
    div <- diversity(x, MARGIN = 1, index="shannon") #diversity within sites
    Nreg <- colSums(x)
    gamma <- diversity(Nreg, index="shannon")#diversity across sites
    mean_alpha <- mean(div) 
    transect_dat4[(transect_dat4$id==dis[i,]$id)&(transect_dat4$Year==dis[i,]$Year)&
    (transect_dat4$Month==dis[i,]$Month)&(transect_dat4$Day==dis[i,]$Day), ]$diversity<-(exp(gamma)/exp(mean_alpha)) # hill number sorenson, bigger numbers mean more diverse (max number is total number of sub plots)
  }
  # Dominant species along each transect (species/cover found most frequently across the transect)
transect_dat5<-transect_dat4%>%
  group_by(id,Year,Month,Day,Time,data)%>%
  
    # count number of points along the transect that each species was found
  summarise(across(ailanthus_pct:wrack_pct,sum),
            diversity=mean(diversity,na.rm=T))%>%
  ungroup()
    # list the species with the highest count as the dominant species
transect_dat5$dom_species<-colnames(
  transect_dat5[,c("ailanthus_pct","alt_pct","alt_short_pct","alt_tall_pct","patens_pct","phrag_pct","distichlis_pct","gerardii_pct","iva_pct",
                   "low_marsh_pct","high_marsh_pct","brackish_border_pct","saltmarsh_border_pct","trees_pct",
                   "water_pct","upland_pct","bare_pct","wrack_pct","thatch_pct")]
)[
  apply(transect_dat5[,c("ailanthus_pct","alt_pct","alt_short_pct","alt_tall_pct","patens_pct","phrag_pct","distichlis_pct","gerardii_pct","iva_pct",
                      "low_marsh_pct","high_marsh_pct","brackish_border_pct","saltmarsh_border_pct","trees_pct",
                      "water_pct","upland_pct","bare_pct","wrack_pct","thatch_pct")],1,which.max)
]




## plot the mean vegetation around species nests
t<-demo_veg%>%
  mutate(species=str_sub(id,start=5,end=8))%>%
  group_by(species)%>%
  summarize(patens_pc=mean(patens_pc,na.rm=T),
            alt_pc=mean(alt_pc,na.rm=T),
            dstch_pc=mean(dstch_pc,na.rm=T),
            jncsg_pc=mean(jncsg_pc,na.rm=T),
            water_pc=mean(water_pc,na.rm=T),
            bare_pc=mean(bare_pc,na.rm=T))%>%
  ungroup()%>%
  dplyr::filter(species%in%c("SALS","SESP","NESP","CLRA","WILL"))%>%
  pivot_longer(-1,names_to = "veg.spp",values_to = "pct.cover")
t$species<-factor(t$species, levels=c("CLRA","SESP","SALS","WILL","NESP"))
t$veg.spp<-factor(t$veg.spp, levels=c("water_pc","bare_pc","alt_pc",'patens_pc','dstch_pc','jncsg_pc'))

ggplot(t,aes(x=species,y=pct.cover,fill=veg.spp))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()

