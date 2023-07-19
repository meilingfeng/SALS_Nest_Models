library(tidyverse)
library(terra)
library(sf)


########################################
# Format Vegetation Characteristics
########################################


### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. Read in vegetation data
# There is veg data from 3 different surveys:
# a) transect survey points (only for conservation/restoration sites, measured in a transect line from upland to coast, centered on rapid veg points)
transect_dat<-read.csv(paste0(dat_path,"Survey Database/Point_Intercept_Vegetation_TRUEFALSE_20221209.csv"), stringsAsFactors = T)%>%
  # add consistent columns names to each dataset (id and year)
  rename(id=Point_ID)%>%
  mutate(year=as.factor(Year))
# b) rapid vegetation data (measured within a 50m radius circle)
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/Rapid_Vegetation_20230303.csv"), stringsAsFactors = T)%>%
  rename(id=BirdPtID)%>%
  mutate(year=as.factor(year(mdy(Date))))
# c) demographic random vegetation data (measured within 1m square)
demo_dat<-read.csv(paste0(dat_path,"Demographic Database/Veg_2011-2020.csv"))%>%
  rename(id=VegPointID)%>%
  mutate(year=as.factor(year(mdy(SurveyDate))))%>%
  mutate(across(c(TallestMP4, TallestMP3, TallestMP2, TallestMP1,TallestCent,contains(c("Avg","Comp","T1","T2","T3","T4","Dist","Percent","EntOrient","Thatch"))),
                ~as.numeric(.x)))

## 3. Add the location data for the transect and rapid veg points
transect_dat<-read.csv(paste0(dat_path,"Survey Database/All_points_attribute_table_20230426.csv"), stringsAsFactors = T)%>%
  rename(id=region_ID)%>%
  dplyr::select("id","Long"="POINT_X","Lat"="POINT_Y")%>%
  right_join(transect_dat,by="id")
rapid_dat<-read.csv(paste0(dat_path,"Survey Database/All_points_attribute_table_20230426.csv"), stringsAsFactors = T)%>%
  rename(id=region_ID)%>%
  dplyr::select("id","Long"="POINT_X","Lat"="POINT_Y")%>%
  right_join(rapid_dat,by="id")



## 3. define function to calculate mode
find_mode <- function(x) {
  u <- unique(na.omit(x))
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}





### Format variables for each dataset
#-----------------------------------------------
## a) Demographic Vegetation Data

# Surrounding Vegetation data
demo_veg<-demo_dat%>%
  mutate(#add note of what dataset this is (good to specify this if we combine across datasets due to different survey methods)
         data="demo_veg")%>%
  ## 1. Format original variables
  # Make all the groups in the categorical variables (species lists) consistent (lump rare species into "Other" category)
  mutate(across(c(TallestMP1Sp, TallestMP2Sp, TallestMP3Sp, TallestMP4Sp, TallestCentSp),
                ~case_when(grepl("stichlis",.x)~"distichlis",
                           grepl("atens",.x)~"patens",
                           grepl("lterniflora",.x)~"alterniflora",
                           grepl("gerardii",.x)~"gerardii",
                           grepl("ragmites",.x)~"phrag",
                           grepl("Iva",.x)~"shrub",
                           grepl("Open water|Bare Ground|Wrack",.x)~"unvegetated",
                           grepl("Not rec|NOT REC",.x)~NA,
                           !grepl("stichlis|atens|lterniflora|ragmites|Not rec|NOT REC|gerardii|Open water|Bare Ground|Wrack|Iva",.x)~"Other")))%>%
  ## 2. Calculate Final Variables
  #take average of transect data (percent coverage along 4 transects) from 2011, then calculate dominant species
  rowwise()%>%
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



# Nest Characteristics data
#--------------------------------------------------------
# Do nest construction decisions influence nest survival? Do more parental decisions than just choosing a good location contribute to nest success?
# Run an analysis with just nest points, looking at nest construction variables

demo_nest<-demo_dat%>%
  mutate(data="demo_nest")%>%
  filter(PointType=="Nest")%>%
  # 1. Select Important Variables
  dplyr::select(id, Site, Lat, 
                #Categorical Covariates for Nests
                NestLoc, # nest is under thatch (dead vegetation), live vegetation, both, or none (exposed)
                WovenNestCanopy, # Does the nest have a woven canopy over it? (None, Partial, Complete)
                WovenNestCanopyLiveDead, # If it has a canopy, it is made with live or thatch/dead material or both?
                WovenNestCanopyVegTypes, # If it has a canopy, what plant species is it made of?
                VertStructureLiveDead, # Is the vertical structure (stems) the nest is attached to live or dead material or both?
                VertStructureVegTypes, # plant species in vertical structure
                NumVertStructureAttachments, # the the vertical structure made of single stems or multiple stems?
                #Continuous Covariates for Nests
                EntOrientation,# degree bearing of the nest entrance
                DistLipToCent, #Distance from lip to nest center
                DistLipToGround, #Distance from lip of nest to ground
                DistBottomToGround, # Distance from bottom of nest to ground (maybe subtract these to get nest depth?)
                PercentVisible # Percent canopy cover over the nest
  )%>%
  #filter(complete.cases(.))%>%
  
  # 2. Aggregate the categories for categorical variables
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
         
    # 3. convert woven canopy vegetation to variables indicating the presence of each species within the canopy
         # lysimachia, pectinata, triglochin, limonium only have 1-5 records, so lump into an "other" category for rare species
         canopy_alt=ifelse(grepl("alterniflora|ALTERNIFLORA|alternifora",WovenNestCanopyVegTypes),1,0),
         canopy_patens=ifelse(grepl("patens|PATENS",WovenNestCanopyVegTypes),1,0),
         canopy_juncus=ifelse(grepl("Juncus|JUNCUS",WovenNestCanopyVegTypes),1,0),
         canopy_distichlis=ifelse(grepl("Distichlis",WovenNestCanopyVegTypes),1,0),
         canopy_phrag=ifelse(grepl("Phragmites",WovenNestCanopyVegTypes),1,0),
         canopy_wrack=ifelse(grepl("Wrack",WovenNestCanopyVegTypes),1,0),
         canopy_shrub=ifelse(grepl("Iva|Baccharis",WovenNestCanopyVegTypes),1,0),
         canopy_other=ifelse(grepl("Salicornia|Triglochin|Limonium|Carex|Schoenoplectus|Solidago|pectinata|Eleocharis|baltica|Plantago|Potentilla|Lysimachia",WovenNestCanopyVegTypes),1,0),
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


#Rapid veg
rapid_dat2<-rapid_dat%>%
  #remove duplicate records
  distinct(.keep_all = TRUE)%>%
  #remove records that are nearly duplicated with less than an hour time difference of the survey start time
  filter((id!="234976_EXC08")&
           (id!="234976_p23"&LowMarshCC!="NULL")&
           (id!="2238900_DCP001"&Time!="11:04:00")&
           (id!="238900_DCP007")&Time!="08:45:00")%>%
  # 1. Select Important Variables
  dplyr::select(id, Lat, year,Date,SHARPTide,
                #Categorical variables:
                #Percent cover within 50 m categories for low and high marsh, saltmarsh border, brackish terrestrial border, invasive species, etc
                LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, PannesChannelsCC,UplandCC,Wrack,OpenWater,
                # number of dead snags (dead trees) >1m (DBH or height?) within 50m radius
                DeadSnags,
                #Continuous variables:
                #Percent cover of individual species
                DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10,
                Sp1Percent,Sp2Percent,Sp3Percent,Sp4Percent,Sp5Percent,Sp6Percent,Sp7Percent,Sp8Percent,Sp9Percent,Sp10Percent)%>%
  # 2. Create variables: 
  #reclassify categorical percent cover variables into 25% intervals
  mutate(across(c(LowMarshCC,HighMarshCC,SaltMarshTBorderCC,BrackishTBorderCC,InvasivesCC, PannesChannelsCC,UplandCC,Wrack,OpenWater),
                ~case_when(  # 0 is not present
                           .x=="0"~as.factor("0"),
                           # 0.5 is less than a percent
                           .x=="0.5"~as.factor("<1"),
                           # groups together 1,2,3 -> 1-25%
                           .x%in%c("1","2","3")~as.factor("1-25"),
                           # 4 -> 26-50%
                           # 5 -> 50-75%
                           # 6-> 76-100%
                           .x=="4"~as.factor("26-50"),
                           .x=="5"~as.factor("51-75"),
                           .x=="6"~as.factor("76-100"),
                           .x=="NULL"|is.na(.x)~NA)))
#factor(LowMarshCC, levels=c('0', '<1', '1-25', '26-50','51-75','76-100'))
rapid_dat3<-rapid_dat2%>%
  mutate(across(c(DomSp1,DomSp2,DomSp3,DomSp4,DomSp5,DomSp6,DomSp7,DomSp8,DomSp9,DomSp10),
                ~case_when(grepl("stichlis",.x,)~"distichlis",
                           grepl("atens",.x)~"patens",
                           grepl("lterniflora-tall",.x)~"alterniflora_tall",
                           grepl("lterniflora-short",.x)~"alterniflora_short",
                           grepl("ragmites",.x)~"phrag",
                           grepl("gerardii",.x)~"gerardii",
                           grepl("lterniflora",.x)&!grepl("lterniflora-tall|lterniflora-short",.x)~"alterniflora",
                           grepl("Conifer|Juniperus|Prunus|Pinus|Quercus|Acer",.x)&.x!="Angiosperm/Conifer shrub" ~"trees",
                           grepl("shrub|Bocconia|Baccharis|Iva|Morella|Myrica|Elaeagnus umbellata",.x)~"shrubs", #Elaeagnus umbellata is invasive
                           grepl("Bare Ground",.x)&.x!="Bare Ground (road/pavement)"~"bare",#includes organic matter, sand, mud
                           grepl("Water|water",.x)&!grepl("(ditch/creek)|(pool/panne)",.x)~"open_water", #(impoundment, Lemna is duckweed, usually just in standing water)
                           grepl("(ditch/creek)|Channel|(pool/panne)",.x)~"pannes/channels",
                           grepl("(road/pavement)|Mowed Grass|(pavement)|(impoundment)",.x)~"modified_habitat", #also ditches and channels?
                           grepl("Upland|Mowed Grass",.x)&.x!="Upland (pavement)"~"non_saltmarsh_habitat", #dont include pavement, include mowed grass?
                           grepl("Non-saltmarsh herbaceous spp.|Pycnanthemum|Lythrum|Alisma|latifolia|Sagittaria|Hibiscus|Helenium amarum|Peltandra|Lespedeza|Symphyotrichum|Erechtites|Ludwigia|Verbena|Glaux|Phytolacca|Pluchea|Angiosperm|Morella|Argentina|Tripleurospermum|Suaeda|Apocynum cannabinum|Ptilimnium|Triglochin|Atriplex|Toxicodendron|Rubus|Unknown|Nuphar|Impatiens|Amaranthus|Pontederia|Sesuvium|Mikania|Celastrus orbiculatus|Lythrum salicaria",.x)&.x!="Angiosperm/Conifer shrub"~"other_species",#divide into terrestrial and aquatic species? Celastrus orbiculatus and Lythrum salicaria is invasive
                           grepl("Juncus",.x)&!grepl("gerardii",.x)~"non_gerardii_rushes", #(juncus roemerianus,arcticus littoralis) , exclude gerardii?
                           grepl("Sedge|Bolboschoenus|Schoenoplectus|Carex|Scirpus|CyperusEleocharis",.x)~"all sedges", #(Schoenoplectus pungens,tabernaemontani)
                           grepl("Poaceae sp.|Leersia|Spartina||Poa sp.|Panicum|Elymus|Festuca|Puccinellia|Ammophila|Thinopyrum|Zizania|Setaria|Echinochloa|Leptochloa",.x)&!grepl("atens|lterniflora",.x)~"non_alt_patens_grasses", #exclude alterniflora and patens? Thinopyrum pycnanthum is invasive
                           grepl("Spartina",.x)~"cord_grasses", #(pectinata,cynosuroides,patens,alterniflora)
                           grepl("Solidago",.x)~"goldenrods", #(sempervirens)
                           grepl("Plantago",.x)~"plantains", #(maritima)
                           grepl("Typha",.x)~"cat_tails", #(angustifolia,latifolia)
                           grepl("Limonium",.x)~"limonium", #(carolinianum,nashii)
                           grepl("Polygonum",.x)~"polygonum", #buckwheat and knotgrasses (pensylvanicum,punctatum,perfoliatum)
                           grepl("Salicornia",.x)~"salicornia", #(depressa,bigelovii)
                           grepl("Wrack",.x)~"wrack",
                           grepl("Thatch",.x)~"thatch",
                           grepl("Lemna|Algae|S. distichum",.x)~NA)))

 rapid_dat4<- pivot_longer(rapid_dat3,cols=contains("Dom"),names_to = "transect", values_to = "species")%>%
   dplyr::select(-contains("Percent"))%>%
   mutate(transect=substr(transect, nchar(transect), nchar(transect)),
          transect=case_when(transect=="0"~paste0("1",transect),
                             transect!=0~transect))
 rapid_dat5<- pivot_longer(rapid_dat3,cols=contains("Percent"),names_to = "transect", values_to = "percent")%>%
   dplyr::select(id,transect,percent,year,Date,SHARPTide)%>%
   mutate(transect=substr(transect, 3, 4),
          transect=case_when(substr(transect,2,2)=="P"~substr(transect,1,1),
                             substr(transect,2,2)!="P"~transect))%>%
   right_join(rapid_dat4,by=c("id","transect","year","Date","SHARPTide"))
#Sci name and percent of any species covering >5% of 50m buffer: DomSp1-10,Sp1-10Percent
unique(rapid_dat$DomSp1)

#Transect
#ScientificName, 1-10,Year, Point_ID
#Remove rows with FALSE for 1-10, remove duplicate rows