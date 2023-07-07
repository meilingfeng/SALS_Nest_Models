

library(tidyverse)
library(terra)
library(sf)


###########################################################################################
# look for associations between:
# 1. vegetation characteristics and nesting probabilities (nest placement and nest success)
#    Do females select for certain veg characteristics associated with nest success?
# 2. nest characteristics and nest success probability
#    Do the choices females make after selecting a nesting site aid in nest success?
###########################################################################################

# if we don't find evidence for these, are extreme flooding events driving success and masking any historical habitat drivers?
# do we see a trend in the quality of nesting habitat over time? Changes in nesting niches over time? 
  # is a change in niche due to adaptation, still within the range limits of the species, or a result of available habitat?
  # How much salt marsh/high marsh have we lost over this period? Are sparrows losing suitable habitat?



### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


# What habitat variables might be important for choosing a nest location?
  # Run an analysis with all data points, looking at veg
dat_all<-read.csv(paste0(path_out,"Final_outputs/Nest_Predictions_at_Veg_Points.csv"), stringsAsFactors = T)%>%
  mutate(year=as.factor(Year))%>%
  dplyr::select(id, Site, 
                #Categorical Co-variates for All Veg data
                NestLoc, # nest is under thatch (dead vegetation), live vegetation, both, or none (exposed)
                WvnNstC, # Does the nest have a woven canopy over it? (None, Partial, Complete)
                WvnNCLD, # If it has a canopy, it is made with live or thatch/dead material or both?
                Wvn
  )

# What variables might be important for nest survival after the female chooses a nest location?
  # Run an analysis with just nest points, looking at nest construction

dat_nest<-read.csv(paste0(path_out,"Final_outputs/Nest_Predictions_at_Veg_Points.csv"), stringsAsFactors = T)%>%
  mutate(year=as.factor(Year))%>%
  filter(PontTyp=="Nest")%>%
  complete.cases()%>%
  dplyr::select(id, Site, 
                #Categorical Covariates for Nests (What is VrtStVT,VrtStAT?)
                NestLoc, # nest is under thatch (dead vegetation), live vegetation, both, or none (exposed)
                WvnNstC, # Does the nest have a woven canopy over it? (None, Partial, Complete)
                WvnNCLD, # If it has a canopy, it is made with live or thatch/dead material or both?
                WvnNCVT, # If it has a canopy, what plant species is it made of?
                VrtStLD, # Is the vertical structure (stems) the nest is attached to live or dead material or both?
                NmVrtSA, # the the vertical structure made of single stems or multiple stems?
                #Continuous Covariates for Nests (What is DstLpTC?)
                Lat,
                DstLpTG, #Distance from lip of nest to ground
                DstBtTG, # Distance from bottom of nest to ground (maybe subtract these to get nest depth?)
                PrcntVs # Percent canopy cover over the nest
                )%>%
  # aggregate the categories for categorical variables
  mutate(NmVrtSA= ifelse(NmVrtSA%in%c("Multiple stems", "Multiple Stems"), "multiple", "single"),
         WvnNstC= case_when(WvnNstC%in%c("Complete","COMPLETE")~"Complete",
                            WvnNstC%in%c("partial","Partial")~"Partial",
                            WvnNstC=="None"~ "None"),
         WvnNCLD=ifelse(WvnNCLD=="Thatch/dead", "Thatch/Dead", WvnNCLD),
         WvnNCVT=case_when(WvnNCVT%in%c("Spartina patens","SPARTINA PATENS")~"patens",
                           WvnNCVT%in%c("Spartina alterniflora","SPARTINA ALTERNIFLORA")~"alterniflora",
                           WvnNCVT%in%c("Spartina patens;Distichlis spicata","Distichlis spicata;Spartina patens")~"patens;distichlis",
                           WvnNCVT%in%c("Spartina alterniflora;Spartina patens","Spartina patens;Spartina alterniflora","SPARTINA ALTERNIFLORA;SPARTINA PATENS")~"patens;alterniflora",
                           WvnNCVT%in%c("Spartina alterniflora;Juncus gerardii","Juncus gerardii;Spartina alterniflora","Juncus geradii;Spartina alterniflora")~"juncus;alterniflora",
                           WvnNCVT%in%c("Spartina alterniflora;Distichlis spicata","SPARTINA PATENS")~"patens",)
    
  )             

unique(dat$WvnNCVT)

### what types of veg variables do we have

### associations between continuous and categorical variables - 
names(which(sapply(dat, class) == "factor"))

### associations between two continuous variables -
names(which(sapply(dat, class) != "factor"))

# correlation

# scatter plots