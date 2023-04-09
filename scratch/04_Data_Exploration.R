
library(tidyverse)
library(viridis)
library(tmap)
library(sf)
library(spdep) #Used for determining spatial autocorrelation https://rspatial.org/raster/analysis/3-spauto.html


#Project Context
# -------------------------------------------
#Salt marsh habitat is rapidly shrinking. Between rising sea levels and human development preventing migration inland, these important habitats are being squeezed between two distubances.
#Coastal tidal marsh covers about 9% of the earth, yet supports >25% of the human population in addition to its economic importance for the remaining population
#(farmland productivity, fisheries, travel and trade, carbon sequestration, pollutant filtering, and storm buffering).

#On the east coast of North America, these marshes are also home to an endemic species. the Salt Marsh Sparrow, which has been predicted to go extinct by mid century.
#Studies on SALS populations and distributions have shown that high marsh habitat, characteristic of less frequent flooding and unique vegetative communities,
#is associated with SALS occurrence and abundance. However, nests are only attended by females and males do not form territories. This 
#behavior, backed with initial findings, has led us to believe nest locations may occupy different habitat conditions than SALS occurrence.

#Conservation managers would like to know where locations with good habitat for nest placement and nesting success exist. By identifying these reference areas,
#we want to be able to describe manageable habitat conditions that promote successful reproduction in SALS. Previous studies have shown that fine-resolution and local scale
#habitat conditions impact nest occurrence, such as vegetation structure (ie, stem density).

#We are approaching this question using fine-resolution remote sensing data because this data can both capture these local habitat conditions and cover a broad region
#such as the Atlantic Coastline.


## Set file path to data
# -------------------------------------------
dat_path<-"D:/Research/SHARP/Data/"


## Load in Data
# -------------------------------------------
#This dataset contains nests for Salt Marsh Sparrows detected across the Northeast and Mid-Atlantic coastline.
#These nest locations were collected by a network of research collaborators and graduate student projects. Multiple independent studies were compiled across years and study sites.
#All nests were searched for within a stratified random sampling of marsh sites.However, a standardized survey approach was not established until 2010.
#Additionally, since these data were collected by independent projects, the temporal and spatial availability of the data are not consistent.


#Load nest fates with environmental covariates
#fates<-read.csv(paste0(dat_path,"nest_fate_model_dataset.csv"))
fates<-st_read(paste0(dat_path,"Nest_Locations/nest_locations_01_3_23.shp"))%>%
  # filter records to just SALS or species of interest
  filter(Species=="SALS"&
           # also filter records missing coordinate information or that have coordinate errors
           crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         #create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                        fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(name=id,latitude=Lat,longitude=Long,fate,Year,site=site_cd,Region)%>%
  st_drop_geometry()%>%
  filter(!is.na(fate))
  #convert missing data to 0, these are land cover classes that were not found within nest boundaries.
fates[is.na(fates)]<-0


#Located nests were monitored several days to determine the ultimate fate of the fledglings. Fledglings that survived and left the nest
#are classified as successful (fate =1 ) and those that did not (due to flooding, depredation, etc.) were given a fate value of 0. 
#the id field indicates the unique id for each nest, attributed to a specific year, marsh site, and species. For this analysis, I have already filtered the nest
#data to just include SALS nests.

#Let's look at the data structure and availability

nrow(fates) 
#There are 2,795 observations

length(unique(fates$site))
length(unique(fates$region)) 
#There are 51 unique marsh sites across 4 geographic regions (1-Maine/NH, 3=CT, 4 = NY, 5=NJ/VA)

#How are observations distributed across marsh sites?
n_nest_per_site<-summarise(group_by(fates,site),n_nests=n())%>%
  left_join(unique(fates[,c("site","Region")]),by="site")
hist(n_nest_per_site$n_nests) 
  #most sites have 50 or less nest observations, but a couple have more than 300 observations.
ggplot(n_nest_per_site, aes(x=as.factor(region), y=n_nests))+
  geom_boxplot()+
  labs(y= "Number of nest observations", x="Geographic Region")
  #different geographic regions have different representation within the data (on average a lot more nest observations in region 5.) 
  #Region 3 has fewer observations per site on average, with only a couple sites that have huge sampling representation.)

#How many sites per region?



#Majority of sites only have 1 year of sampling
#Site AT, ER, HM, ID, MN, OC all have 6 years of sampling
#Sites EL and CL have the most years of sampling (8 and 7 years)
n_yr_per_site<-summarise(group_by(fates,site),n_years=length(unique(Year)))
n_yr_per_site
ggplot(n_yr_per_site)+
  geom_bar(aes(x=n_years))+
  labs(y="Number of Sites",x="Number of Years Sampled")+
  theme_bw()


#Which years sampled the most, for which sites?
#seems like most sites with longer monitoring efforts began in 2011 and continue through 2020
#Most sites sampled only once occured before 2010, but a few one-off sites were also sampled 2011-2020
ggplot(left_join(distinct(fates[,c("Year","site")], .keep_all = T),n_yr_per_site,by="site"))+
  geom_bar(aes(x=Year,fill=as.factor(n_years)))+
  scale_fill_viridis_d()+
  labs(fill="N Years\nSite Sampled", x="Year","Number of Sites Sampled")+
  theme_bw()

n_nest_per_siteyr<-summarise(group_by(fates,site,yr),n_nests=n())
n_nest_per_siteyr

#Looks like temporal samples from 2006, 2009, 2012, and 2015, 2018 might be good.


#There are two predictor variables of nest success. Each has been resampled within a 15 meter buffer of each nest location.
#The first is the Unvegetated-Vegetated Ratio (UVVR). This is the ratio of vegetation to water or bare land. It has been used as a measure of marsh resilience.
#UVVR was summarized as the average ratio value within each nest buffer area.

#The second is a marsh vegetation community layer. This data layer indicates whether the vegetation community is one of 8
#dominant vegetation/land cover classes (High Marsh, Low Marsh, Pool, Phragmites, Stream, Mud, Upland, or Terrestrial Border). This variable has been divided into individual
#variables for each class, represented as the proportion of each land cover class within the nest buffer area.

#distributions of covariates
hist(fates$uvvr_mean) #right skewed
hist(fates$HIMARSH)#left skewed
hist(fates$LOMARSH) #right skewed
hist(fates$POOL)#right skewed
hist(fates$STRM)#right skewed
hist(fates$UPLND)#right skewed? 1 inflated **what types of distributions would these be?
hist(fates$TERRBRD)#bimodal? **what types of distribution?
hist(fates$PHRG)#**what types of distribution?




# Given the nest sampling bias (independent projects with inconsistent monitoring), I would like to know how to assess the data for spatial and temporal autocorrelation?
#What tests are used to detect these data characteristics?
#What implications do autocorelations have for model selection and interpretation?
#How do can I remove autocorrelation from my data if I find it exists?

#Additionally, I would like to explore ideas for what type of model I should use based on the characteristics of my data and desired outcomes (a perdition surface of nesting success across Atlantic Coast tidal marshes).


#Autocorrleation is a measure of similarity between nearby observations
# why are we concerned about it? The assumption that all observations are independent and random?
# Are we concerned about it in the response or the predictors?

#Temporal autocorrelation - measure the degree of association between consecutive observations over time
#Computing auto-correlation over time


library(lme4)

mod<-glm(fate~uvvr_mean+HIMARSH+LOMARSH+POOL+PHRG+MUD+UPLND+TERRBRD,data = fates,family = "binomial")
summary(mod)
