library(tidyverse)
library(sf)
#https://rspatial.org/sdm/1_sdm_introduction.html

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. parameters
# what prediction resolution?
reso<-30

# use random background points or veg plots as absences?
ab_type<-"b"
#ab_type<-"v" 

# pick which predictors and response to use for models
predictors <- "uvvr"
#predictors <- "no uvvr" #removes observations mising uvvr data


# veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))

## 3. Load observations with predictors
# list of predictor files for each zone as a raster stack
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))

# point predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_local.csv"))%>%
  filter(bp%in%c("p",ab_type))%>%
  mutate(presence=ifelse(bp=="p",1,0))%>%
  left_join(veg_codes,by="veg_class")

# buffered predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"))%>%
  dplyr::select(id,HIMARSH,LOMARSH,POOL,PHRG,STRM,MUD,TERRBRD)%>% #remove UPLND, keep phrag since directly management related
  right_join(dat,by="id")%>%
  mutate(veg_code=as.factor(veg_code),
         #create binary variables for each veg class
         Highmarsh=as.factor(ifelse(veg_class=="HIMARSH",1,0)))


# filter nests to thinned data
filt<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_10_20_fail_thin_1m_dist_errors_removed.shp"))
dat<-filter(dat,id%in%filt$id | bp%in%c("b","v"))

## 4. Remove missing values from predictors
#data with all vegetation class and NAIP predictors
dat_comp<-dat[complete.cases(dat[ , c("veg_code","cor_txt")]),]
dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]<-dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]%>%
  replace(is.na(.),0)
#associated predictors
form<-y~ndvi+cor_txt+Highmarsh+HIMARSH+ent_txt+pca#+latitude 

#data with all predictors including UVVR
dat_comp_all<-dat_comp[complete.cases(dat_comp[ , c("uvvr_mean")]),]
#associated predictors
form_all<-y~ndvi+cor_txt+Highmarsh+HIMARSH+ent_txt+pca+uvvr_mean#+latitude
#form_all<-y~ndvi+cor_txt+veg_code+HIMARSH+ent_txt+pca+uvvr_mean#+latitude


## 5. pick which predictors to use for models

# don't use UVVR data (increases data availability)
if(predictors=="no uvvr"){
  pres_dat<-dat_comp
  surv_dat<-dat_comp[!is.na(dat_comp$fate),] #filter just the nests with known fates
  mod_form<-form
}

# use UVVR data (limits data availability by more than half)
if(predictors=="uvvr"){
  pres_dat<-dat_comp_all
  surv_dat<-dat_comp_all[!is.na(dat_comp_all$fate),] #filter just the nests with known fates
  mod_form<-form_all
}





## 6. Check for balance among sample groups 

#Remove Years with less than 5 observations
yr_n<-summarise(group_by(pres_dat,Year),count=n())
pres_dat<-filter(pres_dat,is.na(Year)|Year!=2010)

#remove sites with less than 5 observations
site_n<-summarise(group_by(pres_dat,site),count=n())

pres_dat<-filter(pres_dat,(is.na(site)|pres_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|pres_dat$Year%in% yr_n[yr_n$count>=5,]$Year))
surv_dat<-filter(surv_dat,(is.na(site)|surv_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|surv_dat$Year%in% yr_n[yr_n$count>=5,]$Year))

#balance in responses
pres_ratio<-table(pres_dat$presence)
pres_ratio[2]/pres_ratio[1]

surv_ratio<-table(surv_dat$fate)
surv_ratio[2]/surv_ratio[1]



## 7. divide into testing and training data (80 train/20 test)
pres_split<-pres_dat %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))

surv_split<-surv_dat %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))


#check ratio of data
map_int(pres_split, nrow)
map_int(surv_split, nrow)

#check balance in each split dataset
pres_ratio2<-table(pres_split[["train"]]$presence) #retains same balance
pres_ratio2[2]/pres_ratio2[1]
surv_ratio2<-table(surv_split[["train"]]$fate) #retains same balance
surv_ratio2[2]/surv_ratio2[1]

#split into separate objects
pres_train<-pres_split[["train"]]%>%mutate(y=presence)
pres_test<-pres_split[["test"]]%>%mutate(y=presence)
surv_train<-surv_split[["train"]]%>%mutate(y=fate)
surv_test<-surv_split[["test"]]%>%mutate(y=fate)













