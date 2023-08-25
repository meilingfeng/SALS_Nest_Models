library(tidyverse)
library(sf)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(dismo)
library(gbm)
source("C:/Users/mefen/OneDrive/Documents/Github/SHARP/Functions/gridSample_sf.R")
#https://rspatial.org/sdm/1_sdm_introduction.html

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## 2. parameters
# what prediction surface resolution?
reso<-30

# use random background points or veg plots as absences?
#ab_type<-"v"
ab_type<-"b" 

# veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))


## 3. Load observations with predictors
# list of predictor surface files for each zone
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","cor_txt", "tideres", "uvvr_diff") 

# point predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_local.csv"))%>%
  filter(bp%in%c("p",ab_type))%>%
  mutate(presence=ifelse(bp=="p",1,0))%>%
  left_join(veg_codes,by="veg_class")%>%
  dplyr::select(-pca,-ndvi,-cor_txt,-ent_txt)

# buffered predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"))%>%
  dplyr::select(id,HIMARSH,ndvi,pca,cor_txt)%>% #remove UPLND
  right_join(dat,by="id")%>%
  mutate(veg_code=as.factor(veg_code),
         #create binary variable for if nests intersect High Marsh habitat
         Highmarsh=as.factor(ifelse(veg_class=="HIMARSH",1,0)))

# fill in 0's for land cover proportions with NAs 
dat[,c("HIMARSH")]<-dat[,c("HIMARSH")]%>%
  replace(is.na(.),0)


# Remove observation records with missing predictor values (regression modeling methods cannot have missing values)
dat_comp<-dat[complete.cases(dat[,all_terms]),]
dat[!(dat$id%in%dat_comp$id),]



### Check Class Balance
#---------------------------------------------------------------------------------------------
# balancing can help machine learning models 


## 1. divide observations into nest presence and nest suvival datasets
pres_dat<-dat_comp
surv_dat<-dat_comp[!is.na(dat_comp$fate),] #filter just the nests with known fates



## 2. look at the prevalence of:

# nest fledges to nest fails
addmargins(table(surv_dat$fate));sum(surv_dat$fate)/nrow(surv_dat)

# prevalence of presence to background
addmargins(table(pres_dat$presence)); sum(pres_dat$presence)/nrow(pres_dat)




# plot proportion of fledges (nest success) across sites
success_site<-dat%>%group_by(site)%>%
  filter(!is.na(fate))%>%
  summarize(pr_success=sum(fate,na.rm=T)/n(),
            lat=mean(latitude,na.rm=T))

hist(success_site$pr_success, breaks=seq(0,1,0.1))
plot(success_site$pr_success~success_site$lat,xlab="Site latitude",ylab="Proportion of nests fledging")
#nest success appears to be more variable at higher latitudes
mean(success_site$pr_success)
#nest success across sites is also 32% on average



## 3. Check sampling bias

# look at spatial distribution and temporal distribution
# does the data need thinning?
# best to thin the more prevalent class (nest fails) to improve calibration - thin to balance

# How are observations distributed across marsh sites?

n_nest_per_site<-summarise(group_by(dat,site),n_nests=n())%>%filter(!(is.na(site)))
hist(n_nest_per_site$n_nests) 
length(n_nest_per_site[n_nest_per_site$n_nests<251,]$n_nests)/nrow(n_nest_per_site)
#96% of sites have 250 or less nest observations, but a couple have more than 300 observations.

#create lat bins
n_nest_per_site<-dat%>%
  mutate(lat_bin = cut(latitude, breaks=8))%>%
  dplyr::select(site,lat_bin)%>%
  distinct()%>%
  right_join(n_nest_per_site,by="site")
ggplot(n_nest_per_site, aes(x=lat_bin, y=n_nests))+
  geom_boxplot()+
  labs(y= "Number of nest observations", x="Geographic Region")+
  geom_jitter(width = 0.2)
# some difference in sampling intensity throughout region, but nothing seems overly under sampled
# lat 42 has the least sites



#16 sites have at least 5 years of sampling
#Site AT, ER, HM, ID, MN, OC all have 6 years of sampling
#Sites EL and CL have the most years of sampling (8 and 7 years)
n_yr_per_site<-summarise(group_by(dat,site),n_years=length(unique(Year)))
n_yr_per_site
ggplot(n_yr_per_site)+
  geom_bar(aes(x=n_years))+
  labs(y="Number of Sites",x="Number of Years Sampled")+
  theme_bw()


#Which years sampled the most, for which sites?
#seems like most sites with longer monitoring efforts began in 2011 and continue through 2020
#Most sites sampled only once occured before 2010, but a few one-off sites were also sampled 2011-2020
ggplot(left_join(distinct(dat[,c("Year","site")], .keep_all = T),n_yr_per_site,by="site"))+
  geom_bar(aes(x=Year,fill=as.factor(n_years)))+
  scale_fill_viridis_d()+
  labs(fill="N Years\nSite Sampled", x="Year",y="Number of Sites Sampled")+
  theme_bw()

n_nest_per_siteyr<-summarise(group_by(filter(dat,bp!="b"),site,Year),n_nests=n())
arrange(n_nest_per_siteyr,desc(n_nests))

# 2011-2015 has the most sites continuously sampled




# Thin to balance nest failures 
# create a RasterLayer with the extent of nest locations
points<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_dist_err_removed.shp"))
r <- rast(points)

# project to UTM so we can specify resolution in meters
r<-terra::project(r,"EPSG:26918")
nests_m<-st_transform(points,"EPSG:26918")

# set the resolution of the cells to 30 m
res(r) <- 30
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, ext(r)+10)

# thin nest fails down to 1 nest per 30 meters
fails<-nests_m%>%
  filter(fate==0)

non_fails<-nests_m%>%
  filter(fate!=0|is.na(fate))

nest_thin <- gridSample(fails, r, n=1)

fails$Easting<-sf::st_coordinates(fails)[,1]
fails$Northing = sf::st_coordinates(fails)[,2]
fails_ss<- filter(fails,(Easting%in%nest_thin[,1]&Northing%in%nest_thin[,2]))%>%
  distinct(geometry,.keep_all = T)

#number of nest records removed N=687
nrow(fails)-nrow(fails_ss)

#recombine subset fails and fledges
nests_ss<-rbind(fails_ss%>%dplyr::select(-Easting,-Northing),non_fails)


#balance survival data by filtering for the thinned record ids
surv_dat<-surv_dat%>%
  filter(fate==1|id%in%fails_ss$id)

table(surv_dat$fate);sum(surv_dat$fate)/nrow(surv_dat) # get close to 0.5 prevelance 



# thin the veg data if using as background points
if(ab_type=="v"){
  
  # create a RasterLayer with the extent of veg locations
  points<-st_read(paste0(path_out,"Final_outputs/veg_locations/veg_locations_12_29_22.shp"))%>%
    filter(PontTyp=="Random")
  r <- rast(points)
  
  # project to UTM so we can specify resolution in meters
  r<-terra::project(r,"EPSG:26918")
  veg_m<-st_transform(points,"EPSG:26918")
  
  # set the resolution of the cells to 30 m
  res(r) <- 30
  # expand (extend) the extent of the RasterLayer a little
  r <- extend(r, ext(r)+10)
  
  # thin veg points down to 1 point per 30 meters
  fails<-veg_m
  veg_thin <- gridSample(fails, r, n=1)
  
  fails$Easting<-sf::st_coordinates(fails)[,1]
  fails$Northing = sf::st_coordinates(fails)[,2]
  veg_ss<- filter(fails,(Easting%in%veg_thin[,1]&Northing%in%veg_thin[,2]))%>%
    distinct(geometry,.keep_all = T)
  
  #number of nest records removed N=293
  nrow(fails)-nrow(veg_ss)
  
  #balance survival data by filtering for the thinned record ids
  pres_dat<-pres_dat%>%
    filter(id%in%veg_ss$veg_id | bp=="p")
  
  table(pres_dat$presence);sum(pres_dat$presence)/nrow(pres_dat) 
  
}



## Summary of data availability after filtering and thinning/balancing
#-----------------------------------------------------------------------------------
nrow(pres_dat) 
#There are 4728 nest observations
nrow(surv_dat)
#There are 1611 nest fate observations

length(unique(pres_dat$site))
#32 sites




## Divide into testing and training data (80 train/20 test) - k-fold data partitioning
#---------------------------------------------------------------------------------------------

# if doing presence absence, split both. If doing presence only (bioclim), only split the presence data, background data only used for testing

#number of partitions to make. eg, 5 splits into 5 groups each of which will be used once as a test dataset (20% testing 80% training)
k<-5
set.seed(123)
pres_dat$group<-kfold(pres_dat,k)
surv_dat$group<-kfold(surv_dat,k)


#check balance in each split dataset
table(pres_dat$presence,pres_dat$group) #retains same balance
table(surv_dat$fate,surv_dat$group)



# set a constant response variable name
pres_dat$y<-pres_dat$presence

surv_dat$y<-surv_dat$fate

#*try using background points for success models***
#t<-rbind(pres_dat[pres_dat$bp=="b",],surv_dat[surv_dat$fate==1,])
#t_f<-t[t$y==0,]
#t_s<-t[t$y==1,]
#t_sum<-summarize(group_by(t_s,region),n=n())
#t_f<-rbind(t_f[sample(which(t_f$region==1),as.numeric(t_sum[1,2])),],
#                    t_f[sample(which(t_f$region==3),as.numeric(t_sum[2,2])),],
#                    t_f[sample(which(t_f$region==4),as.numeric(t_sum[3,2])),],
#                    t_f[sample(which(t_f$region==5),as.numeric(t_sum[4,2])),])
#surv_dat<-rbind(t_f,t_s)
#***

table(pres_dat$y)
table(surv_dat$y)
