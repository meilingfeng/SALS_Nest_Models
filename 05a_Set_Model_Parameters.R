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
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"





## 2. parameters
# what prediction resolution?
reso<-30

# use random background points or veg plots as absences?
#ab_type<-"v"
ab_type<-"b" 

# pick which predictors and response to use for models
#predictors <- "uvvr"
#predictors <- "no uvvr" #removes observations mising uvvr data


# veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))



## 3. Load observations with predictors
# list of predictor files for each zone as a raster stack
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","ent_txt","cor_txt", "tideres", "uvvr_diff","latitude") #all environmental predictor names

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
         #create binary variable for if nests intersect High Marsh habitat
         Highmarsh=as.factor(ifelse(veg_class=="HIMARSH",1,0)))


#4326


### Filtering data
#--------------------------------------------------------------------------------------------
## 4. Remove nests that are in non-marsh habitat and contain no marsh in their buffer
table(dat$veg_class) # lots of nests in stream and upland classifications
             # remove if...
#dat<-filter(dat, 
             # it's a nest point and not background
            #!((bp=="p")&
             # if point is in non-marsh
 #            !((veg_class%in%c("UPLND","STRM","PHRG","TERRBRD","MUD","POOL")) &
 #            # and nest buffer contains no vegetated marsh habitat
 #              ((HIMARSH!=0 & LOMARSH!=0) | (!is.na(HIMARSH)&!is.na(LOMARSH)))))

table(dat$presence)
#table(dat2$veg_class) 
# (originally 4160, now 3181 records)

# remove rapid demo sites
#don't seem to be any (would be nest without site name)



## 5. Remove missing values from predictors
# only keep nests within the marsh vegetation layer
dat_comp<-dat[complete.cases(dat[,all_terms]),]
# fill in 0's for land cover proportions with NAs 
dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]<-dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]%>%
  replace(is.na(.),0)



## 6. divide into nest presence and nest suvival datasets
pres_dat<-dat_comp
surv_dat<-dat_comp[!is.na(dat_comp$fate),] #filter just the nests with known fates






## 7. Check for balance among sample groups 
# Remove Years with less than 5 observations
yr_n<-summarise(group_by(pres_dat,Year),count=n())

#remove sites with less than 5 observations
site_n<-summarise(group_by(pres_dat,site),count=n())

pres_dat<-filter(pres_dat,(is.na(site)|pres_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|pres_dat$Year%in% yr_n[yr_n$count>=5,]$Year))
surv_dat<-filter(surv_dat,(is.na(site)|surv_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|surv_dat$Year%in% yr_n[yr_n$count>=5,]$Year))


# Check Class Balance
# balancing can help machine learning models - both above 0.2 prevalance 

#look at prevalence of nest fledge to nest fail
addmargins(table(surv_dat$fate));sum(surv_dat$fate)/nrow(surv_dat)

# prevalence of presence to background
addmargins(table(pres_dat$presence)); sum(pres_dat$presence)/nrow(pres_dat)



success_site<-dat%>%group_by(site)%>%
  filter(!is.na(fate))%>%
  summarize(pr_success=sum(fate,na.rm=T)/n(),
            lat=mean(latitude,na.rm=T))

hist(success_site$pr_success, breaks=seq(0,1,0.1))
plot(success_site$pr_success~success_site$lat,xlab="Site latitude",ylab="Proportion of nests fledging")
#nest success appears to be more variable at higher latitudes
mean(success_site$pr_success)
#nest success across sites is also 32% on average



# Check sampling bias

# look at spatial distribution and temporal distribution
# does the data need thinning?
# best to thin the more prevalent class (nest fails to retain higher calibration) - thin to balance

# Are observers returning to locations where they previously saw nests?
# since we are using multi-year data, some nests might be frequently seen in similar locations
# might be good to thin some of these out
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

n_nest_per_siteyr<-summarise(group_by(dat,site,Year),n_nests=n())
n_nest_per_siteyr

#Looks like temporal samples from 2011-2015 might have the most sites continuously sampled




# Thin to balance nest failures (could have nests in nearly the same locations in different years, prioritize the successful nests at these locations)

#** Do we also want to thin all the nest records since we are pooling across time and might have nests built in nearly the same spot (duplicate records)

# create a RasterLayer with the extent of nest locations
points<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_dist_err_removed.shp"))
r <- rast(points)

# project to UTM so we can specify resolution in meters
r<-terra::project(r,"EPSG:26918")
nests_m<-st_transform(points,"EPSG:26918")

# set the resolution of the cells to 5 m
res(r) <- 5
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, ext(r)+10)

# thin nest fails down to 1 nest per 5 meters
fails<-nests_m%>%
  filter(fate==0)

non_fails<-nests_m%>%
  filter(fate!=0|is.na(fate))

nest_thin <- gridSample(fails, r, n=1)

fails$Easting<-sf::st_coordinates(fails)[,1]
fails$Northing = sf::st_coordinates(fails)[,2]
fails_ss<- filter(fails,(Easting%in%nest_thin[,1]&Northing%in%nest_thin[,2]))%>%
  distinct(geometry,.keep_all = T)

#number of nest records removed N=211
nrow(fails)-nrow(fails_ss)

#recombine subset fails and fledges
nests_ss<-rbind(fails_ss%>%dplyr::select(-Easting,-Northing),non_fails)

# look at which nests were removed
ggplot()+
  geom_sf(aes(color=as.factor(fate)),data=nests_m[nests_m$latitude>41.25&nests_m$latitude<41.3&nests_m$longitude>-72.56,])+
  geom_sf(data=nests_ss[nests_ss$latitude>41.25&nests_ss$latitude<41.3&nests_ss$longitude>-72.56,],color="red",pch="x",cex=2)



#balance survival data
surv_dat<-surv_dat%>%
  filter(fate==1|id%in%fails_ss$id)

table(surv_dat$fate);sum(surv_dat$fate)/nrow(surv_dat)





## 8. divide into testing and training data (80 train/20 test) - k-fold data partitioning

# if doing presence absence, split both. If doing presence only (bioclim), only split the presence data, background data only used for testing

#number of partitions to make. eg, 5 splits into 5 groups each of which will be used once as a test dataset (20% testing 80% training)
k<-5
pres_dat$group<-kfold(pres_dat,k)
surv_dat$group<-kfold(surv_dat,k)


#check balance in each split dataset
table(pres_dat$presence,pres_dat$group) #retains same balance
table(surv_dat$fate,surv_dat$group)



# set a constant response variable name
pres_dat$y<-pres_dat$presence

surv_dat$y<-surv_dat$fate







table(pres_dat$y)
table(surv_dat$y)
