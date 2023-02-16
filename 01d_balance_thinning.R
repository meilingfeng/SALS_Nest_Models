
library(sf)
library(tidyverse)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(dismo)
source("C:/Users/mefen/OneDrive/Documents/Github/UCONN/SHARP/Functions/gridSample_sf.R")


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## 1. Format Nest Observations
# -------------------------------------

# load nest observation shapefile
nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/nest_locations_01_3_23.shp"))%>%
  # filter records to just SALS or species of interest
  filter(Species=="SALS"&
           # also filter records missing coordinate information or that have coordinate errors
           crd_typ!=1&mssng_l_!=1&mssng_c!=1)%>%
  # convert all coordinates to degrees (NAD83) and create Long and Lat columns
  st_transform("EPSG:4269")%>%
  mutate(Long = sf::st_coordinates(.)[,1],
         Lat = sf::st_coordinates(.)[,2],
         #create binary nest success variable for remaining records
         fate=case_when(fate=="FLEDGED" ~ 1,
                        fate%in%c("FLOODED","DEPREDATED","FAIL UNKNOWN","INACTIVE","NEVER HAD EGGS","NEVER HAD") ~ 0))%>%
  dplyr:: select(id,latitude=Lat,longitude=Long,fate,Year,site=site_cd)%>%
  distinct(id,.keep_all=T)


#general data availability
nrow(nests) 
#There are 3039 nest observations
nrow(nests[!(is.na(nests$fate)),])
#There are 2807 nest fate observations

length(unique(nests$site))
#53 sites

sort(unique(nests$Year))
#2006-2020, we may want to just use data from 2010 and later to assume standardized sampling protocol

nests<-nests%>%filter(Year>2009)
nrow(nests) 
#There are 2476 nest observations since 2010
nrow(nests[!(is.na(nests$fate)),])
#There are 2301 nest fate observations since 2010

length(unique(nests$site))
#31 sites

#write nests from past 2010 to assume consistent sampling
st_write(nests,paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020.shp"))


## 2. Check Class Balance
#--------------------------------------
# does the data need class balancing?
# balancing can help machine learning models

#look at prevelance of nest fledge to nest fail
addmargins(table(nests$fate))
#807/2301 nest successes


nests%>%
  st_drop_geometry()%>%
  summarize(pr_success=sum(fate,na.rm=T)/n())
#nest success relative frequency for the whole dataset is 0.32 (32%)

success_site<-nests%>%group_by(site)%>%
  st_drop_geometry()%>%
summarize(pr_success=sum(fate,na.rm=T)/n(),
          lat=mean(latitude,na.rm=T))

hist(success_site$pr_success, breaks=seq(0,1,0.1))
plot(success_site$pr_success~success_site$lat,xlab="Site latitude",ylab="Proportion of nests fledging")
#nest success appears to be more variable at higher latitudes
mean(success_site$pr_success)
#nest success across sites is also 32% on average



## 3. Check sampling bias
#-----------------------------------------------------
# look at spatial distribution and temporal distribution
# does the data need thinning?
# best to thin the more prevalent class (nest fails to retain higher calibration) - thin to balance

# Are observers returning to locations where they previously saw nests?
# since we are using multi-year data, some nests might be frequently seen in similar locations
# might be good to thin some of these out
# How are observations distributed across marsh sites?
n_nest_per_site<-summarise(group_by(st_drop_geometry(nests),site),n_nests=n())
hist(n_nest_per_site$n_nests) 
sum(n_nest_per_site[n_nest_per_site$n_nests<251,]$n_nests)/nrow(nests)
# 90% of sites have 250 or less nest observations, but a couple have more than 300 observations.

#create lat bins
nests<-nests%>%
  mutate(lat_bin = cut(latitude, breaks=8))
n_nest_per_site<-left_join(n_nest_per_site,distinct(st_drop_geometry(nests[,c("site","lat_bin")])),by="site")
ggplot(n_nest_per_site, aes(x=lat_bin, y=n_nests))+
  geom_boxplot()+
  labs(y= "Number of nest observations", x="Geographic Region")+
  geom_jitter(width = 0.2)
# some difference in sampling intesity throughout region, but nothing seems overly under sampled
# lat 42 has the least sites



#16 sites have at least 5 years of sampling
#Site AT, ER, HM, ID, MN, OC all have 6 years of sampling
#Sites EL and CL have the most years of sampling (8 and 7 years)
n_yr_per_site<-summarise(group_by(st_drop_geometry(nests),site),n_years=length(unique(Year)))
n_yr_per_site
ggplot(n_yr_per_site)+
  geom_bar(aes(x=n_years))+
  labs(y="Number of Sites",x="Number of Years Sampled")+
  theme_bw()


#Which years sampled the most, for which sites?
#seems like most sites with longer monitoring efforts began in 2011 and continue through 2020
#Most sites sampled only once occured before 2010, but a few one-off sites were also sampled 2011-2020
ggplot(left_join(distinct(st_drop_geometry(nests[,c("Year","site")]), .keep_all = T),n_yr_per_site,by="site"))+
  geom_bar(aes(x=Year,fill=as.factor(n_years)))+
  scale_fill_viridis_d()+
  labs(fill="N Years\nSite Sampled", x="Year",y="Number of Sites Sampled")+
  theme_bw()

n_nest_per_siteyr<-summarise(group_by(nests,site,Year),n_nests=n())
n_nest_per_siteyr

#Looks like temporal samples from 2011-2015 might have the most sites continuously sampled




## 4. Thin to balance the nest failures 
#------------------------------------------------------
#** Do we also want to thin all the nest records since we are pooling across time and might have nests built in nearly the same spot (duplicate records)

# create a RasterLayer with the extent of nest locations
r <- rast(nests)

# project to UTM so we can specify resolution in meters
r<-terra::project(r,"EPSG:26918")
nests_m<-st_transform(nests,"EPSG:26918")

# set the resolution of the cells to 1 m
res(r) <- 1
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, ext(r)+10)

# thin nest fails down to 1 nest per meter
fails<-nests_m%>%
  filter(fate==0)

non_fails<-nests_m%>%
  filter(fate!=0|is.na(fate))

nest_thin <- gridSample(fails, r, n=1)

fails$Easting<-sf::st_coordinates(fails)[,1]
fails$Northing = sf::st_coordinates(fails)[,2]
fails_ss<- filter(fails,(Easting%in%nest_thin[,1]&Northing%in%nest_thin[,2]))%>%
  distinct(geometry,.keep_all = T)

#number of nest records removed N=160
nrow(fails)-nrow(fails_ss)

#recombine subset fails and fledges
nests_ss<-rbind(fails_ss%>%select(-Easting,-Northing),non_fails)

# look at which nests were removed
ggplot()+
  geom_sf(aes(color=as.factor(fate)),data=nests_m[nests_m$latitude>41.25&nests_m$latitude<41.3&nests_m$longitude>-72.56,])+
  geom_sf(data=nests_ss[nests_ss$latitude>41.25&nests_ss$latitude<41.3&nests_ss$longitude>-72.56,],color="red",pch="x",cex=2)

#write fail thinned nests from past 2010
st_write(nests_ss,paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_fail_thin_1m.shp"))
