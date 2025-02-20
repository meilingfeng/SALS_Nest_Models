
library(tidyverse)

##  file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


## all species to run through
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP")


## calculate nest initiation dates and days since the last spring tide that they were initiated


# 1. read in spring tide dates
#-----------------------------------
tide<-read.csv(paste0(dat_path,"Environmental Predictors/spring_tide_dates_2011_2024.csv"))%>%
  filter(Month>=3&Month<=9)%>%
  transmute(date=ymd(paste(Year,Month,substr(MaxDate,7,8),sep = "/")))%>%
  #create time intervals between spring tides
  arrange(date)%>%
  mutate(lag=lead(date))%>%
  filter(!is.na(lag))%>%
  mutate(int=interval(ymd(date),ymd(lag))%/%months(1))%>%
  #make sure the intervals don't span between breeding seasons
  filter(int<=1)%>%
  dplyr::select(-int)
  

  # create a list of dates between each spring tide date 
tide.dates<-list()
for(i in 1:nrow(tide)){
  temp<-seq(tide[i,]$date, tide[i,]$lag, by="days")
  # pair each date with each spring tide date
  tide.dates[[i]]<-data.frame(spring.tide.date=rep(tide[i,]$date,length(temp)),date=temp)
}
tide.dat<-do.call("rbind",tide.dates)


tide.dat<-tide.dat%>%
  # subtract the spring tide date from all "between" dates to get days since the last spring tide for each date
  mutate(time_since_tide=date-spring.tide.date)%>%
  distinct(date,.keep_all = T)%>%
  mutate(time_since_tide=ifelse(date%in%tide$date,0,time_since_tide))%>%
  dplyr::select(-spring.tide.date)%>%
  mutate(time_since_tide=as.numeric(time_since_tide))





# 2. Calculate nest initiation dates using nest monitoring data
#---------------------------------------------------------------------
first_egg_dates<-read.csv(paste0(dat_path,"Demographic Database/Nest_Monitoring_Feb25.csv"))%>%
  filter(Visit.Date!="")

# a) get the discovery dates for each nest
nest_discovery<-first_egg_dates%>%
  arrange(SHARP.Nest,Visit.Date)%>%
  group_by(SHARP.Nest)%>%
  summarise(across(starts_with(c("Visit.Date","Eggs","Nestlings","Chick.Age")),first))

# b) get total eggs observed at the nest over nest monitoring period from nest fate data and add to the monitoring data
fates<-read.csv(paste0(dat_path,"Demographic Database/NestFates_2001-2024.csv"))%>%
  dplyr::select(SHARP.Nest=SHARPNestID,Max.Num.Eggs,fate=UltimateNestFate)%>%
  distinct(SHARP.Nest,.keep_all = T)

nest_discovery<-nest_discovery%>%
  left_join(fates,by="SHARP.Nest")

# c) get site and year from cleaned nest location data, remove records without valid coordinate info

nests<-read.csv(paste0(path_out,"Final_outputs/Nest_locations/corrected_nest_coords_01_29_25.csv"))%>%
  dplyr::select(SHARP.Nest=id,Species,site.code,Year,missing.location.rec,missing.coords,out.bounds,wrong.site,iso.rec)

#   start looping through species here (different species may have different nest discovery times)
for (j in 1:length(speciesnames)){

nest_discovery2<-nest_discovery%>%
  left_join(nests,by="SHARP.Nest")%>%
  filter(missing.location.rec==0&missing.coords==0,out.bounds==0,wrong.site==0,iso.rec==0)%>%
  filter(Species==speciesnames[j])%>%
  dplyr::select(-c("missing.location.rec","missing.coords","out.bounds","wrong.site","iso.rec","Eggs.Warm","Species"))

# d) replace all missing chick age, egg, and nestling count values with 0
nest_discovery2[,3:21]<-replace(nest_discovery2[,3:21],
                               nest_discovery2[,3:21]==""|nest_discovery2[,3:21]=="NOT REC"|nest_discovery2[,3:21]=="na"|is.na(nest_discovery2[,3:21]),
                               "0")
nest_discovery2<-nest_discovery2%>%
  mutate(across(3:21,as.numeric))

# e) mark at what stage a nest was found (eggs, chicks, or neither)
nest_discovery3<-nest_discovery2%>%  
  mutate(discovery_status=case_when(
    #When there are some but not all eggs in the nest on the discovery date
    Eggs>0&Eggs<Max.Num.Eggs~"laying",
    #When there are nestlings on the discovery date
    Nestlings>0~"hatched",
    #When the nest was discovered during incubation
    Nestlings==0&Eggs>=Max.Num.Eggs&Eggs!=0~"incubation",
    #When the nest was discovered during nest construction
    Nestlings==0&Eggs==0&Max.Num.Eggs>0~"construction",
    #When there were no eggs or nestlings observed in the nest when discovered or monitored (inactive)
    Nestlings==0&Eggs==0&Max.Num.Eggs==0~"inactive"))


# f) calculate nest initiation dates (first egg date minus 4 days, per Shriver et al. 2007)
nest_discovery4<-nest_discovery3%>%
  rowwise()%>%
  mutate(
    # calculate chick age for nests found during hatching
    max_age=as.numeric(max(c_across(contains("Chick.Age")),na.rm = T)),
    initiation_date=case_when(
      # for nests discovered during laying, subtract 1 day per egg to get first egg date
      discovery_status=="laying"~mdy(Visit.Date)-(Eggs+4),
      # for nests discovered during hatching, subtract age of chicks, an 11 day incubation period (inclusive of hatch day), and 1 day per egg to get first egg date
      discovery_status=="hatched"~mdy(Visit.Date)-(max_age+11+Eggs+Nestlings+4),
      # for nests discovered under construction, assign the discovery date as the initiation date
      discovery_status=="construction"~mdy(Visit.Date),
      # for inactive nests or nests discovered during incubation, assign discovery date for now
      discovery_status%in%c("incubation","inactive")~mdy(Visit.Date))
  )

      # for nests discovered during incubation, assign them the average time between initiation and discovery, per site and year 
      # (to account for dependency on crew)
avg_discovery_time<-nest_discovery4%>%
  filter(!discovery_status%in%c("incubation","inactive"))%>%
  mutate(disc_time=interval(ymd(initiation_date),mdy(Visit.Date))%/%days(1))%>%
  group_by(site.code,Year)%>%
  summarise(avg_disc_time=floor(mean(disc_time)))

nest_discovery5<-left_join(nest_discovery4,avg_discovery_time,by=c("site.code","Year"))%>%
  filter(!is.na(avg_disc_time))%>%
  mutate(initiation_date2=ifelse(discovery_status%in%c("incubation","inactive"),ymd(initiation_date)-avg_disc_time,ymd(initiation_date)),
         initiation_date2=as_date(initiation_date2))
      # any nests found during incubation that do not have a site x year average, just back count the number of eggs like the "laying" observations
nest_discovery6<-left_join(nest_discovery4,avg_discovery_time,by=c("site.code","Year"))%>%
  filter(is.na(avg_disc_time))%>%
  mutate(initiation_date2=ymd(initiation_date)-(Eggs+4))

nest_discovery7<-rbind(nest_discovery5,nest_discovery6)%>%
  ungroup()%>%
  dplyr::select(id=SHARP.Nest, init_date=initiation_date2)



#then add the associated time since moon based on each first egg date
time_since_tide<-nest_discovery7%>%
  left_join(tide.dat,by=c("init_date"="date"))%>%
  distinct(id,.keep_all = T)


#write to file

write.csv(time_since_tide,paste0(path_out,"Intermediate_outputs/Tides/time_since_tide_",speciesnames[j],".csv"),row.names = F)
}
