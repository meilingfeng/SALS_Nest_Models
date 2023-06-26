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



## 3. Check for balance among sample groups  #more important for timeseries modeling
# Remove Years with less than 5 observations
yr_n<-summarise(group_by(pres_dat,Year),count=n())

#remove sites with less than 5 observations
site_n<-summarise(group_by(pres_dat,site),count=n())

pres_dat<-filter(pres_dat,(is.na(site)|pres_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|pres_dat$Year%in% yr_n[yr_n$count>=5,]$Year))
surv_dat<-filter(surv_dat,(is.na(site)|surv_dat$site %in% site_n[site_n$count>=5,]$site)|(is.na(Year)|surv_dat$Year%in% yr_n[yr_n$count>=5,]$Year))
