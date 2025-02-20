
## 1. Plot the data
#---------------------------------------------------------------------------------
## a. Check spatial/temporal sampling bias

# look at spatial distribution

# How are observations distributed across marsh sites?
n_nest_per_site<-summarise(group_by(surv_dat,site),n_nests=n())%>%filter(!(is.na(site)))
hist(n_nest_per_site$n_nests) 
length(n_nest_per_site[n_nest_per_site$n_nests<101,]$n_nests)/nrow(n_nest_per_site)
#76% of sites have 100 or fewer nest observations

#create lat bins and look at sampling intensity across latitude
n_nest_per_site<-surv_dat%>%
  mutate(lat_bin = cut(latitude, breaks=8))%>%
  dplyr::select(site,lat_bin)%>%
  distinct()%>%
  right_join(n_nest_per_site,by="site")
ggplot(n_nest_per_site, aes(x=lat_bin, y=n_nests))+
  geom_boxplot()+
  labs(y= "Number of nest observations", x="Geographic Region")+
  geom_jitter(width = 0.2)
# some difference in nest observations across latitude, less nests near range edges.


# look at temporal distribution across sites

#16 sites have at least 5 years of sampling
#Site AT, OC, EL, CL all have 8+ years of sampling
#Sites HM and ER have the most years of sampling (12 years)
n_yr_per_site<-summarise(group_by(surv_dat,site),n_years=length(unique(Year)))
n_yr_per_site
ggplot(n_yr_per_site)+
  geom_bar(aes(x=n_years))+
  labs(y="Number of Sites",x="Number of Years Sampled")+
  theme_bw()


#Which years sampled the most, for which sites?
#most sites covered between 2011-2015, more single year sites appear past 2020
ggplot(left_join(distinct(surv_dat[,c("Year","site")], .keep_all = T),n_yr_per_site,by="site"))+
  geom_bar(aes(x=as.factor(Year),fill=as.factor(n_years)))+
  scale_fill_viridis_d()+
  labs(fill="N Years\nSite Sampled", x="Year",y="Number of Sites Sampled")+
  theme_bw()

n_nest_per_siteyr<-summarise(group_by(surv_dat,site,Year),n_nests=n())
arrange(n_nest_per_siteyr,desc(n_nests))



# b. Visualize some assumptions we have about the system
# assume nest site selection is not changing over time - short time period and a bird with high site fidelity
ggplot(pres_dat[pres_dat$y==1,],aes(x=as.factor(Year),y=ndvi))+#test different habitat variables
  geom_boxplot()+
  geom_jitter(aes(color=site),alpha=0.3)
#there appears to be some differences over the years, but mostly due to changes in which sites were sampled those years

# assume that later into the spring tide cycle results in more nest flooding, depending on the height of the spring tide too
ggplot(surv_dat, aes(y=y,x=time_since_tide,group=site))+
  geom_point(aes(color=site))+
  geom_smooth(aes(color=site),method="lm",se=F)
#generally yes for many sites, predation and low tides may also be introducting noise

# also expect that sparrows nest closer to new moon later in the season as they become syncronized, but this may depend on the height of MHHWs
plot_dat<-surv_dat%>%mutate(MonthDay = format(as.Date(init_date), "%m-%d"))%>%
  arrange(MonthDay,Year)
ggplot(plot_dat,aes(x=MonthDay,y=time_since_tide,group=site))+
  geom_point(aes(color=site))+
  geom_smooth(aes(color=site),method = "lm",se=F)


# assume that site selection and survival may change between sites
# plot proportion of fledges (nest success) across sites (is nesting success different across sites?)
success_site<-surv_dat%>%group_by(site)%>%
  filter(!is.na(fate))%>%
  summarize(pr_success=sum(fate,na.rm=T)/n(),
            lat=mean(latitude,na.rm=T))

hist(success_site$pr_success, breaks=seq(0,1,0.1))
plot(success_site$pr_success~success_site$lat,xlab="Site latitude",ylab="Proportion of nests fledging")
#nest success appears to be more variable at higher latitudes, but similar means
mean(success_site$pr_success)
#nest success across sites is 34% on average


#plot relationships with covariates
ggplot(pres_dat,aes(x=LOMARSH,y=y))+#test different habitat variables
  geom_point()+
  geom_smooth(method="gam",se=F)+
  facet_wrap(~Region)

# day of year influences site selection (effects of habitat variables)
plot_dat<-pres_dat%>%mutate(Month = format(as.Date(init_date), "%m"))%>%
  arrange(Month,Year)
ggplot(plot_dat,aes(x=LOMARSH,y=y,group=Month))+#test different habitat variables
  geom_point(aes(color=Month))+
  geom_smooth(aes(color=Month),method="gam",se=F)
#pca has some weird looking zero values
#


#overall, different sites (locations) appear to have different habitat selection, but selection at sites appears constant over time
# account for latitude but not year in nest presence model
# relationships between nest presence and habitat variables change over the breeding season months- maybe as vegetation grows in?
# account for the interaction between day of year and vegetation related habitat variables.


#most sites show that nest success decreases as time from the last spring tide increases.
#include time since moon as a variable for nest success models
#also include year to see if success is changing over time
#latitude because mean nest success varies across sites. 


