library(tidyverse)
library(car)#vif
library(PerformanceAnalytics)#correlation plots
library(patchwork)
library(sf) #spatial data
library(gstat) #variograms
library(tseries) #autocorrelation
library(glmmTMB)
library(DHARMa)
library(lme4)

### Set up
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

pres_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_pres_dat.csv"))
surv_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_surv_dat.csv"))


## 1. Plot the data
#---------------------------------------------------------------------------------
## a. Check spatial/temporal sampling bias

# look at spatial distribution

  # How are observations distributed across marsh sites?
n_nest_per_site<-summarise(group_by(surv_dat,site),n_nests=n())%>%filter(!(is.na(site)))
hist(n_nest_per_site$n_nests) 
length(n_nest_per_site[n_nest_per_site$n_nests<251,]$n_nests)/nrow(n_nest_per_site)
  #96% of sites have 250 or less nest observations

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
  # some difference in sampling intensity throughout region, but nothing seems overly under sampled
  # lat 42 has the least sites


# look at temporal distribution across sites

  #16 sites have at least 5 years of sampling
  #Site AT, ER, HM, ID, MN, OC all have 6 years of sampling
  #Sites EL and CL have the most years of sampling (8 and 7 years)
n_yr_per_site<-summarise(group_by(surv_dat,site),n_years=length(unique(Year)))
n_yr_per_site
ggplot(n_yr_per_site)+
  geom_bar(aes(x=n_years))+
  labs(y="Number of Sites",x="Number of Years Sampled")+
  theme_bw()


  #Which years sampled the most, for which sites?
  #seems like most sites with longer monitoring efforts began in 2011 and continue through 2020
  #Most sites sampled only once occured before 2010, but a few one-off sites were also sampled 2011-2020
ggplot(left_join(distinct(surv_dat[,c("Year","site")], .keep_all = T),n_yr_per_site,by="site"))+
  geom_bar(aes(x=as.factor(Year),fill=as.factor(n_years)))+
  scale_fill_viridis_d()+
  labs(fill="N Years\nSite Sampled", x="Year",y="Number of Sites Sampled")+
  theme_bw()

n_nest_per_siteyr<-summarise(group_by(filter(surv_dat,bp!="b"),site,Year),n_nests=n())
arrange(n_nest_per_siteyr,desc(n_nests))

  # 2011-2015 has the most sites continuously sampled




# b. Visualize some assumptions we have about the system
  # assume nest site selection is not changing over time - short time period and a bird with high site fidelity
ggplot(pres_dat[pres_dat$y==1,],aes(x=as.factor(Year),y=ndvi))+#test different habitat variables
  geom_boxplot()+
  geom_jitter(aes(color=site),alpha=0.3)
  #there appears to be some differences over the years, but mostly due to changes in which sites were sampled those years

  
  # assume that later into the moon cycle results in more nest flooding, depending on the height of the spring tide too
ggplot(surv_dat, aes(y=y,x=time_since_moon,group=site))+
  geom_point(aes(color=site))+
  geom_smooth(aes(color=site),method="lm",se=F)
  #generally yes for many sites, predation and low tides may also be introducting noise

  # also expect that sparrows nest closer to new moon later in the season as they become syncronized, but this may depend on the height of MHHWs
plot_dat<-surv_dat%>%mutate(Month = format(as.Date(date), "%m"))%>%
  arrange(Month,Year)
ggplot(plot_dat,aes(x=Month,y=time_since_moon))+
  geom_boxplot()+
  geom_jitter(aes(color=site),alpha=0.3)+
  facet_wrap(~Year)

plot_dat<-surv_dat%>%mutate(MonthDay = format(as.Date(date), "%m-%d"))%>%
  arrange(MonthDay,Year)
ggplot(plot_dat,aes(x=MonthDay,y=time_since_moon,group=site))+
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
#nest success appears to be more variable at higher latitudes, but similar on average
mean(success_site$pr_success)
#nest success across sites is 50% on average

ggplot(pres_dat,aes(x=elevation,y=y))+#test different habitat variables
  geom_point()+
  geom_smooth(method="gam",se=F)+
  facet_wrap(~region)

# day of year influences site selection (effects of habitat variables)
plot_dat<-pres_dat%>%mutate(Month = format(as.Date(date), "%m"))%>%
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




## 2. Look for outliers 
#-------------------------------------------------------------------
#list all predictor variables in the full model
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 
pres_preds<-pres_dat[,c(all_terms,"doy")]
surv_preds<-surv_dat[,c(all_terms,"time_since_moon","Year")]


#Boxplot, histogram, scatterplot are graphical tools for outlier detection
for(i in 1:ncol(surv_preds)){
  par(mfrow = c(1, 3),mar=c(5,4,7,2)) #bottom, left, top, right
  hist(surv_preds[,i], main = "Histogram")
  boxplot(surv_preds[,i], main = "Boxplot")
  qqnorm(surv_preds[,i], main = "Normal Q-Q plot")
  mtext(paste0(all_terms[i]), side = 3, line = - 2, outer = TRUE)
}
dev.off()
#sqrt or log tideres, low marsh, ndvi, uvvrmean?




## 3. Is there collinearity among the covariates?
#-------------------------------------------------------------
# a) plot pairwise scatterplots, (alternatively look at correlation coefficients, or a PCA biplot of all covariates)
#PerformanceAnalytics::chart.Correlation(pres_preds, histogram = TRUE, method = "spearman")
#elevation, pca, and ndvi correlated ~0.6 max

# b) Look at variance inflation factor (larger if the covariate R2 is large, most variation in the covariate is explained by variation in all the other covariates)
# remove variables with the highest VIF until all remaining variable have below a select threshold
# even moderate collinearity can mask weak ecological signals so using a low threshold like 3 or 2 is important.
mod.p<-glm(as.formula(paste0("y~",paste(c(all_terms,"doy"),collapse = "+"))), pres_dat,family=binomial(link="logit"))
vif(mod.p)
mod.s<-glm(as.formula(paste0("y~",paste(c(all_terms,"time_since_moon","Year"),collapse = "+"))), surv_dat,family=binomial(link="logit"))
vif(mod.s)
#keeping all variables seems fine, no VIFs above 3



## 4. Should we consider interactions?
#---------------------------------------------------
# Hypothesis: The quality (resiliency) of high marsh is more attractive for nest building and leads to more successful nests. 
# Interactions: Areas with more high marsh AND 
              # lower mean UVVR (more vegetated/resilient)
              # 0-negative UVVR change (stable resilience or improved vegetation cover)
              # certain value of reflectance that indicates an unmeasured habitat characteristic
              # higher NDVI (healthier veg)
# within the high marsh, do locations with (successful) nests have higher NDVI, etc?

#high marsh used for nesting seems to have distinct characteristics from the rest of the high marsh area
# adding more high marsh habitat also adds higher PCA values for nesting sites, sparrows are picking high marsh habitat with certain reflective qualities.
p1<-ggplot(pres_dat,aes(x=HIMARSH,y=pca,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Site", linetype= "Nest Site",y = "Raw Reflectance (PCA)", x="Proportion High Marsh") + 
  theme_classic(base_size = 12)

# successful nests are in high marsh with higher levels of NDVI
p2<-ggplot(pres_dat,aes(x=HIMARSH,y=tideres,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Site", linetype= "Nest Site",y = "Tidal Restriction Index", x="Proportion High Marsh") + 
  theme_classic(base_size=12)

# successful nests are in high marsh with higher levels of NDVI
p3<-ggplot(surv_dat,aes(x=HIMARSH,y=elevation,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Success", linetype= "Nest Success",y = "Elevation", x="Proportion High Marsh") + 
  theme_classic(base_size=12)

#p3+(p1/p2)
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_quality_interactions_",ab_type,".png"), width = 10, height = 8, dpi = "retina")



# conclusions: high marsh is selected nesting habtiat but the suitable range for successful nesting lies within the range of possible elevations in the high marsh.
  # high marsh * elevation for success
# high marsh has higher PCA (brightness) for nesting areas and lower tidal restriction
  # interaction between high marsh and these variables suggest that quality of high marsh is dependent on these other factors.



## 5. Are relationships between response and predictors linear?

# to all models, assume relationships between phenological vegetation metrics change with DoY
# for surival models, assume survival decreases with increasing time since last new moon, also test whether survival is changing over time (Year)

# make ndvi and elevation quadratic and center time variables and latitude
pres_dat<-pres_dat%>%
  mutate(doy_c=scale(doy,center = T, scale = T),
         latitude_c=scale(latitude,center = T, scale = T),
         #log tideres, low marsh, ndvi, uvvrmean, right skewed data
         log.ndvi=log10(ndvi+0.0001),
         log.tideres=log10(tideres+0.0001),
         log.LOMARSH=log10(LOMARSH+0.0001),
         log.uvvr_mean=log10(uvvr_mean+0.0001))

surv_dat<-surv_dat%>%
  mutate(doy_c=scale(doy,center = T, scale = T),
         latitude_c=scale(latitude,center = T, scale = T),
         # also apply Year and time since moon to survival
         time_since_moon_c=scale(time_since_moon,center = T, scale = T),
         Year_c=scale(Year,center = T, scale = T),
         log.ndvi=log10(ndvi+0.0001),
         log.tideres=log10(tideres+0.0001),
         log.LOMARSH=log10(LOMARSH+0.0001),
         log.uvvr_mean=log10(uvvr_mean+0.0001))

# expect high marsh and low marsh to be linear (postive high, negative low)
# ndvi expect a quadratic relationship.
# less than 0.1 is barren, so expect to be higher than that. Sparse veg like grassland is around 0.2-0.5. Dense veg like shrubs and trees have 0.6-0.9.
# expect uvvr to be negative linear.
# small uvvr means more vegetation cover and resilience. Stable marshes usually at 0-0.15 (Ganju, Couvillion, Defne, and Ackerman 2022)
# expect elevation to be quadratic
# too low is flood prone, too high is depredation prone, also don't nest in uplands.
# let tidal restriction be non-linear
# tidal restriction has been shown to be negative on nesting habitat but can sometimes protect from extreme high tides. 
#let pca be non linear
#
pres_long<-pivot_longer(pres_dat,c(all_of(c(all_terms,"log.ndvi","log.tideres","log.uvvr_mean","log.LOMARSH"))),names_to = "var",values_to = "value")
ggplot(pres_long, aes(x = value, y = y)) +
  geom_point()+
  labs(x="Nest Presence")+
  geom_smooth(method="gam", method.args = list( family = "binomial"))+
  facet_wrap(~var,scales = "free")

surv_long<-pivot_longer(surv_dat,c(all_of(c(all_terms,"log.ndvi","log.tideres","log.uvvr_mean","log.LOMARSH"))),names_to = "var",values_to = "value")
ggplot(surv_long, aes(x = value, y = y)) +
  geom_point()+
  labs(x="Nest Survival")+
  geom_smooth(method="gam", method.args = list( family = "binomial"))+
  facet_wrap(~var,scales = "free")
## suggest a strong non-linearlity, although an odd one                    
ggplot(pres_dat, aes(x = log(ndvi+0.01), y = y)) +
  geom_smooth(method="gam", method.args = list( family = "binomial")) +
  theme_bw()
ggplot(pres_dat, aes(x = ndvi, y = y)) +
  geom_smooth(method="gam", method.args = list( family = "binomial")) +
  theme_bw()

#square term for success high marsh?

## 6. Are we meeting the assumptions of a GLM?
#----------------------------------------------------------------------


## Start a list of potential models
mod_list_pres<-list()
mod_list_surv<-list()


# Do marsh vegetation community classifications best describe nesting habitat?
## compare a model using just proportion of high marsh + low marsh + latitude
### low marsh is also included to give context to where high marsh is
###  - edge of the flood zone vs closer to the marsh boundary
### marsh zones often used to delineate suitable habitat - based on elevation relative to tidal amplitude and vegetation communities
mod_list_pres[[1]]<-glm(y~HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[1]]<-glm(y~HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_moon_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))

## compare a model using high marsh and additional variables +latitude (this is the global model)
mod_list_pres[[2]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+poly(elevation,2)+
                           poly(log.ndvi,2)+pca+HIMARSH*log.LOMARSH+
                           latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[2]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+poly(elevation,2)+
                           poly(log.ndvi,2)+pca+HIMARSH*log.LOMARSH+
                           latitude_c+Year_c+time_since_moon_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))



# Does nest habitat quality vary within high marsh? Does including quality improve the variance explained by veg community classifications?

# although high marsh is broadly representative of suitable nesting vegetation and elevation for flood protection,
# there seems to be a smaller range of elevation within this habitat that is best suited for successful nesting.
mod_list_pres[[3]]<-glm(y~HIMARSH*elevation+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))

mod_list_surv[[3]]<-glm(y~HIMARSH*elevation+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_moon_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# sparrows also may use vegetative habitat cues to find high quality nesting habitat within high marsh (gerardii)
mod_list_pres[[4]]<-glm(y~HIMARSH*pca+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[4]]<-glm(y~HIMARSH*pca+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_moon_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# sparrows have also been shown to occupy high marsh habitat that has been degraded by tidal restriction, this may be a driver of high marsh habitat quality.
# holistically represent characteristics of degraded habitat that are directly attributed to a major driver of habitat loss. 
mod_list_pres[[5]]<-glm(y~HIMARSH*log.tideres+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[5]]<-glm(y~HIMARSH*log.tideres+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_moon_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))





plot(mod_list_pres[[5]]) # 4 seems to have the best diagnostics, most normally distributed, equal variance, independence

plot(mod_list_surv[[4]])



##testing for overdispersion 
chisq<-sum(resid(mod_list_pres[[4]],type='pearson')^2) 
chisq/df.residual(mod_list_pres[[2]]) ##not too bad

1-pchisq(chisq,df=df.residual(mod.p))
#not significant

chisq<-sum(resid(mod_list_surv[[1]],type='pearson')^2) 
chisq/df.residual(mod_list_surv[[1]]) ##not too bad

1-pchisq(chisq,df=df.residual(mod.s))
#not significant




#get residuals
pres_dat$res<-mod_list_pres[[2]]$residuals  
surv_dat$res<-mod_list_surv[[2]]$residuals 
pres_dat$prob<-mod_list_pres[[2]]$fitted.values
surv_dat$prob<-mod_list_surv[[2]]$fitted.values



#Good to check dependence before and after modeling
# plot response against spatial and temporal variables
# check that model residuals do not have dependence structure

#formal way to check is with autocorrelation
#plot auto correlation funcitons (ACF) for regular time series (pearsons correlation between time series and lagged time series)
#plot variograms for irregularly spaced time series and spatial data
# straight diagonal line indicates perfect correlation. Want it to reach a sill where the linear trend devolves

# variograms https://csegrecorder.com/articles/view/the-variogram-basics-a-visual-intro-to-useful-geostatistical-concepts 


# spatial autocorrelation - doesn't seem to be any
dat2<-st_as_sf(pres_dat, coords=c("longitude","latitude"))
dat3<-st_as_sf(surv_dat, coords=c("longitude","latitude"))

gram1<-variogram(res~1,data=dat2)  #filter(dat2,site=="ER")
plot(gram1)
gram2<-variogram(res~1,data=dat3)  #filter(dat2,site=="ER")
plot(gram2)


# temporal autocorrelation
dat2.2<-dat2%>%group_by(site,Year)%>%summarise(res=mean(res))%>%
  arrange(site,Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))

dat3.2<-dat3%>%group_by(site,Year)%>%summarise(res=mean(res))%>%
  arrange(site,Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))


ggplot(filter(dat3.2,!(is.na(lag1))),aes(x=res,y=lag2))+
  geom_point()+
  labs(x="Residuals",y="Residuals 1 Year Lag",title="Nest Survival")+
  theme_classic(base_size=12)
  

ggplot(filter(dat2.2,!(is.na(lag1))),aes(x=res,y=lag1))+
  geom_point()+
  labs(x="Residuals",y="Residuals 1 Year Lag",title="Nest Placement")+
  theme_classic(base_size=12)


#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/temporal_autocorr.png"), width = 10, height = 5, dpi = "retina")


#plot residuals over time and lat
ggplot(filter(pres_dat,!is.na(latitude)),aes(x=latitude,y=res))+
  geom_point()+
  theme_classic(base_size = 12)

ggplot(filter(surv_dat,!is.na(latitude)),aes(x=latitude,y=res))+
  geom_point()+
  theme_classic(base_size = 12)


ggplot(filter(pres_dat,!is.na(Year)&Year!=2010),aes(x=as.factor(Year),y=res))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.3,width=0.2)+
  labs(x="Year",y="residuals")+
  theme_classic(base_size = 12)


ggplot(filter(surv_dat,!is.na(Year)&Year!=2010),aes(x=as.factor(Year),y=res))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.3,width=0.2)+
  labs(x="Year",y="residuals")+
  theme_classic(base_size = 12)


stripchart(surv_dat$res~surv_dat$site,pch=1,method="jitter")
stripchart(pres_dat$res~pres_dat$site,pch=1,method="jitter")



