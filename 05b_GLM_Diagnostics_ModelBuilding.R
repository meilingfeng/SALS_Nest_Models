library(tidyverse)
library(car)#vif
library(PerformanceAnalytics)#correlation plots
library(patchwork)
library(sf) #spatial data
library(gstat) #variograms
library(tseries) #autocorrelation
library(lme4)


### Set up
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

pres_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_pres_dat.csv"))
surv_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_fate_dat.csv"))



## 2. Look for outliers 
#-------------------------------------------------------------------
#list all predictor variables in the full model
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 
pres_preds<-pres_dat[,c(all_terms,"doy")]
surv_preds<-surv_dat[,c(all_terms,"doy","time_since_tide","Year")]


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
car::vif(mod.p)
mod.s<-glm(as.formula(paste0("y~",paste(c(all_terms,"time_since_tide","Year"),collapse = "+"))), surv_dat,family=binomial(link="logit"))
car::vif(mod.s)
#keeping all variables seems fine, no VIFs above 3



## 4. Should we consider interactions?
#---------------------------------------------------
# Hypothesis: The quality of high marsh varies. 
# Interactions: Areas with more high marsh AND 
              # higher elevation
              # brighter PCA = less inundated
              # higher NDVI (healthier veg)
# within the high marsh, do locations with (successful) nests have higher NDVI, etc?

#high marsh used for nesting seems to have distinct characteristics from the rest of the high marsh area
# adding more high marsh habitat also adds higher PCA values for nesting sites, sparrows are picking high marsh habitat with certain reflective qualities.
p1<-ggplot(pres_dat,aes(x=HIMARSH,y=pca,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Site", linetype= "Nest Site",y = "Surface Brightness", x="Proportion High Marsh") + 
  theme_classic(base_size = 12)

# nests are in high marsh with higher elevation
p2<-ggplot(pres_dat,aes(x=HIMARSH,y=elevation,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Site", linetype= "Nest Site",y = "Elevation", x="Proportion High Marsh") + 
  theme_classic(base_size=12)

# nests are in high marsh with higher NDVI
p3<-ggplot(pres_dat,aes(x=HIMARSH,y=ndvi,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Site", linetype= "Nest Site",y = "NDVI", x="Proportion High Marsh") + 
  theme_classic(base_size=12)



(p1/p2/p3)+plot_annotation(tag_levels = 'A')+
  plot_layout(axis_titles = "collect")
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_quality_interactions.png"), width = 5, height = 8, dpi = "retina")

# high marsh has higher PCA (brightness) for nesting areas and higher elevation
  # plots suggest that quality of high marsh is dependent on these other factors.



## 5. Are relationships between response and predictors linear?
#----------------------------------------------------------------------------------

# for all models, assume relationships between phenological vegetation metrics change with DoY
# for surival models, assume survival decreases with increasing time since last spring tide, also test whether survival is changing over time (Year)

# make ndvi and elevation quadratic and center time variables and latitude
pres_dat<-pres_dat%>%
  mutate(doy_c=scale(doy,center = T, scale = T),
         latitude_c=scale(latitude,center = T, scale = T),
         pca_c=scale(pca,center = T, scale = T),
         #log tideres, low marsh, ndvi, uvvrmean, right skewed data
         log.ndvi=log10(ndvi+0.0001),
         log.tideres=log10(tideres+0.0001),
         log.LOMARSH=log10(LOMARSH+0.0001),
         log.uvvr_mean=log10(uvvr_mean+0.0001))

surv_dat<-surv_dat%>%
  mutate(doy_c=scale(doy,center = T, scale = T),
         latitude_c=scale(latitude,center = T, scale = T),
         pca_c=scale(pca,center = T, scale = T),
         # also apply Year and time since spring tide to survival
         time_since_tide_c=scale(time_since_tide,center = T, scale = T),
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


#square term for success high marsh?

## 6. Build candidate models and check they meet assumptions of a binomial GLM
#-----------------------------------------------------------------------------
## Start a list of potential models
mod_list_pres<-list()
mod_list_surv<-list()


#NULL MODELS
mod_list_pres[[1]]<-glm(y~latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[1]]<-glm(y~latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))

#ELEVATION MODELS
# does vegetation matter or is elevation enough to signal where nests can avoid flooding?
mod_list_pres[[2]]<-glm(y~poly(elevation,2,raw = T)+latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[2]]<-glm(y~poly(elevation,2,raw=T)+latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))


#VEG COMMUNITY MODELS
# Do marsh vegetation community classifications best describe nesting habitat?
# compare a model using just proportion of high marsh + low marsh + latitude
# low marsh is also included to give context to where high marsh is edge of the flood zone vs closer to the marsh boundary
# marsh zones often used to delineate suitable habitat - based on elevation relative to tidal amplitude and vegetation communities
mod_list_pres[[3]]<-glm(y~HIMARSH*LOMARSH+latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[3]]<-glm(y~HIMARSH*LOMARSH+latitude_c+Year_c+time_since_tide_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))

#WITHIN-MARSH MODELS
# Does nest habitat quality vary within high marsh?
# elevation
mod_list_pres[[4]]<-glm(y~HIMARSH*elevation+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))

mod_list_surv[[4]]<-glm(y~HIMARSH*elevation+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# brightness
mod_list_pres[[5]]<-glm(y~HIMARSH*pca+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[5]]<-glm(y~HIMARSH*pca+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# NDVI
mod_list_pres[[6]]<-glm(y~HIMARSH*poly(log.ndvi,2,raw=T)+HIMARSH*log.LOMARSH+latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[6]]<-glm(y~HIMARSH*poly(log.ndvi,2,raw=T)+HIMARSH*log.LOMARSH+latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# all within-marsh characteristics (additive)
mod_list_pres[[7]]<-glm(y~poly(elevation,2,raw = T)+
                          poly(log.ndvi,2,raw=T)+
                          pca_c+
                          HIMARSH*LOMARSH+
                          latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[7]]<-glm(y~poly(elevation,2,raw=T)+
                          poly(log.ndvi,2,raw=T)+
                          pca_c+
                          HIMARSH*LOMARSH+
                          latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))




#BETWEEN-MARSH MODELS
mod_list_pres[[8]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+HIMARSH*LOMARSH+
                          latitude_c+doy_c, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv[[8]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+HIMARSH*LOMARSH+
                          latitude_c+Year_c+time_since_tide_c+doy_c, 
                        data=surv_dat,
                        family = binomial(link="logit"))



#GLOBAL MODEL
mod_list_pres[[9]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+
                          poly(elevation,2,raw = T)+poly(log.ndvi,2,raw=T)+pca_c+
                          HIMARSH*LOMARSH+
                          latitude_c+doy_c, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv[[9]]<-glm(y~log.uvvr_mean+uvvr_diff+log.tideres+
                          poly(elevation,2,raw=T)+poly(log.ndvi,2,raw=T)+pca_c+
                          HIMARSH*LOMARSH+
                          latitude_c+Year_c+time_since_tide_c+doy_c, 
                         data=surv_dat,
                         family = binomial(link="logit"))


#Diagnostics
mod<-mod_list_surv[[1]]
  #overdispersed?
chisq <- sum(resid(mod, type='pearson')^2)
chisq/df.residual(mod) 

#par(mfrow=c(2,2));plot(mod) # normally distributed, equal variance, independence?



#get residuals
pres_dat$res<-rstandard(mod_list_pres[[9]]) 
surv_dat$res<-rstandard(mod_list_surv[[9]])
pres_dat$prob<-mod_list_pres[[9]]$fitted.values
surv_dat$prob<-mod_list_surv[[9]]$fitted.values



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
dat2.2<-dat2%>%group_by(Year)%>%summarise(res=mean(res))%>%
  arrange(Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))

dat3.2<-dat3%>%group_by(Year)%>%summarise(res=mean(res))%>%
  arrange(Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))



ggplot(filter(dat3.2,!(is.na(lag1))),aes(x=res,y=lag1))+
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



