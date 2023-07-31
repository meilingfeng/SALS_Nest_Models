library(tidyverse)
library(car)#vif
library(PerformanceAnalytics)#correlation plots
library(patchwork)
library(sf) #spatial data
library(gstat) #variograms
library(tseries) #autocorrelation

### Set up
# -------------------------------------------
source("05a_Set_Model_Parameters.R")


## 1. Look for outliers 
#-------------------------------------------------------------------
#list all predictor variables in the full model
pres_preds<-pres_dat[,all_terms]
surv_preds<-surv_dat[,all_terms]


#Boxplot, histogram, scatterplot are graphical tools for outlier detection
for(i in 1:(length(all_terms))){
  par(mfrow = c(1, 3),mar=c(5,4,7,2)) #bottom, left, top, right
  hist(pres_preds[,i], main = "Histogram")
  boxplot(pres_preds[,i], main = "Boxplot")
  qqnorm(pres_preds[,i], main = "Normal Q-Q plot")
  mtext(paste0(all_terms[i]), side = 3, line = - 2, outer = TRUE)
}
dev.off()
#sqrt or log ndvi and uvvr_mean?




## 2. Is there collinearity among the covariates?
#-------------------------------------------------------------
# a) plot pairwise scatterplots, (alternatively look at correlation coefficients, or a PCA biplot of all covariates)
PerformanceAnalytics::chart.Correlation(pres_preds, histogram = TRUE, method = "spearman")
#entropy, pca, and ndvi correlated
#try removing entropy all_terms[-5]

# b) Look at variance inflation factor (larger if the covariate R2 is large, most variation in the covariate is explained by variation in all the other covariates)
# remove variables with the highest VIF until all remaining variable have below a select threshold
# even moderate collinearity can mask weak ecological signals so using a low threshold like 3 or 2 is important.
mod.p<-glm(as.formula(paste0("y~",paste(all_terms,collapse = "+"))), pres_dat,family=binomial(link="logit"))
vif(mod.p)
mod.s<-glm(as.formula(paste0("y~",paste(all_terms,collapse = "+"))), surv_dat,family=binomial(link="logit"))
vif(mod.s)
#keeping entropy seems fine, no VIFs above 3



## 3. Should we consider interactions?
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
p2<-ggplot(surv_dat,aes(x=HIMARSH,y=ndvi,group=as.factor(y)))+
  geom_point(aes(color=as.factor(y)))+
  geom_smooth(method="lm",aes(linetype=as.factor(y)),color="black")+
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  labs(color = "Nest Success", linetype= "Nest Success",y = "NDVI", x="Proportion High Marsh") + 
  theme_classic(base_size=12)

p1+p2
ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_quality_interactions_",ab_type,".png"), width = 10, height = 5, dpi = "retina")



# conclusions: consider HIMARSH*pca (just presence), HIMARSH*NDVI


## 4. Are response observations independent?
#----------------------------------------------------------------------
# inflates type 1 error or false positives and raises p values in regression

#Good to check dependence before and after modeling
# plot response against spatial and temporal variables
# check that model residuals do not have dependence structure

#formal way to check is with autocorrelation
#plot auto correlation funcitons (ACF) for regular time series (pearsons correlation between time series and lagged time series)
#plot variograms for irregularly spaced time series and spatial data
# straight diagonal line indicates perfect correlation. Want it to reach a sill where the linear trend devolves

# variograms https://csegrecorder.com/articles/view/the-variogram-basics-a-visual-intro-to-useful-geostatistical-concepts 

#get residuals
pres_dat$res<-mod.p$residuals  
surv_dat$res<-mod.s$residuals 
pres_dat$prob<-mod.p$fitted.values
surv_dat$prob<-mod.s$fitted.values

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


p1+p2
ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/temporal_autocorr_",ab_type,".png"), width = 10, height = 5, dpi = "retina")


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


## 5. Are relationships between response and predictors linear?
pres_long<-pivot_longer(pres_dat,c(all_of(all_terms)),names_to = "var",values_to = "value")
ggplot(pres_long, aes(x = value, y = y)) +
  geom_point()+
  labs(x="Nest Presence")+
  geom_smooth(method="gam", method.args = list( family = "binomial"))+
  facet_wrap(~var,scales = "free")

surv_long<-pivot_longer(surv_dat,c(all_of(all_terms)),names_to = "var",values_to = "value")
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
#log ndvi

