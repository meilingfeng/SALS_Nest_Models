

library(tidyverse)
library(sf) #spatial data
library(gstat) #variograms
library(lme4) #mixed models
library(tseries) #autocorrelation


### Set up
# -------------------------------------------
source("05a_Set_Model_Parameters.R") #set up data and model parameters
source("05b_Collinearity_Outliers.R") #remove collinearity in models



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


# a) variograms https://csegrecorder.com/articles/view/the-variogram-basics-a-visual-intro-to-useful-geostatistical-concepts 

#get residuals
pres_train$res<-mod_list_pres[[3]]$residuals  
surv_train$res<-mod_list_surv[[3]]$residuals 

# spatial autocorrelation
dat2<-st_as_sf(surv_train, coords=c("longitude","latitude"))
dat3<-st_as_sf(surv_train, coords=c("longitude","latitude"))


gram1<-variogram(res~1,data=dat2)  #filter(dat2,site=="ER")
plot(gram1)
gram2<-variogram(res~1,data=dat3)  #filter(dat2,site=="ER")
plot(gram2)

#fit.gram<-fit.variogram(gram, model = vgm(psill=4,model= "Exp",range=1,nugget=1.5))
#plot(gram,fit.gram)#correlation in residuals seems to diminish after 1 meter distance??



# temporal autocorrelation
dat2.2<-dat2%>%group_by(site,Year)%>%summarise(res=mean(res))%>%
  arrange(site,Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))

plot(dat2.2$res,dat2.2$lag2,xlab="t",ylab="t-1")




#plot residuals over time and lat
plot(pres_train$res~pres_train$latitude)

plot(surv_train$res~surv_train$latitude)

plot(surv_train$res~surv_train$Year)

stripchart(surv_train$res~surv_train$site,pch=1,method="jitter")

#not much pattern across years, but seem to differ between sites for survival and latitude for placement models


#solutions:
#model spatial or temporal relationships
#use lagged response variable as covariate
#mixed effects models
#include the variable in the model
#nest data in hierarchical structure 

# b) compare a structured model to one without structure using a hypothesis test like ANOVA (JUST FOR NEST SURVIVAL)

#glmer_surv<-glmer(paste0(deparse1(mod_list_surv[[3]]$formula),"+ (1|site)"),data=surv_train, family = binomial(link="logit"))
#block_surv<-glm(paste0(deparse1(mod_list_surv[[3]]$formula),"+ site"),data=surv_train, family = binomial(link="logit"))

#anova(mod_list_surv[[3]],glmer_surv,test="Chisq")

#anova(mod_list_surv[[3]],block_surv,test="LRT") #doesn't seem like adding site improves model fit (significantly far from max)

#add model with site as block for survival model list
#mod_list_surv[[4]]<-block_surv
#names(mod_list_surv)[[4]]<- "Site Block"
#model_selection[7,"response"]<-"success"
#model_selection[7,"mod_name"]<-"Low VIF, Site Block"
#model_selection[7,"mod_function"]<-deparse1(mod_list_surv[[4]]$formula)
#model_selection[7,"AIC"]<-mod_list_surv[[4]]$aic




## 5. Should we consider interactions?
#---------------------------------------------------
# Hypothesis: Does the continuity of high marsh matter for nesting or survival? - Yes, seems like more nests when placed in high marsh and completely surrounded by high marsh
ggplot(pres_train,aes(x=as.numeric(y),y=HIMARSH,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2),aes(col=as.factor(Highmarsh)))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placed?", y = "Proportion High Marsh in Buffer",color="Nest in High Marsh") + 
  theme_classic()
ggplot(surv_train,aes(x=as.numeric(y),y=HIMARSH,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2),aes(col=as.factor(Highmarsh)))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Survived?", y = "Proportion High Marsh in Buffer",color="Nest in High Marsh") + 
  theme_classic()

# Hypothesis: Does the quality/resiliency of high marsh matter for nesting or survival?

  # For nest presence? (pca, ent and cor textures)
    #of the locations in high marsh, do locations with nests have higher NDVI, etc?

ggplot(filter(pres_train,Highmarsh==1),aes(x=as.numeric(y),y=uvvr_mean,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "UVVR Mean") + 
  theme_classic()
ggplot(filter(pres_train,Highmarsh==1),aes(x=as.numeric(y),y=ndvi,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "NDVI") + 
  theme_classic()
ggplot(filter(pres_train,Highmarsh==1),aes(x=as.numeric(y),y=ent_txt,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "Entropy") + 
  theme_classic()
ggplot(filter(pres_train,Highmarsh==1),aes(x=as.numeric(y),y=cor_txt,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "Correlation") + 
  theme_classic()
    #pca only seems to have an impact
ggplot(filter(pres_train,Highmarsh==1),aes(x=as.numeric(y),y=pca,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "PC1") + 
  theme_classic()
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_pca_placement_interaction.png"), width = 7, height = 6, dpi = "retina")

    #add a model with this interaction
mod_list_pres[[4]]<-update(mod_list_pres[["High Marsh Presence"]],.~ ndvi + cor_txt + ent_txt + pca*Highmarsh + 
                              uvvr_mean)#+lat_scale 
names(mod_list_pres)[[4]]<-"High Marsh Interaction"
model_selection[8,"response"]<-"placement"
model_selection[8,"mod_name"]<-"High Marsh Interaction"
model_selection[8,"mod_function"]<-deparse1(mod_list_pres[[4]]$formula)
model_selection[8,"AIC"]<-mod_list_pres[[4]]$aic

  # for nest survival? (NDVI)
ggplot(filter(surv_train,Highmarsh==1),aes(x=as.numeric(y),y=ndvi,group=as.numeric(y)))+
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))+
  scale_x_continuous(breaks=c(0,1),labels = c("No","Yes"))+
  labs(x = "Nest Placement", y = "PC1") + 
  theme_classic()
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_ndvi_success_interaction.png"), width = 7, height = 6, dpi = "retina")


    #add a model with this interaction
mod_list_surv[[4]]<-update(mod_list_surv[[3]],.~ ndvi*Highmarsh + cor_txt + ent_txt + uvvr_mean)#+lat_scale
names(mod_list_surv)[[4]]<-"High Marsh Interaction"
model_selection[9,"response"]<-"success"
model_selection[9,"mod_name"]<-"High Marsh Interaction"
model_selection[9,"mod_function"]<-deparse1(mod_list_surv[[4]]$formula)
model_selection[9,"AIC"]<-mod_list_surv[[4]]$aic
