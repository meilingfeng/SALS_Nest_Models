library(tidyverse)
library(glmmTMB)
library(car)
library(bbmle)
library(DHARMa)
library(broom.mixed)
library(ggstats)
library(buildmer)
library(patchwork)


dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#read in veg-nest prediction dataset
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_Nest_Predictions_at_Rapid_Survey_Points.csv"))



#center/scale Latitude
dat[,5]=scale(dat[,5])
#convert the veg percent cover to proportions to match scale
dat[,c(19:32)]<-dat[,c(19:32)]/100



## Test different error distribution fits
#---------------------------------------------------------------------------

# presence
####
# normal beta
pres_mod<-glmmTMB(presence~alt_tall_pct+distichlis_pct+gerardii_pct+
                       patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                       high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                       upland_pct+trees_pct+water_pct+Lat + I(Lat^2),data=dat,family=ordbeta(link = "logit"))
# zero inflation (intercept only)
zipres_mod<-glmmTMB(presence~alt_tall_pct+distichlis_pct+gerardii_pct+
                      patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                      high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                      upland_pct+trees_pct+water_pct+Lat + I(Lat^2),data=dat,
                    ziformula= ~1,  family=ordbeta(link = "logit"))

# zero inflation using latitude
zipres_mod_lat<-glmmTMB(presence~alt_tall_pct+distichlis_pct+gerardii_pct+
                       patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                       high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                       upland_pct+trees_pct+water_pct,data=dat,
                     ziformula= ~Lat + I(Lat^2),  family=ordbeta(link = "logit") )
#try simplifying the deterministic equation
#stepmod_pres_pct<-buildglmmTMB(presence~alt_tall_pct+distichlis_pct+gerardii_pct+
#                                 patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
#                                 high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
#                                 upland_pct+trees_pct+water_pct+Lat+I(Lat^2),data=dat,family=ordbeta(link = "logit"))
#summary(stepmod_pres_pct)# all terms significant based on change in log likelihood, AIC, BIC, and explained deviance - forward and backward effect selection methods

#diagnostic plots with simulated residuals
pres_simres1 <- simulateResiduals(pres_mod)
plot(pres_simres1)
pres_simres2 <- simulateResiduals(zipres_mod)
plot(pres_simres2)
pres_simres3 <- simulateResiduals(zipres_mod_lat)
plot(pres_simres3)


#add fitted values to dataframe
fit_obs<-data.frame(pres_mod=predict(pres_mod,newdata=dat,type="response"),
                    zipres_mod=predict(zipres_mod,newdata=dat,type="response"),
                    zipres_mod_lat=predict(zipres_mod_lat,newdata=dat,type="response"),
                    obs=dat$presence)
#compare with original distribution
hist(dat$presence)
hist(fit_obs$pres_mod)
#how does the estimated mean look?
plogis(pres_mod$fit$par[1])



# survival
####
# normal beta
surv_mod<-glmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
                    patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                    high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                    upland_pct+trees_pct+water_pct+Lat + I(Lat^2),data=dat,family=ordbeta(link = "logit"))


# zero inflation (intercept only)
zisurv_mod<-glmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
                      patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                      high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                      upland_pct+trees_pct+water_pct+Lat + I(Lat^2),data=dat,
                    ziformula= ~1,  family=ordbeta(link = "logit"))

# zero inflation using latitude
zisurv_mod_lat<-glmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
                          patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                          high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                          upland_pct+trees_pct+water_pct,data=dat,
                        ziformula= ~Lat + I(Lat^2),  family=ordbeta(link = "logit") )

#can we drop any variables? All seem to be significant.
#stepmod_surv_pct<-buildglmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
#                                 patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
#                                 high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
#                                 upland_pct+trees_pct+water_pct+Lat+I(Lat^2),data=dat,family=ordbeta(link = "logit"))
#summary(stepmod_survival_pct)
#record the simplified model
#  simplified_surv_pct_Lat2<-betareg(survival ~ alt_tall_pct+ distichlis_pct+gerardii_pct +patens_pct + alt_short_pct+ phrag_pct +low_marsh_pct+ high_marsh_pct + brackish_border_pct  + upland_pct + water_pct  + Lat+I(Lat^2), data = rapid_dat10.model )
#record the simplified model with Lat^2 removed for further comparison
#simplified_surv_pct<-stepmod_survival_pct@model



#diagnostic plots with simulated residuals
surv_simres1 <- simulateResiduals(surv_mod)
plot(surv_simres1)
surv_simres2 <- simulateResiduals(zisurv_mod)
plot(surv_simres2)
surv_simres3 <- simulateResiduals(zisurv_mod_lat)
plot(surv_simres3)


#add fitted values to dataframe
fit_obs<-fit_obs%>%
  mutate(surv_mod=predict(surv_mod,newdata=dat,type="response"),
         zisurv_mod=predict(zisurv_mod2,newdata=dat,type="response"),
         zisurv_mod_lat=predict(zisurv_mod_lat,newdata=dat,type="response"))
#compare with original distribution
hist(dat$survival)
hist(fit_obs$surv_mod)
#how does the estimated mean look?
plogis(surv_mod2$fit$par[1])



#list to contain all the models
pres_list<-list(pres_mod,zipres_mod,zipres_mod_lat)
surv_list<-list(surv_mod,zisurv_mod,zisurv_mod_lat)



## Model Comparison
#---------------------------------------------------------------------
AICtab(pres_list)
AICtab(surv_list)


## Model Inference
#---------------------------------------------------------------------

summary(surv_mod)


#format the coefficients and confidence intervals and plot
#for presence
t.pres <- broom.mixed::tidy(pres_list[[2]], conf.int = TRUE)
t.pres <- transform(t.pres,
                   term=sprintf("%s.%s", component, term))%>%
  filter(effect!="ran_pars")
t.pres$term<-str_replace(t.pres$term,"\\(","")
t.pres$term<-str_replace(t.pres$term,"\\)","")

t.pres2<-t.pres%>%
  mutate(estimate=exp(estimate),
            conf.low=exp(conf.low),
            conf.high=exp(conf.high),
         term=as.factor(substring(term,6,nchar(term))))%>%
  filter(component!="zi"&term!="Intercept")
t.pres2$term<-factor(t.pres2$term, levels=c('patens_pct', 'gerardii_pct', 'distichlis_pct', 'alt_short_pct', 'alt_tall_pct',
                                      'high_marsh_pct','brackish_border_pct','saltmarsh_border_pct','upland_pct',
                                      'phrag_pct','low_marsh_pct','trees_pct','water_pct',
                                      'Lat','ILat^2'),
                     labels=c("S. patens","J. gerardii","D. spicata","Short S. alterniflora","Tall S. alterniflora",
                              "High Marsh","Brackish Border","Saltmarsh Border","Upland",
                              "Phragmites","Low Marsh","Tree Cover","Water",
                              "Latitude","Latitude^2"))

pres_plot<-ggplot(t.pres2, aes(x=estimate, y=term, color=term)) + 
  geom_vline(xintercept = 1)+
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  geom_vline(xintercept=c(0,2),color="red",linetype="longdash")+
  geom_point()+
  scale_x_continuous(breaks = c(0:2,4,6,8,10,12,14))+
  theme_bw(base_size = 12)+
  theme(legend.position = "none")+
  xlab("Odds Ratio of Nest Presence")+
  ylab("Vegetation Type (Proportion Cover)")




# for survival
t.ind <- broom.mixed::tidy(surv_list[[1]], conf.int = TRUE)
t.ind <- transform(t.ind,
                   term=sprintf("%s.%s", component, term))%>%
  filter(effect!="ran_pars")
t.ind$term<-str_replace(t.ind$term,"\\(","")
t.ind$term<-str_replace(t.ind$term,"\\)","")

t.ind2<-t.ind%>%
  mutate(estimate=exp(estimate),
         conf.low=exp(conf.low),
         conf.high=exp(conf.high),
         term=as.factor(substring(term,6,nchar(term))))%>%
  filter(component!="zi"&term!="Intercept")
t.ind2$term<-factor(t.ind2$term, levels=c('patens_pct', 'gerardii_pct', 'distichlis_pct', 'alt_short_pct', 'alt_tall_pct',
                                          'high_marsh_pct','brackish_border_pct','saltmarsh_border_pct','upland_pct',
                                          'phrag_pct','low_marsh_pct','trees_pct','water_pct',
                                          'Lat','ILat^2'),
                    labels=c("S. patens","J. gerardii","D. spicata","Short S. alterniflora","Tall S. alterniflora",
                             "High Marsh","Brackish Border","Saltmarsh Border","Upland",
                             "Phragmites","Low Marsh","Tree Cover","Water",
                             "Latitude","Latitude^2"))

surv_plot<-ggplot(t.ind2, aes(x=estimate, y=term, color=term)) + 
  geom_vline(xintercept = 1)+
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  geom_vline(xintercept=c(0,2),color="red",linetype="longdash")+
  geom_point()+
  theme_bw(base_size = 12)+
  theme(legend.position = "none")+
  xlab("Odds Ratio of Nesting Success")+
  ylab("")


pres_plot+surv_plot+plot_annotation(tag_levels = "A")
ggsave(filename=paste0(path_out,"Final_outputs/Model_Results/veg_coefs_pres_surv.png"), width = 11, height = 4, dpi = "retina")

# ANOVA
Anova(pres_list[[2]],type="III")
Anova(surv_list[[1]],type="III")
