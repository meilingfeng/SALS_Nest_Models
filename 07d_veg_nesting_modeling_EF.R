library(tidyverse)
library(glmmTMB)
library(car)
library(bbmle)
library(DHARMa)
library(broom.mixed)
library(ggstats)
library(buildmer)
library(patchwork)
library(ggcorrplot)
library(factoextra)
library(terra)

## Set up data
#----------------------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#read in predictions at veg sites dataset and the veg survey data
dat_preds<-read.csv(paste0(path_out,"Final_outputs/Survey_Veg_Predictions/rapid_veg_BRT_predictions_CIs_15buff.csv"))%>%
  distinct(id,Year,.keep_all = T)
dat_veg<-read.csv(paste0(path_out,"Intermediate_outputs/Survey_Vegetation/processed_rapid_veg.csv"))%>%
  rename(Year=year)%>%
  distinct(id,Year,.keep_all = T)

dat<-left_join(dat_veg,dat_preds,by=c("id","Year"))%>%#change to left join
  mutate(id=as.factor(id),
         PatchID=as.factor(PatchID),
         ObserverInitials=as.factor(ObserverInitials),
         SHARPTide=as.factor(SHARPTide))

#filter out observations with missing values
dat<-dat%>%
  filter(!is.na(SHARPTide)|!is.na(DeadSnags))

#only select points in breeding range
veg_coords<-st_as_sf(dat,coords=c("Long","Lat"),crs="EPSG:4269")%>%
  st_transform("EPSG:26918")
long_min<-as.numeric(st_bbox(veg_coords)[1])
long_max<-as.numeric(st_bbox(veg_coords)[3])
#Focal species breeding range limits within the study region
range_limits<-read.csv(paste0(path_out,"Intermediate_outputs/ME_VA_range_limits_focal_spp.csv"))
lat_max<-range_limits[range_limits$species=="SALS","ymax"]
#list with 30m res tidal marsh layers
load(paste0(path_out,"/predictor_files_all_zones_30m.rds"))
lat_min<-as.numeric(ext(rast(file_list_all_zones[[8]][[1]]))[3])
veg_coords_clip<-st_crop(veg_coords,st_bbox(c(xmin=long_min,xmax=long_max,ymin=lat_min,ymax=lat_max),crs=st_crs(veg_coords)))

dat<-dat%>%
  filter(id%in%veg_coords_clip$id)

#center Latitude
dat$Latitude<-as.numeric(scale(dat$Latitude,center=T,scale=T))
dat<-dat%>%
  mutate(across(ends_with(c("pct","Snags")),~scale(.,center=F,scale=T)))



#calculate CI ranges and take ratio of minimum range:range to down weight obs with larger ranges (greater uncertainty) 
# (for weighting observations in regression)
dat<-dat%>%
  mutate(w_pres=u.ci_pres-l.ci_pres,
         w_surv=u.ci_surv-l.ci_surv,
         w_pres=min(w_pres)/w_pres,
         w_surv=min(w_surv)/w_surv)


# Check covariance among veg survey predictors
#-----------------------------------------------------------------------------
# Check for covariance among vegetation cover variables
# Run a principal components analysis on the vegetation cover variables to see if they are all needed
dat.cor<-dat%>%
  dplyr::select(ends_with("pct"))

dat.corr<-cor(dat.cor) 
ggcorrplot(dat.corr,ggtheme="theme_bw")
# high marsh species and categories correlated, same with low marsh, and invasive/phrag categories


#takes 10 PCs to explain 90% of variance. First PCs only capture 20% or less.
dat.pca.com<-dat%>%
  dplyr::select(-starts_with(c("alt","disti","gerardi","patens")))%>%
  dplyr::select(ends_with(c("pct")))
dat.pca.sp<-dat%>%
  dplyr::select(-starts_with(c("low","high")))%>%
  dplyr::select(ends_with(c("pct")))

dat.pca <- prcomp(dat.pca.com, scale = TRUE)
fviz_eig(dat.pca)


hist(dat$surv_pred)
ggplot(dat,aes(x=high_marsh_pct,y=pres_pred))+
  geom_point()+
  geom_smooth(method="gam")

#Check collinearity problems of predictor variables (use the broad community classes and individual species in sep models)
library(usdm)
vifcor(dat.cor[,-c(8,20)],th=0.6)
vifstep(dat.cor[,-c(1:5)],th=10)

# any transformations needed? how are the covariates distributed
#maybe log transform rarer categories (seem to be skewed) -all except high_marsh_pct
hist(dat$pool_panne_channel_pct)
#dat<-

## Test different stochastic distribution fits for nesting probability
#---------------------------------------------------------------------------
# weight observations by the width of their CI's. 
# account for multiple observations at the same marsh patch (individual observations are replicates across different locations and years within patch)

#Presence
####
# normal beta
pres_mod<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                    brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                    Latitude+I(Latitude^2)+
                    (1|SHARPTide)+
                    (1|PatchID),
                  data=dat,
                  family=beta_family(link = "logit"))#,
                  #weights=w_pres)


# zero inflation (intercept only)
zipres_mod<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                      brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                      Latitude+I(Latitude^2)+
                      (1|SHARPTide)+
                      (1|PatchID),
                    ziformula = ~1,
                    family=beta_family(link = "logit"),
                    data=dat)#,
                    #weights=w_pres)

# test gamma or log-normal (account for over dispersion- allow variance to increase with mean)
pres_mod_ln<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                       brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                       Latitude+ I(Latitude^2)+
                       (1|SHARPTide)+
                       (1|PatchID),
                      data=dat,
                      family=lognormal(link = "log"))#,
                      #weights=w_pres)

#diagnostics with simulated residuals
pres_simres1 <- simulateResiduals(pres_mod)
plot(pres_simres1)
pres_simres2 <- simulateResiduals(zipres_mod)
plot(pres_simres2)
pres_simres3 <- simulateResiduals(pres_mod_ln)
plot(pres_simres3)


# see if residuals vary with any predictors (should they be included, are they linear or quadratic?)
plotResiduals(pres_simres1,form=dat$Latitude)
plotResiduals(pres_simres1,form=dat$patens_pct)



# check for over dispersion (>1)
testDispersion(pres_simres3)

# check for zero inflation compared to that expected from the stochastic distribution we picked
testZeroInflation(pres_simres3)

#check for spatial autocorrelation
pres_simres1_noYear <- recalculateResiduals(pres_simres3, group = dat$PatchID,
                                            rotation = "estimated")
testSpatialAutocorrelation(pres_simres1_noYear,
                           x =  aggregate(dat$Long, list(dat$PatchID), mean)$x,
                           y = aggregate(dat$Lat, list(dat$PatchID), mean)$x)


#add fitted values to dataframe
fit_obs<-data.frame(pres_mod=predict(pres_mod,newdata=dat,type="response"),
                    zipres_mod=predict(zipres_mod,newdata=dat,type="response"),
                    pres_mod_ln=predict(pres_mod_ln,newdata=dat,type="response"),
                    obs=dat$pres_pred)
#compare with original distribution
hist(fit_obs$obs)
hist(fit_obs$pres_mod)
#how does the estimated mean look?
plogis(pres_mod_ln$fit$par[1])


#it looks like a log-normal distribution is the best stochastic model


# Nest Fate
####
# normal beta
surv_mod<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                    brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                    Latitude+I(Latitude^2)+
                    (1|SHARPTide)+
                    (1|PatchID),
                  data=dat,
                  family=beta_family(link = "logit"))#,
#weights=w_pres)


# zero inflation (intercept only)
zisurv_mod<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                      brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                      Latitude+I(Latitude^2)+
                      (1|SHARPTide)+
                      (1|PatchID),
                    ziformula = ~1,
                    family=beta_family(link = "logit"),
                    data=dat)#,
#weights=w_pres)

# test gamma or log-normal (account for over dispersion- allow variance to increase with mean)
surv_mod_ln<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                       brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                       Latitude+ I(Latitude^2)+
                       (1|SHARPTide)+
                       (1|PatchID),
                     data=dat,
                     family=lognormal(link = "log"))#,
#weights=w_pres)

#diagnostics with simulated residuals
surv_simres1 <- simulateResiduals(surv_mod)
plot(surv_simres1)
surv_simres2 <- simulateResiduals(zisurv_mod)
plot(surv_simres2)
surv_simres3 <- simulateResiduals(surv_mod_ln)
plot(surv_simres3)


# check for over dispersion (>1)
testDispersion(surv_simres1)

# check for zero inflation compared to that expected from the stochastic distribution we picked
testZeroInflation(surv_simres3)

#check for spatial autocorrelation
surv_simres1_noYear <- recalculateResiduals(surv_simres3, group = dat$PatchID,
                                            rotation = "estimated")
testSpatialAutocorrelation(surv_simres1_noYear,
                           x =  aggregate(dat$Long, list(dat$PatchID), mean)$x,
                           y = aggregate(dat$Lat, list(dat$PatchID), mean)$x)


#add fitted values to dataframe
fit_obs<-data.frame(surv_mod=predict(surv_mod,newdata=dat,type="response"),
                    zisurv_mod=predict(zisurv_mod,newdata=dat,type="response"),
                    surv_mod_ln=predict(surv_mod_ln,newdata=dat,type="response"),
                    obs=dat$surv_pred)
#compare with original distribution
hist(fit_obs$obs)
hist(fit_obs$surv_mod_ln)
#how does the estimated mean look?
plogis(surv_mod_ln$fit$par[1])


#list to contain all the models
pres_list<-list(pres_mod,zipres_mod,pres_mod_ln)
surv_list<-list(surv_mod,zisurv_mod,surv_mod_ln)


# Use AIC to determine which stochastic distribution has the best fit
AICtab(pres_list)
AICtab(surv_list)

#log-normal appears on top for both presence and fate models

# Compare different veg models:
#----------------------------------------------------------------------------------

# versions of these models that use quadratic terms
pres_list<-list()
surv_list<-list()

# model that looks at percent cover of broader marsh veg communities 
pres_list[[1]]<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at percent cover of broader veg communities + broader land cover classes
pres_list[[2]]<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                                       brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                                       Latitude+ I(Latitude^2)+
                                       (1|SHARPTide)+
                                       (1|PatchID),
                                     data=dat,
                                     family=lognormal(link = "log"))


# model that looks at percent cover of individual species known to influence nesting
pres_list[[3]]<-glmmTMB(pres_pred~ alt_tall_pct+alt_short_pct+distichlis_pct+gerardii_pct+patens_pct+phrag_pct+
                                       Latitude+ I(Latitude^2)+
                                       (1|SHARPTide)+
                                       (1|PatchID),
                                     data=dat,
                                     family=lognormal(link = "log"))

# model that looks at percent cover of individual species known to influence nesting + broad land cover+community classes
pres_list[[4]]<-glmmTMB(pres_pred~ alt_tall_pct+alt_short_pct+distichlis_pct+gerardii_pct+patens_pct+phrag_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at broad communities + snags + broad cover classes
pres_list[[5]]<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          DeadSnags+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at broad communities + less common species + broad cover classes
pres_list[[6]]<-glmmTMB(pres_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          Schoenoplectus_pct+Morella_pct+Salicornia_pct+Scirpus_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))




# model that looks at percent cover of broader marsh veg communities 
surv_list[[1]]<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at percent cover of broader veg communities + broader land cover classes
surv_list[[2]]<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))


# model that looks at percent cover of individual species known to influence nesting
surv_list[[3]]<-glmmTMB(surv_pred~ alt_tall_pct+alt_short_pct+distichlis_pct+gerardii_pct+patens_pct+phrag_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at percent cover of individual species known to influence nesting + broad land cover+community classes
surv_list[[4]]<-glmmTMB(surv_pred~ alt_tall_pct+alt_short_pct+distichlis_pct+gerardii_pct+patens_pct+phrag_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at broad communities + snags + broad cover classes
surv_list[[5]]<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          DeadSnags+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID),
                        data=dat,
                        family=lognormal(link = "log"))

# model that looks at broad communities + less common species + broad cover classes
surv_list[[6]]<-glmmTMB(surv_pred~ high_marsh_pct+low_marsh_pct+invasives_pct+
                          brackish_border_pct+saltmarsh_border_pct+upland_pct+trees_pct+water_pct+pool_panne_channel_pct+
                          Schoenoplectus_pct+Morella_pct+Salicornia_pct+Scirpus_pct+
                          Latitude+ I(Latitude^2)+
                          (1|SHARPTide)+
                          (1|PatchID/SHARPTide),
                        data=dat,
                        family=lognormal(link = "log"))

## Model Comparison
#---------------------------------------------------------------------
AICtab(pres_list)
AICtab(surv_list)
# tells us whether measurable vegetation characteristics (ranging from simple to more detailed) describe signatures of nest sites

#broad communities+land cover most parsimonious for presence
#breaking down communities into individual species known to influence nest selection was important for nest fate
#snags and less common species did not improve model fit.


#use a gam() with betr family?
# that way we can find an optimal amount of veg
library(mgcv)



## Model Inference
#---------------------------------------------------------------------
# show which veg variables describe characteristics of nest sites (non-zero coeffs)
# then show the relationships between those significant veg variables and nest sites

summary(pres_list[[2]])
summary(surv_list[[5]])
#log transformation= coefficients are multiplicative on the data scale (they are a percent of the mean nesting probability)


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
  xlab("Percent Change in Probability of Nest Presence")+
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
