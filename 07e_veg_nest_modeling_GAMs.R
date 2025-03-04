#use a gam() with betr family?
# that way we can find an optimal amount of veg

library(mgcv)
library(tidyverse)
library(bbmle)
library(sf)
library(terra)
## Set up data
#----------------------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


#make predictions at each plot location, set null variables to standard values that max presence and success
dat_preds<-read.csv(paste0(path_out,"Final_outputs/Survey_Veg_Predictions/rapid_veg_BRT_predictions_CIs_15buff_single_plots.csv"))

# take the mean of veg characteristics at each plot
dat_veg<-read.csv(paste0(path_out,"Intermediate_outputs/Survey_Vegetation/processed_rapid_veg.csv"))%>%
  group_by(id)%>%
  summarise(across(ends_with(c("Snags","pct","PatchID","Lat","Long")),~mean(.,na.rm=T)))

dat<-left_join(dat_veg,dat_preds,by=c("id"))%>%#change to left join
  mutate(id=as.factor(id),
         PatchID=as.factor(PatchID))

#filter out observations with missing values
#dat<-dat%>%
#  filter(!is.na(DeadSnags))

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

#cap the percent cover of individual classes at 100
dat<-dat%>%
  mutate(across(ends_with("pct"), ~ifelse(.>100,100,.)))

#center/scale Latitude, scale percents
#dat$Latitude<-as.numeric(scale(dat$Latitude,center=T,scale=T))
#dat<-dat%>%
#  mutate(across(ends_with(c("pct","Snags")),~as.numeric(scale(.,center=F,scale=T))))



#calculate CI ranges and take ratio of minimum range:range to down weight obs with larger ranges (greater uncertainty) 
# (for weighting observations in regression)
dat<-dat%>%
  mutate(w_pres=u.ci_pres-l.ci_pres,
         w_surv=u.ci_surv-l.ci_surv,
         w_pres=min(w_pres)/w_pres,
         w_surv=min(w_surv)/w_surv)



#Test different stochastic distributions and components to account for data structure
#--------------------------------------------------------------------------------------
# A) consider variation among marsh patches?
#  don't expect relationship between ground characteristics and habitat signatures to vary among marshes
#  most patches also only have 1 survey plot
patch_samples<-count(dat,PatchID)
hist(patch_samples$n)
#almost half of marsh patches only have 1 survey plot
nrow(patch_samples[patch_samples$n>1,])/nrow(patch_samples)




# B) Start with a simple model to test stochastic distributions
# beta best meets the values of our response (0-1)
gam.plain<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr"),
            data=dat, weight=w_surv/mean(w_surv),
            family = betar(link="logit")) 
summary(gam.plain)
par(mfrow=c(2,2)) 
plot(gam.plain)

par(mfrow=c(2,2)) 
gam.check(gam.plain)

# try a different distribution? Gamma or lognormal?
gam.plain.gm<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr"),
                 data=dat, weight=w_surv/mean(w_surv),
                 family = Gamma(link="log")) 
summary(gam.plain.gm)
par(mfrow=c(2,2)) 
plot(gam.plain.gm)

par(mfrow=c(2,2)) 
gam.check(gam.plain.gm)


#Try plotting the predictions to see which distribution's predictions match the data better
newdat<-data.frame(high_marsh_pct=seq(0,100,by=5))
#extract predictions
pred_b<-predict(gam.plain,newdata=newdat,type="response",se.fit=T)
pred_g<-predict(gam.plain.gm,newdata=newdat,type="response",se.fit=T)

newdat$b.fit<-pred_b$fit
newdat$b.se<-pred_b$se.fit
newdat$g.fit<-pred_g$fit
newdat$g.se<-pred_g$se.fit


#plot
ggplot(dat,aes(x=high_marsh_pct,y=surv_pred))+
  geom_ribbon(data=newdat,aes(y=b.fit,
                              ymin=b.fit-1.96*b.se,
                              ymax=b.fit+1.96*b.se,
                              alpha=0.3))+
  geom_line(data=newdat,aes(y=b.fit))+
  geom_point()+
  theme_bw()


ggplot(dat,aes(x=high_marsh_pct,y=surv_pred))+
  geom_ribbon(data=newdat,aes(y=g.fit,
                              ymin=g.fit-1.96*g.se,
                              ymax=g.fit+1.96*g.se,
                              alpha=0.3))+
  geom_line(data=newdat,aes(y=g.fit))+
  geom_point()+
  theme_bw()

#gamma looks like it better accounts for the uncertainty at higher proportions of veg cover




# C) consider spatial autocorr?
# plots within a marsh (spatially close) may share veg characteristics and RS signatures
gam.spat<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr")+
                   s(Long,Lat,k=50),
                 data=dat, weight=w_surv/mean(w_surv),
              family = betar(link="logit")) 

gam.spat.gm<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr")+
                 s(Long,Lat,k=50),
               data=dat, weight=w_surv/mean(w_surv),
               family = Gamma(link="log")) 

summary(gam.spat.gm)
par(mfrow=c(2,2)) 
plot(gam.spat.gm)

par(mfrow=c(2,2)) 
gam.check(gam.spat)
dev.off()

par(mfrow=c(2,2)) 
gam.check(gam.spat.gm)
dev.off()


#Try plotting the predictions to see which matches the data better
newdat<-dat%>%
  dplyr::select(Lat,Long,high_marsh_pct)
#extract predictions
pred_b<-predict(gam.spat,newdata=newdat,type="response",se.fit=T)
pred_g<-predict(gam.spat.gm,newdata=newdat,type="response",se.fit=T)

newdat$b.fit<-pred_b$fit
newdat$b.se<-pred_b$se.fit
newdat$g.fit<-pred_g$fit
newdat$g.se<-pred_g$se.fit


#plot
ggplot(dat,aes(x=high_marsh_pct,y=surv_pred))+
  geom_ribbon(data=newdat,aes(y=b.fit,
                              ymin=b.fit-1.96*b.se,
                              ymax=b.fit+1.96*b.se,
                              alpha=0.3))+
  geom_line(data=newdat,aes(y=b.fit))+
  geom_point()+
  theme_bw()


ggplot(dat,aes(x=high_marsh_pct,y=surv_pred))+
  geom_ribbon(data=newdat,aes(y=g.fit,
                              ymin=g.fit-1.96*g.se,
                              ymax=g.fit+1.96*g.se,
                              alpha=0.3))+
  geom_line(data=newdat,aes(y=g.fit))+
  geom_point()+
  theme_bw()


# the gamma model is predicting values above 1, which the beta distribution accounts for

#does adding the spatial component improve model fit? A gamma distribution also seems to fit better.
AICtab(gam.plain,gam.spat,gam.plain.gm,gam.spat.gm)
#yes



# Test different deterministic models
#----------------------------------------------------------------------------
# Vegetation community model (with and without broader land cover types)
# broad tidal marsh vegetation communities + additional land cover types
com<-gam(pres_pred~s(high_marsh_pct,k=10,bs="cr")+s(low_marsh_pct,k=10,bs="cr")+
                 s(invasives_pct,k=10,bs="cr")+s(saltmarsh_border_pct,k=10,bs="cr")+
                 s(brackish_border_pct,k=10,bs="cr")+
                 s(Long,Lat,k=50),
               data=dat, weight=w_pres/mean(w_pres),
         family = Gamma(link="log"),
         method = "REML") 
com.lc<-gam(pres_pred~s(high_marsh_pct,k=10,bs="cr")+s(low_marsh_pct,k=10,bs="cr")+
                 s(invasives_pct,k=10,bs="cr")+s(saltmarsh_border_pct,k=10,bs="cr")+
                 s(brackish_border_pct,k=10,bs="cr")+s(trees_pct,k=10,bs="cr")+
                 s(upland_pct,k=10,bs="cr")+s(pool_panne_channel_pct,k=10,bs="cr")+
                 s(water_pct,k=10,bs="cr")+
                 s(Long,Lat,k=50),
               data=dat, weight=w_pres/mean(w_pres),
            family = Gamma(link="log"),
            method = "REML") 
# Species specific model (with and without broader land cover types)
# specific species within each vegetation community known to influence nest site selection + additional land cover types
sp<-gam(pres_pred~s(patens_pct,k=10,bs="cr")+s(alt_tall_pct,k=10,bs="cr")+s(alt_short_pct,k=10,bs="cr")+
                s(gerardii_pct,k=10,bs="cr")+s(distichlis_pct,k=10,bs="cr")+s(phrag_pct,k=10,bs="cr")+
                s(Baccharis_pct,k=10,bs="cr")+s(Iva_pct,k=10,bs="cr")+
                s(Long,Lat,k=50),
              data=dat, weight=w_pres/mean(w_pres),
        family = Gamma(link="log"),
        method = "REML") 
sp.lc<-gam(pres_pred~s(patens_pct,k=10,bs="cr")+s(alt_tall_pct,k=10,bs="cr")+s(alt_short_pct,k=10,bs="cr")+
             s(gerardii_pct,k=10,bs="cr")+s(distichlis_pct,k=10,bs="cr")+s(phrag_pct,k=10,bs="cr")+
             s(Baccharis_pct,k=10,bs="cr")+s(Iva_pct,k=10,bs="cr")+
             s(trees_pct,k=10,bs="cr")+s(upland_pct,k=10,bs="cr")+s(pool_panne_channel_pct,k=10,bs="cr")+
             s(water_pct,k=10,bs="cr")+
             s(Long,Lat,k=50),
            data=dat, weight=w_pres/mean(w_pres),
           family = Gamma(link="log"),
           method = "REML") 





AICtab(com,com.lc,sp,sp.lc)


gam.check(sp.lc)
summary(sp.lc)


# Repeat for nest fate models
# Vegetation community model (with and without broader land cover types)
# broad tidal marsh vegetation communities + additional land cover types
com.fate<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr")+s(low_marsh_pct,k=10,bs="cr")+
           s(invasives_pct,k=10,bs="cr")+s(saltmarsh_border_pct,k=10,bs="cr")+
           s(brackish_border_pct,k=10,bs="cr")+
           s(Long,Lat,k=50),
         data=dat, weight=w_surv/mean(w_surv),
         family = Gamma(link="log"),
         method = "REML") 
com.lc.fate<-gam(surv_pred~s(high_marsh_pct,k=10,bs="cr")+s(low_marsh_pct,k=10,bs="cr")+
              s(invasives_pct,k=10,bs="cr")+s(saltmarsh_border_pct,k=10,bs="cr")+
              s(brackish_border_pct,k=10,bs="cr")+s(trees_pct,k=10,bs="cr")+
              s(upland_pct,k=10,bs="cr")+s(pool_panne_channel_pct,k=10,bs="cr")+
              s(water_pct,k=10,bs="cr")+
              s(Long,Lat,k=50),
            data=dat, weight=w_surv/mean(w_surv),
            family = Gamma(link="log"),
            method = "REML") 
# Species specific model (with and without broader land cover types)
# specific species within each vegetation community known to influence nest site selection + additional land cover types
sp.fate<-gam(surv_pred~s(patens_pct,k=10,bs="cr")+s(alt_tall_pct,k=10,bs="cr")+s(alt_short_pct,k=10,bs="cr")+
          s(gerardii_pct,k=10,bs="cr")+s(distichlis_pct,k=10,bs="cr")+s(phrag_pct,k=10,bs="cr")+
            s(Baccharis_pct,k=10,bs="cr")+s(Iva_pct,k=10,bs="cr")+
          s(Long,Lat,k=50),
        data=dat, weight=w_surv/mean(w_surv),
        family = Gamma(link="log"),
        method = "REML") 
sp.lc.fate<-gam(surv_pred~s(patens_pct,k=10,bs="cr")+s(alt_tall_pct,k=10,bs="cr")+s(alt_short_pct,k=10,bs="cr")+
             s(gerardii_pct,k=10,bs="cr")+s(distichlis_pct,k=10,bs="cr")+s(phrag_pct,k=10,bs="cr")+
               s(Baccharis_pct,k=10,bs="cr")+s(Iva_pct,k=10,bs="cr")+
             s(trees_pct,k=10,bs="cr")+s(upland_pct,k=10,bs="cr")+s(pool_panne_channel_pct,k=10,bs="cr")+
             s(water_pct,k=10,bs="cr")+
             s(Long,Lat,k=50),
           data=dat, weight=w_surv/mean(w_surv),
           family = Gamma(link="log"),
           method = "REML") 




AICtab(com.fate,com.lc.fate,sp.fate,sp.lc.fate)


gam.check(com.fate)
summary(com.fate)
# Plot the relationships between each veg covariate and nesting signatures
#------------------------------------------------------------------------------------
par(mfrow=c(2,2)) 
plot(sp)

par(mfrow=c(2,2)) 
plot(sp.lat.fate)

library(gratia)

#get mean response from model
coefs_pres<-sp.lc$coefficients
inter_pres<-as.numeric(coefs_pres[1])
#generate custom plot with smooth outputs
sm_pres <- smooth_estimates(sp.lc)%>%
  pivot_longer(cols=ends_with("pct"),names_to = "var",values_to="val")%>%
  mutate(var=case_when(
    var=="patens_pct"~"S. patens",
    var=="alt_tall_pct"~"S. alterniflora (Tall)",
    var=="alt_short_pct"~"S. alterniflora (Short)",
    var=="gerardii_pct"~"J. gerardii",
    var=="distichlis_pct"~"D. spicata",
    var=="phrag_pct"~"P. australis",
    var=="Baccharis_pct"~"B. halimifolia",
    var=="Iva_pct"~"Iva frutiscens",
    var=="trees_pct"~"Trees",
    var=="upland_pct"~"Upland",
    var=="pool_panne_channel_pct"~"Pools, pannes, channels",
    var=="water_pct"~"Open water"),
  u.l=ifelse(exp(inter_pres+(.estimate+(.se*1.95)))>1,1,exp(inter_pres+(.estimate+(.se*1.95)))),
  l.l=ifelse(exp(inter_pres+(.estimate-(.se*1.95)))>1,1,exp(inter_pres+(.estimate-(.se*1.95)))),
  .estimate=ifelse(exp(inter_pres+.estimate)>1,1,exp(inter_pres+.estimate))
  )

sm_pres %>%
  #add_confint()%>%
  ggplot(aes(y = .estimate, x = val)) +
  geom_ribbon(aes(ymin = l.l, ymax = u.l),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(colour = "forestgreen", linewidth = 1.5) +
  labs(
    y = "Partial effect",
    x="Percent Cover")+
  facet_wrap(~var,scales = "free")+
  theme_classic(base_size = 12)
ggsave(filename = paste0(path_out,"Final_outputs/Model_Results/veg_analysis_gams_pres.jpeg"),
       width = 9,height=6,dpi = "retina")

#repeat for fates
coefs_surv<-com.fate$coefficients
inter_surv<-as.numeric(coefs_surv[1])
#generate custom plot with smooth outputs
sm_surv <- smooth_estimates(com.fate)%>%
  pivot_longer(cols=ends_with("pct"),names_to = "var",values_to="val")%>%
  mutate(var=case_when(
    var=="high_marsh_pct"~"High marsh",
    var=="low_marsh_pct"~"Low marsh",
    var=="invasives_pct"~"Invasives",
    var=="saltmarsh_border_pct"~"Saltmarsh border",
    var=="brackish_border_pct"~"Brackish terrestrial border"),
    u.l=ifelse(exp(inter_surv+(.estimate+(.se*1.95)))>1,1,exp(inter_surv+(.estimate+(.se*1.95)))),
    l.l=ifelse(exp(inter_surv+(.estimate-(.se*1.95)))>1,1,exp(inter_surv+(.estimate-(.se*1.95)))),
    .estimate=ifelse(exp(inter_surv+.estimate)>1,1,exp(inter_surv+.estimate))
  )
unique(sm_surv$var)

sm_surv %>%
  #add_confint()%>%
  ggplot(aes(y = .estimate, x = val)) +
  geom_ribbon(aes(ymin = l.l, ymax = u.l),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(colour = "forestgreen", linewidth = 1.5) +
  labs(
    y = "Partial effect",
    x="Percent Cover")+
  facet_wrap(~var,scales = "free")+
  theme_classic(base_size = 12)+
  ylim(0.4,0.7)
ggsave(filename = paste0(path_out,"Final_outputs/Model_Results/veg_analysis_gams_surv.jpeg"),
       width = 8,height=5,dpi = "retina")
