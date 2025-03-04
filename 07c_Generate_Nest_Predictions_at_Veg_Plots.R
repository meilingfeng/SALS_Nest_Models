library(tidyverse)
library(terra)
library(sf)
library(dismo)

########
# generate nest predictions at SHARP survey veg sites
#############

out_list<-list()

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


#make predictions at each plot location, set null variables to standard values that max presence and success
# load vegetation survey data with covariates
veg_data<-read.csv(paste0(path_out,"Intermediate_outputs/Survey_Vegetation/veg_buff_RS_covariates_15buff.csv"))%>%
  dplyr::select(id,Latitude=Lat,Day.of.the.Year=doy,Mean.UVVR=uvvr_mean,NDVI=ndvi,Surface.Brightness=pca,Days.Since.Spring.Tide=time_since_tide,Year,
                Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH,Tidal.Restriction=tideres,Change.in.UVVR=uvvr_diff,Elevation=elevation)%>%
  #keep only one record of each location
  distinct(id,.keep_all = T)%>%
  #set null variables to max presence/fate based on PDP plots from BRTs 
  # (Success increases later in the season and decreases over years. Days since tide was set to 1 in prior step.)
  mutate(Day.of.the.Year=220,
         Year=2011)

#don't need to filter for complete observations, BRTs can predict with missing data.



### Predict
#---------------------------------------------------
# get nest predictions and uncertainty at each site

#get CI for predictions by boostrapping the nest data (training and refitting a model 1000 times)

#Load and tidy nest data
  all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 

  speciesnames<-c("SALS")

  pres_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_pres_dat.csv"))
  surv_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/SALS_nest_fate_dat.csv"))
  
  #number of bootstrapping iterations
  k<-100 #This leaves about 50 obs out for each pres iteration, and 20 observations for nest surv
  set.seed(123)
  pres_dat$group<-kfold(pres_dat,k)
  surv_dat$group<-kfold(surv_dat,k)
  count(pres_dat,group)
  count(surv_dat,group)
  # select all predictors and response
  pres_dat2<-pres_dat%>%dplyr::select("id","y","group",,"latitude","doy",all_of(all_terms))%>%
    rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
           Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH, Latitude=latitude,Day.of.the.Year=doy)
  surv_dat2<-surv_dat%>%dplyr::select("id","y","group","latitude","doy","time_since_tide","Year",,all_of(all_terms))%>%
    rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
           Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH, Latitude=latitude,Day.of.the.Year=doy, Days.Since.Spring.Tide=time_since_tide)
  

#Create objects to hold predictions for each iteration
    brt_predict_pres<-list()
    brt_predict_surv<-list()
    
  
#Bootstrap models with subsets of nest data and predict to veg sites    
for(i in 1:k){
      #use subset observations to train models
      pres_dat_train<-pres_dat2[pres_dat2$group!=i,]%>%dplyr::select(-group,-id)
      surv_dat_train<-surv_dat2[surv_dat2$group!=i,]%>%dplyr::select(-group,-id)
      
      
      # 1. train models
      #---------------------------------------------
      set.seed(123)
      #First for nest presence
      # typically want at least 1000 trees, so decrease learning rate if you get less
      brt_pres<-gbm.step(data=pres_dat_train, gbm.x = 2:length(pres_dat_train), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.05)# default bag fraction is 0.75
      #start at lr of 0.05, decrease until optimal trees hits 1000
      if(brt_pres$gbm.call$best.trees<1000){
        brt_pres<-gbm.step(data=pres_dat_train, gbm.x = 2:length(pres_dat_train), gbm.y=1, 
                           family = "bernoulli",
                           tree.complexity = 5,
                           learning.rate=0.01)
      }
      if(brt_pres$gbm.call$best.trees<1000){
        brt_pres<-gbm.step(data=pres_dat_train, gbm.x = 2:length(pres_dat_train), gbm.y=1, 
                           family = "bernoulli",
                           tree.complexity = 5,
                           learning.rate=0.005)
      }
      if(brt_pres$gbm.call$best.trees<1000){
        brt_pres<-gbm.step(data=pres_dat_train, gbm.x = 2:length(pres_dat_train), gbm.y=1, 
                           family = "bernoulli",
                           tree.complexity = 5,
                           learning.rate=0.001)
      }
      
      
      # Then nest survival
      brt_surv<-gbm.step(data=surv_dat_train, gbm.x = 2:length(surv_dat_train), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.01)
      
      # start at lr of 0.05, decrease until optimal trees hits 1000
      if(brt_surv$gbm.call$best.trees<1000){
        brt_surv<-gbm.step(data=surv_dat_train, gbm.x = 2:length(surv_dat_train), gbm.y=1, 
                          family = "bernoulli",
                           tree.complexity = 5,
                           learning.rate=0.005)
      }
      if(brt_surv$gbm.call$best.trees<1000){
        brt_surv<-gbm.step(data=surv_dat_train, gbm.x = 2:length(surv_dat_train), gbm.y=1, 
                           family = "bernoulli",
                           tree.complexity = 5,
                           learning.rate=0.001)
      }
      
      
      
      
      # 3. predict at veg survey sites
      #--------------------------------------
      
      brt_predict_pres[[i]]<- predict(brt_pres, veg_data[,-c(1)],  
                                      n.trees=brt_pres$gbm.call$best.trees, type="response")
      brt_predict_surv[[i]]<- predict(brt_surv, veg_data[,-c(1)], 
                                 n.trees=brt_surv$gbm.call$best.trees,type='response')
      # create df of test observations and predictions to evaluate model performance
      brt_predict_pres[[i]] <- data.frame(id=veg_data$id, pred=brt_predict_pres[[i]])%>%
        rename_with(~paste0(.x,i, recycle0 = TRUE), "pred")
      brt_predict_surv[[i]] <- data.frame(id=veg_data$id, pred=brt_predict_surv[[i]])%>%
        rename_with(~paste0(.x,i, recycle0 = TRUE), "pred")
    }

brt_ci_pres<- brt_predict_pres[[1]]
brt_ci_pres<-cbind(brt_ci_pres,do.call(cbind,lapply(brt_predict_pres[-1],function(x) x[2])))
brt_ci_pres<-brt_ci_pres%>%
  pivot_longer(cols = starts_with("pred"), names_to = "iteration",values_to = "pred")%>%
  group_by(id)%>%
  summarise(mean_pres_pred=mean(pred,na.rm=T),
            u.ci_pres=quantile(pred,probs = 0.975),
            l.ci_pres=quantile(pred,probs=0.025)
  )
    
brt_ci_surv<- brt_predict_surv[[1]]
brt_ci_surv<-cbind(brt_ci_surv,do.call(cbind,lapply(brt_predict_surv[-1],function(x) x[2])))
brt_ci_surv<-brt_ci_surv%>%
  pivot_longer(cols = starts_with("pred"), names_to = "iteration",values_to = "pred")%>%
  group_by(id)%>%
  summarise(mean_surv_pred=mean(pred,na.rm=T),
            u.ci_surv=quantile(pred,probs = 0.975),
            l.ci_surv=quantile(pred,probs=0.025)
  )
veg_data2<-veg_data%>%
  left_join(brt_ci_surv,by="id")%>%
  left_join(brt_ci_pres,by="id")


#Then make predictions using the full model
load(paste0(path_out,"Final_outputs/Nest_Predictions/SALS_final_BRT_mods.RDS"))


# 3. predict at veg survey sites
#--------------------------------------

veg_data2$pres_pred<- predict(brt_pres, veg_data2[,-c(1,14:19)],  
                                n.trees=brt_pres$gbm.call$best.trees, type="response")
veg_data2$surv_pred<- predict(brt_surv, veg_data2[,-c(1,14:19)], 
                                n.trees=brt_surv$gbm.call$best.trees,type='response')

write.csv(veg_data2,paste0(path_out,"Final_outputs/Survey_Veg_Predictions/rapid_veg_BRT_predictions_CIs_15buff_single_plots.csv"),row.names = F)

  

    
    
   