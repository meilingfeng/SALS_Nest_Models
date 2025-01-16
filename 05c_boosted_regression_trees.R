library(tidyverse)
library(car)
library(gbm)
library(raster)
library(terra)
library(dismo)
library(PresenceAbsence)
library(patchwork)


## 1. Set up
#--------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

#Load and tidy data
all_terms<-c("uvvr_mean","ndvi","pca","HIMARSH","LOMARSH", "tideres", "uvvr_diff","elevation") 

#if running all species, set sals_only to FALSE.
if(sals_only==T){
speciesnames<-c("SALS")
}else{
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
}

for (s in 1:length(speciesnames)){
  pres_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[s],"_nest_pres_dat.csv"))
  surv_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[s],"_nest_surv_dat.csv"))
      
if(build==T){
#Create objects to hold BRT parameters for each fold
n.tree.p<-c()
n.tree.s<-c()
lr.p<-c()
lr.s<-c()
contrib_pres<-list()
contrib_surv<-list()
d.pres.brt<-list()
d.surv.brt<-list()


# select all predictors and response
pres_dat2<-pres_dat%>%dplyr::select("id","y","group",all_of(all_terms))%>%
  rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
         Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH)
surv_dat2<-surv_dat%>%dplyr::select("id","y","group",all_of(all_terms))%>%
  rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
         Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH)


for(i in 1:k){
  #divide into testing and training
  pres_dat_train<-pres_dat2[pres_dat2$group!=i,]%>%dplyr::select(-group,-id)
  surv_dat_train<-surv_dat2[surv_dat2$group!=i,]%>%dplyr::select(-group,-id)
  
  pres_dat_test<-pres_dat2[pres_dat2$group==i,]%>%dplyr::select(-group,-id)
  surv_dat_test<-surv_dat2[surv_dat2$group==i,]%>%dplyr::select(-group,-id)


# 2. determine optimal number of trees
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
              learning.rate=0.001)

    # start at lr of 0.05, decrease until optimal trees hits 1000
#if(brt_surv$gbm.call$best.trees<1000){
#  brt_surv<-gbm.step(data=surv_dat_train, gbm.x = 2:length(surv_dat_train), gbm.y=1, 
 #                    family = "bernoulli",
#                     tree.complexity = 5,
#                     learning.rate=0.001)
#}



# list the optimal number of trees and lr for the final models
n.tree.p[i]<-brt_pres$gbm.call$best.trees
n.tree.s[i]<-brt_surv$gbm.call$best.trees

lr.p[i]<-brt_pres$gbm.call$learning.rate
lr.s[i]<-brt_surv$gbm.call$learning.rate

# 2. Variable contribution/importance
#------------------------------------

summary(brt_pres)
summary(brt_surv)
contrib_pres[[i]]<-brt_pres$contributions
contrib_surv[[i]]<-brt_surv$contributions


# 3. predict to testing data
#--------------------------------------
preds.pres <- predict.gbm(brt_pres, pres_dat_test,
                     n.trees=brt_pres$gbm.call$best.trees, type="response")
preds.surv <- predict.gbm(brt_surv, surv_dat_test,
                          n.trees=brt_surv$gbm.call$best.trees, type="response")



# create df of test observations and predictions to evaluate model performance
d.pres.brt[[i]] <- data.frame(id=pres_dat2[pres_dat2$group==i,]$id,
                     obs=pres_dat_test$y, 
                     pred=preds.pres)
d.surv.brt[[i]] <- data.frame(id=surv_dat2[surv_dat2$group==i,]$id,
                     obs=surv_dat_test$y, 
                     pred=preds.surv)
}

}


  
  
  
# 4. Create final predictions using all data
#----------------------------------------------------------------
if(predict.surf==T){
    
  pres_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[s],"_nest_pres_ml_dat.csv"))#use all the records since BRTs handle missing data
  surv_dat<-read.csv(paste0(path_out,"Intermediate_outputs/Nest_Datasets/",speciesnames[s],"_nest_surv_ml_dat.csv"))
  
    #Final model for nest presence
    set.seed(123)
    pres_dat2<-pres_dat%>%dplyr::select("id","y",all_of(all_terms))%>%
      rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
             Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH)
    brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.05)# default bag fraction is 0.75
    #start at lr of 0.05, decrease until optimal trees hits 1000
    if(brt_pres$gbm.call$best.trees<1000||length(brt_pres)==0){
      brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1,
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.01)
    }
    if(brt_pres$gbm.call$best.trees<1000||length(brt_pres)==0){
      brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.005)
    }
    if(brt_pres$gbm.call$best.trees<1000||length(brt_pres)==0){
      brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.001)
    }
    
    
    # Final model for nest survival
    set.seed(123)
    surv_dat2<-surv_dat%>%dplyr::select("id","y",all_of(all_terms))%>%
      rename(Elevation=elevation,NDVI=ndvi,Tidal.Restriction=tideres, Surface.Brightness=pca, Change.in.UVVR=uvvr_diff, Mean.UVVR=uvvr_mean,
             Proportion.High.Marsh=HIMARSH,Proportion.Low.Marsh=LOMARSH)
    brt_surv<-gbm.step(data=surv_dat2[,-1], gbm.x = 2:length(surv_dat2[,-1]), gbm.y=1, 
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.005)
    
    # start at lr of 0.05, decrease until optimal trees hits 1000
    if(brt_surv$gbm.call$best.trees<1000||length(brt_surv)==0){
      brt_surv<-gbm.step(data=surv_dat2[,-1], gbm.x = 2:length(surv_dat2[,-1]), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.001)
    }
      if(brt_surv$gbm.call$best.trees<1000||length(brt_surv)==0){
        brt_surv<-gbm.step(data=surv_dat2[,-1], gbm.x = 2:length(surv_dat2[,-1]), gbm.y=1, 
                           family = "bernoulli",
                           tree.complexity = 2,
                           learning.rate=0.001)
      }

    
    
    
# list the optimal number of trees and lr for the final models
brt_pres$gbm.call$best.trees
brt_surv$gbm.call$best.trees

brt_pres$gbm.call$learning.rate    
brt_surv$gbm.call$learning.rate

    
    
# save final models and their thresholds
    #plot fitted values against predictor
gbm.plot(brt_pres, n.plots=8, plot.layout=c(2, 4), write.title = FALSE, smooth=T) #, continuous.resolution=100, type="response"

gbm.plot(brt_surv, n.plots=8, plot.layout=c(2, 4), write.title = FALSE, smooth=T)

    
      # get thresholds
        # predictions
    preds.pres <- predict.gbm(brt_pres, pres_dat2,
                              n.trees=brt_pres$gbm.call$best.trees, type="response")
    preds.surv <- predict.gbm(brt_surv, surv_dat2,
                              n.trees=brt_surv$gbm.call$best.trees, type="response")
    
        # create df of test observations and predictions to evaluate model performance
    pres.brt.predict <- data.frame(id=pres_dat2$id,
                                  obs=pres_dat2$y, 
                                  pred=preds.pres)
    surv.brt.predict <- data.frame(id=surv_dat2$id,
                                  obs=surv_dat2$y, 
                                  pred=preds.surv)
    thr.p.brt<-optimal.thresholds(pres.brt.predict,opt.methods = "MaxPCC")$pred
    thr.s.brt<-optimal.thresholds(surv.brt.predict,opt.methods = "MaxPCC")$pred
    
    save(brt_pres,brt_surv,thr.p.brt,thr.s.brt,file=paste0(path_out,"Final_outputs/Nest_Predictions/",speciesnames[s],"_final_BRT_mods.RDS"))
    
    
# predict to spatial surface
#empty list to hold prediction surfaces
brt_predict_pres<-list()
brt_predict_surv<-list()


mat_p<-list()
mat_s<-list()
mat_p_z1<-list()
mat_s_z1<-list()

load(paste0(path_out,"predictor_files_all_zones_30m.rds"))

#Predict across each zone
for (i in 1:length(file_list_all_zones)){
  
  # a) load raster predictor layers for a particular zone
  predictors<-rast(unlist(file_list_all_zones[[i]]))
  
  # b) name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
  names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","elevation","uvvr_diff","uvvr_mean","precip","tideres","HIMARSH","LOMARSH")

  # c) select just the layers that were used as predictor variables in the model
  mod_preds<-predictors[[names(predictors)%in%c("Highmarsh",all_terms)]]

  # d) set areas outside marsh to NA (mask all predictors with marsh area)
  mask<-predictors["Highmarsh"]
  mask[mask==0|mask==9|mask==7]<-NA
  preds_mask<-list()
  for(j in 1:nlyr(mod_preds)){
    preds_mask[[j]]<-mask(mod_preds[[j]],mask)
  }
  
  
  #turn back from list to raster stack
  preds_mask2<-rast(preds_mask)
  names(preds_mask2)<-names(mod_preds)
  
  
  # e) PREDICT
  brt_predict_pres[[i]]<- mask(predict(preds_mask2[[-1]], brt_pres, 
                                       n.trees=brt_pres$gbm.call$best.trees, type="response"),
                               mask)
  brt_predict_surv[[i]]<- mask(predict(preds_mask2[[-1]], brt_surv, 
                                       n.trees=brt_surv$gbm.call$best.trees, type="response"),
                               mask)

  
  # f) Write predicted values to rasters
  writeRaster(brt_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/",speciesnames[s],"/",speciesnames[s],"z",i,"_pres_BRTpreds_30m.tif"),overwrite=T)
  writeRaster(brt_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/",speciesnames[s],"/",speciesnames[s],"z",i,"_surv_BRTpreds_30m.tif"),overwrite=T)
  

}
}
}
