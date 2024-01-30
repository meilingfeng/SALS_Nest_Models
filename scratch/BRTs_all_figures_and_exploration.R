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

#Load and tidy data
if(!exists("pres_dat")){
  source("05a_Set_Model_Parameters.R")
}

if(build==T){
  #Create objects to hold BRT parameters for each fold
  n.tree.p<-c()
  n.tree.s<-c()
  lr.p<-c()
  lr.s<-c()
  contrib_pres<-list()
  contrib_surv<-list()
  var_plots_pres<-list()
  var_plots_surv<-list()
  int.pres<-list()
  int.surv<-list()
  int.plots.pres<-list()
  int.plots.surv<-list()
  d.pres.brt<-list()
  d.surv.brt<-list()
  
  
  # select all predictors and response
  pres_dat2<-pres_dat%>%dplyr::select("id","y","group",all_of(all_terms))
  surv_dat2<-surv_dat%>%dplyr::select("id","y","group",all_of(all_terms))
  
  
  set.seed(123)
  for(i in 1:k){
    #divide into testing and training
    pres_dat_train<-pres_dat2[pres_dat2$group!=i,]%>%dplyr::select(-group,-id)
    surv_dat_train<-surv_dat2[surv_dat2$group!=i,]%>%dplyr::select(-group,-id)
    
    pres_dat_test<-pres_dat2[pres_dat2$group==i,]%>%dplyr::select(-group,-id)
    surv_dat_test<-surv_dat2[surv_dat2$group==i,]%>%dplyr::select(-group,-id)
    
    
    # 2. determine optimal number of trees
    #---------------------------------------------
    
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
                       learning.rate=0.005)
    
    # start at lr of 0.05, decrease until optimal trees hits 1000
    if(brt_surv$gbm.call$best.trees<1000){
      brt_surv<-gbm.step(data=surv_dat_train, gbm.x = 2:length(surv_dat_train), gbm.y=1, 
                         family = "bernoulli",
                         tree.complexity = 5,
                         learning.rate=0.001)
    }
    
    
    
    # list the optimal number of trees and lr for the final models
    n.tree.p[i]<-brt_pres$gbm.call$best.trees
    n.tree.s[i]<-brt_surv$gbm.call$best.trees
    
    lr.p[i]<-brt_pres$gbm.call$learning.rate
    lr.s[i]<-brt_surv$gbm.call$learning.rate
    
    # 2. plot fitted functions
    #------------------------------------
    #plot fitted values against predictor
    gbm.plot(brt_pres, n.plots=8, plot.layout=c(2, 4), write.title = FALSE, smooth=T) #, continuous.resolution=100, type="response"
    var_plots_pres[[i]]<-recordPlot()
    gbm.plot(brt_surv, n.plots=8, plot.layout=c(2, 4), write.title = FALSE, smooth=T)
    var_plots_pres[[i]]<-recordPlot()
    
    #gbm.plot.fits(brt_pres)
    #gbm.plot.fits(brt_surv)
    
    #Variable contribution/importance
    summary(brt_pres)
    summary(brt_surv)
    contrib_pres[[i]]<-brt_pres$contributions
    contrib_surv[[i]]<-brt_surv$contributions
    
    
    #look for interactions (tree complexity)
    find.int.pres <- gbm.interactions(brt_pres)
    find.int.surv <- gbm.interactions(brt_surv)
    #lists the top interactions and plots them
    # for presence
    int.pres[[i]]<-find.int.pres$rank.list 
    var1<-find.int.pres$rank.list[1,"var1.index"]
    var2<-find.int.pres$rank.list[1,"var2.index"]
    gbm.perspec(brt_pres, var1, var2) #3d plots the interactions with fitted values
    int.plots.pres[[i]]<-recordPlot()
    # for survival
    int.surv[[i]]<-find.int.surv$rank.list 
    var1<-find.int.surv$rank.list[1,"var1.index"]
    var2<-find.int.surv$rank.list[1,"var2.index"]
    gbm.perspec(brt_surv, var1, var2)
    int.plots.surv[[i]]<-recordPlot()
    
    
    # predict to testing data
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



if(predict.surf==T){
  
  #Final model for nest presence
  set.seed(123)
  pres_dat2<-pres_dat%>%dplyr::select("id","y",all_of(all_terms))
  brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate=0.05)# default bag fraction is 0.75
  #start at lr of 0.05, decrease until optimal trees hits 1000
  if(brt_pres$gbm.call$best.trees<1000){
    brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1,
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.01)
  }
  if(brt_pres$gbm.call$best.trees<1000){
    brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.005)
  }
  if(brt_pres$gbm.call$best.trees<1000){
    brt_pres<-gbm.step(data=pres_dat2[,-1], gbm.x = 2:length(pres_dat2[,-1]), gbm.y=1, 
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.001)
  }
  
  
  # Final model for nest survival
  set.seed(123)
  surv_dat2<-surv_dat%>%dplyr::select("id","y",all_of(all_terms))
  brt_surv<-gbm.step(data=surv_dat2[,-1], gbm.x = 2:length(surv_dat2[,-1]), gbm.y=1, 
                     family = "bernoulli",
                     tree.complexity = 5,
                     learning.rate=0.005)
  
  # start at lr of 0.05, decrease until optimal trees hits 1000
  if(brt_surv$gbm.call$best.trees<1000){
    brt_surv<-gbm.step(data=surv_dat2[,-1], gbm.x = 2:length(surv_dat2[,-1]), gbm.y=1, 
                       family = "bernoulli",
                       tree.complexity = 5,
                       learning.rate=0.001)
  }
  
  
  
  # save final models and their thresholds
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
  thr.p.brt<-optimal.thresholds(pres.brt.predict,opt.methods = "MaxSens+Spec")$pred
  thr.s.brt<-optimal.thresholds(surv.brt.predict,opt.methods = "MaxSens+Spec")$pred
  
  save(brt_pres,brt_surv,thr.p.brt,thr.s.brt,file=paste0(path_out,"Final_outputs/Nest_Predictions/final_BRT_mods.RDS"))
  
  
  if(!file.exists(paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",1,"_pres_BRTpreds_30m",ab_type,".tif"))){    
    # predict to spatial surface
    #empty list to hold prediction surfaces
    brt_predict_pres<-list()
    brt_predict_surv<-list()
    
    
    mat_p<-list()
    mat_s<-list()
    mat_p_z1<-list()
    mat_s_z1<-list()
    
    #Predict across each zone
    for (i in 1:length(file_list_all_zones)){
      
      # a) load raster predictor layers for a particular zone
      predictors<-rast(unlist(file_list_all_zones[[i]]))
      
      # b) name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
      names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","elevation","uvvr_diff","uvvr_mean","precip","tideres","HIMARSH")
      
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
      writeRaster(brt_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_BRTpreds_30m",ab_type,".tif"),overwrite=T)
      writeRaster(brt_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_BRTpreds_30m",ab_type,".tif"),overwrite=T)
      
      
      # g) Write predicted values to csv
      # if(i!=1){
      # pres_out<-rast(list(preds_mask2,brt_predict_pres[[i]]))
      # names(pres_out[[nlyr(pres_out)]])<-"predictions"
      # surv_out<-rast(list(preds_mask2,brt_predict_surv[[i]]))
      # names(surv_out[[nlyr(surv_out)]])<-"predictions"
      
      #  mat_p[[i]]<-as.data.frame(terra::as.matrix(pres_out,wide=F))%>%
      #    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      #  mat_s[[i]]<-as.data.frame(terra::as.matrix(surv_out,wide=F))%>%
      #    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      #  }
      
      #zone 1 is too big, process it in parts
      #  if(i==1){
      
      #   n_div<-5
      #   dat_zone_list<-list()
      #   row_start_p<-c()
      #   row_end_p<-c()
      #   row_start_s<-c()
      #   row_end_s<-c()
      
      #   pres_out<-rast(list(preds_mask2,brt_predict_pres[[i]]))
      #   names(pres_out[[nlyr(pres_out)]])<-"predictions"
      
      
      #   surv_out<-rast(list(preds_mask2,brt_predict_surv[[i]]))
      #   names(surv_out[[nlyr(surv_out)]])<-"predictions"
      
      #divide zone into 4
      #   for(j in 1:n_div){
      #     row_start_p[j]<-((round(nrow(pres_out)/n_div))*(j-1))+1
      #     row_end_p[j]<-(round(nrow(pres_out)/n_div))*j
      #     pres_out_sub<-pres_out[c(row_start_p[j]:row_end_p[j]),c(1:ncol(pres_out)),drop=F]
      
      #    row_start_s[j]<-((round(nrow(surv_out)/n_div))*(j-1))+1
      #    row_end_s[j]<-(round(nrow(surv_out)/n_div))*j
      #    surv_out_sub<-surv_out[c(row_start_s[j]:row_end_s[j]),c(1:ncol(surv_out)),drop=F]
      
      #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
      #    mat_p_z1[[j]]<-as.data.frame(terra::as.matrix(pres_out_sub,wide=F))%>%
      #      mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      
      #     mat_s_z1[[j]]<-as.data.frame(terra::as.matrix(surv_out_sub,wide=F))%>%
      #       mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      #  }
      
      
      #bind all the zone sections into one dataframe
      #   mat_p[[i]]<-do.call("rbind",mat_p_z1)
      #  mat_s[[i]]<-do.call("rbind",mat_s_z1)
      #}
      
      
      #write each zone prediction dataset to matrix file
      #  write.csv(mat_p[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Placement/Z",i,"_BRTprediction_30m_",ab_type,"_placement.csv"),row.names = F)
      #  write.csv(mat_s[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Success/Z",i,"_BRTprediction_30m_",ab_type,"_success.csv"),row.names = F)
    }
    
  }
}