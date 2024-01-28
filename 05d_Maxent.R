library(dismo)
library(ggpubr)
library(patchwork)
# test if you can use maxent 
maxent()

### Set up
# -------------------------------------------
#Load and tidy occurrence data
if(!exists("pres_dat")){
  source("05a_Set_Model_Parameters.R")
} 

if(build==T){
# select all predictors and response
pres_dat2<-pres_dat%>%dplyr::select("id","y","group",all_of(all_terms))
surv_dat2<-surv_dat%>%dplyr::select("id","y","group",all_of(all_terms))

# object to hold the predictions from each fold
d.pres.mxt<-list()
d.surv.mxt<-list()

# For each fold, divide occurrence data into testing and training datasets
for(i in 1:k){
  
  pres_dat_train<-pres_dat2[pres_dat2$group!=i,]%>%dplyr::select(-group,-id) #select all groups but 1
  surv_dat_train<-surv_dat2[surv_dat2$group!=i,]%>%dplyr::select(-group,-id)
  
  pres_dat_test<-pres_dat2[pres_dat2$group==i,]%>%dplyr::select(-group,-id) #select the one group out for testing
  surv_dat_test<-surv_dat2[surv_dat2$group==i,]%>%dplyr::select(-group,-id)


### for each zone...





### Fit model
#------------------------------------------------------------
  me_p <- maxent(dplyr::select(pres_dat_train,all_of(all_terms)), # predictors (to fit categorical variables, indicate them using ",factors='variablename'")
                 dplyr::select(pres_dat_train,y)) # occurrences
  
  me_s <- maxent(dplyr::select(surv_dat_train,all_of(all_terms)), 
                 dplyr::select(surv_dat_train,y))
  

 
  
### Predict to testing data
#--------------------------------------------------
# Store model predictions for each fold
  # create df of test observations and predictions to evaluate model performance
  d.pres.mxt[[i]] <- data.frame(id=pres_dat2[pres_dat2$group==i,]$id,
                            obs=pres_dat_test$y, 
                            pred=predict(me_p, dplyr::select(pres_dat_test,all_of(all_terms))))
  d.surv.mxt[[i]] <- data.frame(id=surv_dat2[surv_dat2$group==i,]$id,
                            obs=surv_dat_test$y, 
                            pred=predict(me_s, dplyr::select(surv_dat_test,all_of(all_terms))))
  
  
  
## Model Parameters
#------------------------------------------------------------------------
  var_imp_p<-list()
  var_imp_s<-list()
  me_acc_p<-list()
  me_acc_s<-list()
  auc_p<-list()
  auc_s<-list()
  

  
# plot showing importance of each variable
  var_imp_s[[i]]<-plot(me_s)
  var_imp_p[[i]]<-plot(me_p)
# response curves
 # response(me_s)
 # response(me_p)

}
  
}


if(predict.surf==T){  
if(!file.exists(paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",1,"_pres_MXTpreds_30m",ab_type,".tif"))){  
  
### Predict to spatial grid
#--------------------------------------------------
  mxt_predict_pres<-list()
  mxt_predict_surv<-list()
for (j in 1:length(file_list_all_zones)){
    
    # a) load raster predictor layers for a particular zone
    predictors<-rast(unlist(file_list_all_zones[[j]]))
    
    # b) name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
    names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","tideres","HIMARSH", "LOMARSH")
    
    # c) select just the layers that were used as predictor variables in the model
    mod_preds<-predictors[[names(predictors)%in%c("Highmarsh",all_terms)]]
    
    # d) set areas outside marsh to NA (mask all predictors with marsh area)
    mask<-predictors["Highmarsh"]
    mask[mask==0|mask==9|mask==7|mask==8]<-NA
    preds_mask<-list()
    for(n in 1:nlyr(mod_preds)){
      preds_mask[[n]]<-mask(mod_preds[[n]],mask)
    }
    
    
    #turn back from list to raster stack
    preds_mask2<-rast(preds_mask)
    names(preds_mask2)<-names(mod_preds)
    
    
    
    # e) PREDICT 
    mxt_predict_pres[[j]]<- mask(predict(preds_mask2[[-1]], me_p, progress=''),
                                 mask)
    mxt_predict_surv[[j]]<- mask(predict(preds_mask2[[-1]], me_s, progress=''),
                                 mask)
    
    
    
    
    # g) Write predicted values to csv
    if(j!=1){
      pres_out<-rast(list(preds_mask2,mxt_predict_pres[[j]]))
      names(pres_out[[nlyr(pres_out)]])<-"predictions"
      mat_p[[j]]<-as.data.frame(terra::as.matrix(pres_out,wide=F))%>%
        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      
      surv_out<-rast(list(preds_mask2,mxt_predict_surv[[j]]))
      names(surv_out[[nlyr(surv_out)]])<-"predictions"
      mat_s[[j]]<-as.data.frame(terra::as.matrix(surv_out,wide=F))%>%
        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
    }
    
    #zone 1 is too big, process it in parts
    if(j==1){
      
      n_div<-4
      dat_zone_list<-list()
      row_start<-c()
      row_end<-c()
      
      
      pres_out<-rast(list(preds_mask2,mxt_predict_pres[[j]]))
      names(pres_out[[nlyr(pres_out)]])<-"predictions"
      
      
      surv_out<-rast(list(preds_mask2,mxt_predict_surv[[j]]))
      names(surv_out[[nlyr(surv_out)]])<-"predictions"
      
      #divide zone into 4
      for(z in 1:n_div){
        row_start_p[z]<-((round(nrow(pres_out)/n_div))*(z-1))+1
        row_end_p[z]<-(round(nrow(pres_out)/n_div))*z
        pres_out_sub<-pres_out[c(row_start_p[z]:row_end_p[z]),c(1:ncol(pres_out)),drop=F]
        
        row_start_s[z]<-((round(nrow(surv_out)/n_div))*(z-1))+1
        row_end_s[z]<-(round(nrow(surv_out)/n_div))*z
        surv_out_sub<-surv_out[c(row_start_s[z]:row_end_s[z]),c(1:ncol(surv_out)),drop=F]
        
        #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
        mat_p_z1[[z]]<-as.data.frame(terra::as.matrix(pres_out_sub,wide=F))%>%
          mutate(id=paste0(seq(1,nrow(.),1),"z",i))
        
        mat_s_z1[[z]]<-as.data.frame(terra::as.matrix(surv_out_sub,wide=F))%>%
          mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      }
      
      
      #bind all the zone sections into one dataframe
      mat_p[[j]]<-do.call("rbind",mat_p_z1)
      mat_s[[j]]<-do.call("rbind",mat_s_z1)
    }
    
    
    #write each zone prediction dataset to matrix file
    write.csv(mat_p[[j]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Placement/Z",j,"_MXTprediction_30m_",ab_type,"_placement.csv"),row.names = F)
    write.csv(mat_s[[j]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Success/Z",j,"_MXTprediction_30m_",ab_type,"_success.csv"),row.names = F)
    
    #Write rasters
    writeRaster(mxt_predict_pres[[j]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",j,"_pres_MXTpreds_30m",ab_type,".tif"),overwrite=T)
    writeRaster(mxt_predict_surv[[j]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",j,"_surv_MXTpreds_30m",ab_type,".tif"),overwrite=T)
    
    
} 
}
}