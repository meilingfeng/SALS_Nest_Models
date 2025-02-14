library(tidyverse)
library(car)#vif
library(AICcmodavg)
library(bbmle)
library(broom)
library(ggstats)

### Set up
# -------------------------------------------
source("C:/Users/mefen/OneDrive/Documents/Github/SALS_Nest_Models/05b_GLM_Diagnostics_ModelBuilding.R")

## Model interpretation
#-----------------------------------------------------------------------------

#Things to test with ANOVA:

# Are we missing information about nesting habitat characteristics in classified RS variables? How much variance does raw reflectance (PC1) explain?
# Does vegetation describe flood risk better than elevation? - expect elevation to be more important north (Ruskin)
  # how much variance does elevation explain? Compared to vegetation metrics?

pres_anova<-as.data.frame(Anova(mod_list_pres[[9]],type=3))
write.csv(pres_anova,paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/pres_anova.csv"))


surv_anova<-as.data.frame(Anova(mod_list_surv[[9]],type=3))
write.csv(surv_anova,paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/surv_anova.csv"))

# Elevation played less of a role in nest success compared to nest site selection.
# Low marsh vegetative community explained a significant amount of variation in nest success.
# unsurprisingly, time since the last spring tide and day of the year also explain much of the variation in nesting success.
# One of the main drivers of nest failure is flooding, which was supported by the importance of spring tide timing and DOY.
# The increased importance of vegetation communities relative to elevation as seen in the site selection models
# could suggest that vegetation is currently a more accurate predictor of flood conditions on the marsh.
# An important change in elevation for a tidal nesting specialist like SALS is a matter of centimeters above or below the high tide line.



#check out coefficients
summary(mod_list_pres[[9]])
summary(mod_list_surv[[9]])



#format the coefficients and confidence intervals and plot

#ggcoef_table(mod_list_pres[[9]],exponentiate = T)
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/glm_pres_coefs.png"), width = 10, height = 7, dpi = "retina")

#ggcoef_table(mod_list_surv[[9]],exponentiate = T)
#ggsave(filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/glm_surv_coefs.png"), width = 10, height = 7, dpi = "retina")



## 2. Model Comparison 
#------------------------------------------

# model comparison table
#want to report AIC or AICc if small sample or BIC if large sample (better at parsimony than AIC)
#also AIC weight, how many times that model is likely to give the best predictions on new data
#also redicual deviance and df for each model provide rough goodness of fit.
#https://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/

model.names<-c("Null","Elevation",
               "Vegetation Communnities","Within-Marsh: Elevation","Within-Marsh: Brightness", "Within-Marsh: NDVI","Within-Marsh: Global",
               "Between-Marsh","Global")
#For presence
summ.table <- do.call(rbind, lapply(mod_list_pres, broom::glance))
mod_tab_p<-ICtab(mod_list_pres, type="AIC", weights=T, delta=T,nobs=nrow(pres_dat),sort=F,mnames=model.names)
mod_tab_p[["Resid.Dev"]]<-summ.table[["deviance"]]
mod_tab_p[["Null.Dev"]]<-summ.table[["null.deviance"]]
mod_tab_p[["Response"]]<-rep("Presence",nrow(summ.table))
mod_tab_p<-as.data.frame(mod_tab_p)%>%
  mutate(across(-Response,\(x) round(x, 2)),
         Pseudo_R2= 1-(round(Resid.Dev,4)/round(Null.Dev,4)))

# for success
summ.table <- do.call(rbind, lapply(mod_list_surv, broom::glance))
mod_tab_s<-ICtab(mod_list_surv, type="AIC", weights=T, delta=T,nobs=nrow(surv_dat),sort=F,mnames=model.names)
mod_tab_s[["Resid.Dev"]]<-summ.table[["deviance"]]
mod_tab_s[["Null.Dev"]]<-summ.table[["null.deviance"]]
mod_tab_s[["Response"]]<-rep("Fate",nrow(summ.table))
  #calculate McFaddens R2
mod_tab_s<-as.data.frame(mod_tab_s)%>%
  mutate(across(-Response,\(x) round(x, 2)),
         Pseudo_R2= 1-(round(Resid.Dev,4)/round(Null.Dev,4)))

mod_tab<-rbind(mod_tab_p,mod_tab_s)%>%
  arrange(Response,dAIC)



if(!file.exists(paste0(path_out,"Final_outputs/Model_Results/model_selection_table_2_13_25.csv"))){
write.csv(mod_tab,paste0(path_out,"Final_outputs/Model_Results/model_selection_table_2_13_25.csv"))
}


# select the top ranked model based on AIC
  # presence models
mod_pres<-mod_list_pres[[#in the list of models
  which( #select the index of the model that
    model.names%in%row.names(mod_tab[mod_tab$dAIC==0 & mod_tab$Response=="Presence",])# is has the model name as the row with the lowest AIC for presence
    )
  ]]
form_pres<-deparse1(mod_pres$formula)# use deparse1 to turn the formula into a string

  # repeat for the success models
mod_surv<-mod_list_surv[[which(model.names%in%str_sub(row.names(mod_tab[mod_tab$dAIC==0 & mod_tab$Response=="Fate",]),end=-2))]]
form_surv<-deparse1(mod_surv$formula)




if(build==T){
## 2. Predict to test data
#--------------------------------------

d.pres.glm<-list()
d.surv.glm<-list()

for (i in 1:k) {
  train <- pres_dat[pres_dat$group != i,]
  test <- pres_dat[pres_dat$group == i,]
  mod.p <- glm(form_pres, data = train, family = binomial(link = "logit"))
  d.pres.glm[[i]] <- data.frame(id=test$id,
                                obs=test$y, 
                                pred=predict(mod.p,test, type="response"))
}

for (i in 1:k) {
  train <- surv_dat[surv_dat$group != i,]
  test <- surv_dat[surv_dat$group == i,]
  mod.s <- glm(form_surv, data = train, family = binomial(link = "logit"))
  d.surv.glm[[i]] <- data.frame(id=test$id,
                                obs=test$y, 
                                pred=predict(mod.s,test, type="response"))
}

}



if(predict.surf==T){
## 3. Predict to Spatial Surface
#-------------------------------------------
final_mod_pres <- glm(form_pres, data = pres_dat, family = binomial(link = "logit"))
final_mod_surv <- glm(form_surv, data = surv_dat, family = binomial(link = "logit"))

if(!file.exists(paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",1,"_pres_GLMpreds_30m",ab_type,".tif"))){
# a) list predictor variables in the selected models
terms_pres<-unlist(strsplit(as.character(substring(form_pres,5)),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()%>%
  unique()

terms_surv<-unlist(strsplit(as.character(substring(form_surv,5)),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()%>%
  unique()


# b) create empty lists to hold prediction surfaces
glm_predict_pres<-list()
glm_predict_surv<-list()


mat_p<-list()
mat_s<-list()
mat_p_z1<-list()
mat_s_z1<-list()

# c) loop predictions through each zone
for (i in 1:length(file_list_all_zones)){
  
  # a) load raster predictor layers for a particular zone
  predictors<-rast(unlist(file_list_all_zones[[i]]))
  
  # b) name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
  names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","elevation","uvvr_diff","uvvr_mean","precip","tideres","HIMARSH","LOMARSH")

  # c) select just the layers that were used as predictor variables in the model
  mod_preds_p<-predictors[[names(predictors)%in%terms_pres]]
  mod_preds_s<-predictors[[names(predictors)%in%terms_surv]]
  
  # d) create a latitude layer if it was selected in the final model
  if("latitude"%in%terms_pres){
    lat<-rast(crs=crs(mod_preds_p[[1]]),extent=ext(mod_preds_p[[1]]),resolution=res(mod_preds_p[[1]]),vals=xyFromCell(terra::project(mod_preds_p[[1]],"EPSG:4617"),1:ncell(mod_preds_p[[1]]))[,2])
    # add to predictor list
    mod_preds_p<-rast(list(mod_preds_p,lat))
    names(mod_preds_p[[nlyr(mod_preds_p)]])<-"latitude"
  }
  if("latitude"%in%terms_surv){
    lat<-rast(crs=crs(mod_preds_s[[1]]),extent=ext(mod_preds_s[[1]]),resolution=res(mod_preds_s[[1]]),vals=xyFromCell(terra::project(mod_preds_s[[1]],"EPSG:4617"),1:ncell(mod_preds_s[[1]]))[,2])
    # add to predictor list
    mod_preds_s<-rast(list(mod_preds_s,lat))
    names(mod_preds_s[[nlyr(mod_preds_s)]])<-"latitude"
  }
  
  # e) set upland and open water to NA (mask all predictors)
  mask<-predictors["Highmarsh"]
  mask[mask==0|mask==9|mask==7]<-NA
  preds_mask_s<-list()
  preds_mask_p<-list()
  for(j in 1:nlyr(mod_preds_s)){
    preds_mask_s[[j]]<-mask(mod_preds_s[[j]],mask)
  }
  
  for(j in 1:nlyr(mod_preds_p)){
    preds_mask_p[[j]]<-mask(mod_preds_p[[j]],mask)
  }
  
  
  
  #turn back from list to raster stack
  preds_mask_p2<-rast(preds_mask_p)
  names(preds_mask_p2)<-names(mod_preds_p)
  
  preds_mask_s2<-rast(preds_mask_s)
  names(preds_mask_s2)<-names(mod_preds_s)
  
  #PREDICT
  glm_predict_pres[[i]]<-predict(preds_mask_p2,final_mod_pres, type="response")
  glm_predict_surv[[i]]<-predict(preds_mask_s2,final_mod_surv, type="response")
  
  #Write predicted values to rasters
  writeRaster(glm_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_predsGLM_30m",ab_type,".tif"),overwrite=T)
  writeRaster(glm_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_predsGLM_30m",ab_type,".tif"),overwrite=T)
  
  #Write predicted values to csv
#  pres_out<-rast(list(preds_mask_p2,glm_predict_pres[[i]]))
#  names(pres_out[[nlyr(pres_out)]])<-"predictions"
#  mat_p[[i]]<-as.data.frame(terra::as.matrix(pres_out,wide=F))%>%
#    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
  
#  surv_out<-rast(list(preds_mask_s2,glm_predict_surv[[i]]))
#  names(surv_out[[nlyr(surv_out)]])<-"predictions"
#  mat_s[[i]]<-as.data.frame(terra::as.matrix(surv_out,wide=F))%>%
#    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
  
  #zone 1 is too big, process it in parts
#  if(i==1){
    
#    n_div<-4
#    dat_zone_list<-list()
#    row_start<-c()
#    row_end<-c()
    
    
#    pres_out<-rast(list(preds_mask_p2,glm_predict_pres[[i]]))
#    names(pres_out[[nlyr(pres_out)]])<-"predictions"
    
    
#    surv_out<-rast(list(preds_mask_s2,glm_predict_surv[[i]]))
#    names(surv_out[[nlyr(surv_out)]])<-"predictions"
    
#    #divide zone into 4
#    for(j in 1:n_div){
#      row_start_p[j]<-((round(nrow(pres_out)/n_div))*(j-1))+1
#      row_end_p[j]<-(round(nrow(pres_out)/n_div))*j
#      pres_out_sub<-pres_out[c(row_start_p[j]:row_end_p[j]),c(1:ncol(pres_out)),drop=F]
      
#      row_start_s[j]<-((round(nrow(surv_out)/n_div))*(j-1))+1
#      row_end_s[j]<-(round(nrow(surv_out)/n_div))*j
#      surv_out_sub<-surv_out[c(row_start_s[j]:row_end_s[j]),c(1:ncol(surv_out)),drop=F]
      
      #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
#      mat_p_z1[[j]]<-as.data.frame(terra::as.matrix(pres_out_sub,wide=F))%>%
#        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
#      
#      mat_s_z1[[j]]<-as.data.frame(terra::as.matrix(surv_out_sub,wide=F))%>%
#        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
#    }
    
    
    #bind all the zone sections into one dataframe
#    mat_p[[1]]<-do.call("rbind",mat_p_z1)
#    mat_s[[1]]<-do.call("rbind",mat_s_z1)
#  }
  
  
  #write each zone prediction dataset to matrix file
#  write.csv(mat_p[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Placement/Z",i,"_predictionGLM_30m_",ab_type,"_placement.csv"),row.names = F)
#  write.csv(mat_s[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Success/Z",i,"_predictionGLM_30m_",ab_type,"_success.csv"),row.names = F)
  
}
}
}