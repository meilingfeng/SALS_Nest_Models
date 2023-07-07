library(tidyverse)
library(car)#vif
library(AICcmodavg)

### Set up
# -------------------------------------------
if(!exists("pres_dat")){
  source("05a_Set_Model_Parameters.R")
}


## 1. Select model 
#------------------------------------------


## Start a list of potential models
# A. Habitat characteristics not captured by High Marsh

# compare a model using just proportion of high marsh to models just using NDVI, Entropy, UVVR
# compare models using additive combinations of each: highmarsh+NDVI. Highmarsh+UVVR. Highmarsh+PC, 
# Compare model using combinations of non-highmarsh: NDVI+UVVR, NDVI+PCA, PCA+UVVR, PCA+UVVR+NDVI
mod_list_pres1<-list()
mod_list_surv1<-list()



# How much variation does each variable explain relative to each other?
# NDVI
mod_list_pres1[[1]]<-glm(y~ndvi, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[1]]<-glm(y~ndvi, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# PCA
mod_list_pres1[[2]]<-glm(y~pca, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[2]]<-glm(y~pca, 
                        data=surv_dat,
                        family = binomial(link="logit"))
# UVVR
mod_list_pres1[[3]]<-glm(y~uvvr_mean, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[3]]<-glm(y~uvvr_diff, 
                         data=surv_dat,
                         family = binomial(link="logit"))
# texture
mod_list_pres1[[4]]<-glm(y~ent_txt, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[4]]<-glm(y~ent_txt, 
                         data=surv_dat,
                         family = binomial(link="logit"))
mod_list_pres1[[5]]<-glm(y~cor_txt, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[5]]<-glm(y~cor_txt, 
                         data=surv_dat,
                         family = binomial(link="logit"))
# tide restriction
mod_list_pres1[[6]]<-glm(y~tideres, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[6]]<-glm(y~tideres, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# HIGH MARSH
mod_list_pres1[[7]]<-glm(y~HIMARSH, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[7]]<-glm(y~HIMARSH, 
                        data=surv_dat,
                        family = binomial(link="logit"))

#Does adding NDVI or UVVR or PCA improve high marsh?
# HIGH MARSH +NDVI
mod_list_pres1[[8]]<-glm(y~HIMARSH+ndvi, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[8]]<-glm(y~HIMARSH+ndvi, 
                         data=surv_dat,
                         family = binomial(link="logit"))
# HIGH MARSH + UVVR
mod_list_pres1[[9]]<-glm(y~HIMARSH+uvvr_mean, 
                          data=pres_dat,
                          family = binomial(link="logit"))
mod_list_surv1[[9]]<-glm(y~HIMARSH+uvvr_diff, 
                          data=surv_dat,
                          family = binomial(link="logit"))
# HIGH MARSH + PCA
mod_list_pres1[[10]]<-glm(y~HIMARSH+pca, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[10]]<-glm(y~HIMARSH+pca, 
                         data=surv_dat,
                         family = binomial(link="logit"))
# HIGH MARSH + texture
mod_list_pres1[[11]]<-glm(y~HIMARSH+ent_txt, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[11]]<-glm(y~HIMARSH+cor_txt, 
                         data=surv_dat,
                         family = binomial(link="logit"))
mod_list_pres1[[12]]<-glm(y~HIMARSH+ent_txt, 
                          data=pres_dat,
                          family = binomial(link="logit"))
mod_list_surv1[[12]]<-glm(y~HIMARSH+cor_txt, 
                          data=surv_dat,
                          family = binomial(link="logit"))
#HIGH MARSH + TIDE RESTRICTIONS
mod_list_pres1[[13]]<-glm(y~HIMARSH+ent_txt, 
                          data=pres_dat,
                          family = binomial(link="logit"))
mod_list_surv1[[13]]<-glm(y~HIMARSH+cor_txt, 
                          data=surv_dat,
                          family = binomial(link="logit"))

# ALL
mod_list_pres1[[14]]<-glm(y~ndvi+uvvr_mean+pca+ent_txt+cor_txt+tideres+HIMARSH, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[14]]<-glm(y~ndvi+uvvr_diff+pca+cor_txt+HIMARSH,
                         data=surv_dat,
                         family = binomial(link="logit"))



# store model selection factors in a table
mod_tab_habstr<-data.frame(Response=rep(c("Presence","Success"),length(mod_list_pres1)),
                           Model_Comparison=rep("Additional Habitat Characteristics",2*length(mod_list_pres1)),
                           Model_Name=rep(NA,2*length(mod_list_pres1)),
                           Function=rep(NA,2*length(mod_list_pres1)),
                           AIC=rep(NA,2*length(mod_list_pres1)),
                           dAIC=rep(NA,2*length(mod_list_pres1)))

mod_tab_habstr[1,"Model_Name"]<-"Vegetation Health (NDVI)"
mod_tab_habstr[1,"Function"]<-deparse1(mod_list_pres1[[1]]$formula)
mod_tab_habstr[1,"AIC"]<-mod_list_pres1[[1]]$aic

mod_tab_habstr[2,"Model_Name"]<-"Vegetation Health (NDVI)"
mod_tab_habstr[2,"Function"]<-deparse1(mod_list_surv1[[1]]$formula)
mod_tab_habstr[2,"AIC"]<-mod_list_surv1[[1]]$aic

mod_tab_habstr[3,"Model_Name"]<-"Raw Reflectance (PCA)"
mod_tab_habstr[3,"Function"]<-deparse1(mod_list_pres1[[2]]$formula)
mod_tab_habstr[3,"AIC"]<-mod_list_pres1[[2]]$aic

mod_tab_habstr[4,"Model_Name"]<-"Raw Reflectance (PCA)"
mod_tab_habstr[4,"Function"]<-deparse1(mod_list_surv1[[2]]$formula)
mod_tab_habstr[4,"AIC"]<-mod_list_surv1[[2]]$aic

mod_tab_habstr[5,"Model_Name"]<-"Marsh Resilience (UVVR)"
mod_tab_habstr[5,"Function"]<-deparse1(mod_list_pres1[[3]]$formula)
mod_tab_habstr[5,"AIC"]<-mod_list_pres1[[3]]$aic

mod_tab_habstr[6,"Model_Name"]<-"Marsh Resilience (UVVR)"
mod_tab_habstr[6,"Function"]<-deparse1(mod_list_surv1[[3]]$formula)
mod_tab_habstr[6,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_habstr[7,"Model_Name"]<-"Texture (Dissimiarlity)"
mod_tab_habstr[7,"Function"]<-deparse1(mod_list_pres1[[4]]$formula)
mod_tab_habstr[7,"AIC"]<-mod_list_pres1[[4]]$aic

mod_tab_habstr[8,"Model_Name"]<-"Texture (Dissimilarity)"
mod_tab_habstr[8,"Function"]<-deparse1(mod_list_surv1[[4]]$formula)
mod_tab_habstr[8,"AIC"]<-mod_list_surv1[[4]]$aic

mod_tab_habstr[9,"Model_Name"]<-"Texture (Similarity)"
mod_tab_habstr[9,"Function"]<-deparse1(mod_list_pres1[[5]]$formula)
mod_tab_habstr[9,"AIC"]<-mod_list_pres1[[5]]$aic

mod_tab_habstr[10,"Model_Name"]<-"Texture (Similarity)"
mod_tab_habstr[10,"Function"]<-deparse1(mod_list_surv1[[5]]$formula)
mod_tab_habstr[10,"AIC"]<-mod_list_surv1[[5]]$aic

mod_tab_habstr[11,"Model_Name"]<-"Tidal Restrictions"
mod_tab_habstr[11,"Function"]<-deparse1(mod_list_pres1[[6]]$formula)
mod_tab_habstr[11,"AIC"]<-mod_list_pres1[[6]]$aic

mod_tab_habstr[12,"Model_Name"]<-"Tidal Restrictions"
mod_tab_habstr[12,"Function"]<-deparse1(mod_list_surv1[[6]]$formula)
mod_tab_habstr[12,"AIC"]<-mod_list_surv1[[6]]$aic

mod_tab_habstr[13,"Model_Name"]<-"High Marsh"
mod_tab_habstr[13,"Function"]<-deparse1(mod_list_pres1[[7]]$formula)
mod_tab_habstr[13,"AIC"]<-mod_list_pres1[[7]]$aic

mod_tab_habstr[14,"Model_Name"]<-"High Marsh"
mod_tab_habstr[14,"Function"]<-deparse1(mod_list_surv1[[7]]$formula)
mod_tab_habstr[14,"AIC"]<-mod_list_surv1[[7]]$aic

mod_tab_habstr[15,"Model_Name"]<-"High Marsh + NDVI"
mod_tab_habstr[15,"Function"]<-deparse1(mod_list_pres1[[8]]$formula)
mod_tab_habstr[15,"AIC"]<-mod_list_pres1[[8]]$aic

mod_tab_habstr[16,"Model_Name"]<-"High Marsh + NDVI"
mod_tab_habstr[16,"Function"]<-deparse1(mod_list_surv1[[8]]$formula)
mod_tab_habstr[16,"AIC"]<-mod_list_surv1[[8]]$aic

mod_tab_habstr[17,"Model_Name"]<-"High Marsh + UVVR"
mod_tab_habstr[17,"Function"]<-deparse1(mod_list_pres1[[9]]$formula)
mod_tab_habstr[17,"AIC"]<-mod_list_pres1[[9]]$aic

mod_tab_habstr[18,"Model_Name"]<-"High Marsh + UVVR"
mod_tab_habstr[18,"Function"]<-deparse1(mod_list_surv1[[9]]$formula)
mod_tab_habstr[18,"AIC"]<-mod_list_surv1[[9]]$aic


mod_tab_habstr[19,"Model_Name"]<-"High Marsh + PCA"
mod_tab_habstr[19,"Function"]<-deparse1(mod_list_pres1[[10]]$formula)
mod_tab_habstr[19,"AIC"]<-mod_list_pres1[[10]]$aic

mod_tab_habstr[20,"Model_Name"]<-"High Marsh + PCA"
mod_tab_habstr[20,"Function"]<-deparse1(mod_list_surv1[[10]]$formula)
mod_tab_habstr[20,"AIC"]<-mod_list_surv1[[10]]$aic


mod_tab_habstr[21,"Model_Name"]<-"High Marsh + Texture (Dissimilarity)"
mod_tab_habstr[21,"Function"]<-deparse1(mod_list_pres1[[11]]$formula)
mod_tab_habstr[21,"AIC"]<-mod_list_pres1[[11]]$aic

mod_tab_habstr[22,"Model_Name"]<-"High Marsh + Texture (Dissimilarity)"
mod_tab_habstr[22,"Function"]<-deparse1(mod_list_surv1[[11]]$formula)
mod_tab_habstr[22,"AIC"]<-mod_list_surv1[[11]]$aic

mod_tab_habstr[23,"Model_Name"]<-"High Marsh + Texture (Similarity)"
mod_tab_habstr[23,"Function"]<-deparse1(mod_list_pres1[[12]]$formula)
mod_tab_habstr[23,"AIC"]<-mod_list_pres1[[12]]$aic

mod_tab_habstr[24,"Model_Name"]<-"High Marsh + Texture (Similarity)"
mod_tab_habstr[24,"Function"]<-deparse1(mod_list_surv1[[12]]$formula)
mod_tab_habstr[24,"AIC"]<-mod_list_surv1[[12]]$aic

mod_tab_habstr[25,"Model_Name"]<-"High Marsh + Tidal Restrictions"
mod_tab_habstr[25,"Function"]<-deparse1(mod_list_pres1[[13]]$formula)
mod_tab_habstr[25,"AIC"]<-mod_list_pres1[[13]]$aic

mod_tab_habstr[26,"Model_Name"]<-"High Marsh + Tidal Restrictions"
mod_tab_habstr[26,"Function"]<-deparse1(mod_list_surv1[[13]]$formula)
mod_tab_habstr[26,"AIC"]<-mod_list_surv1[[13]]$aic

mod_tab_habstr[27,"Model_Name"]<-"Full Additional Habitat"
mod_tab_habstr[27,"Function"]<-deparse1(mod_list_pres1[[14]]$formula)
mod_tab_habstr[27,"AIC"]<-mod_list_pres1[[14]]$aic

mod_tab_habstr[28,"Model_Name"]<-"Full Additional Habitat"
mod_tab_habstr[28,"Function"]<-deparse1(mod_list_surv1[[14]]$formula)
mod_tab_habstr[28,"AIC"]<-mod_list_surv1[[14]]$aic


mod_tab_habstr<-group_by(mod_tab_habstr,Response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(Response,Model_Comparison,dAIC)





# B. Quality within high marsh (interactions) - sparrow presence does not mean nest presence
  # within high marsh habitat, certain habitat qualities are more suitable for nesting
  # eg healthier vegetation (NDVI) grows taller and denser
  # compare highmarsh only, ndvi only, highmarsh + ndvi to high marsh*ndvi
mod_list_pres2<-list()
mod_list_surv2<-list()

#do the interaction models do better than the additive models?


# high marsh and PC reflect
mod_list_pres2[[1]]<-glm(y~HIMARSH*pca, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv2[[1]]<-glm(y~HIMARSH*pca, 
                        data=surv_dat,
                        family = binomial(link="logit"))


#Do adding the interactions improve all habitat models?
mod_list_pres2[[2]]<-glm(y~uvvr_mean+uvvr_diff+cor_txt+tideres+ent_txt+HIMARSH*pca, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv2[[2]]<-glm(y~uvvr_diff+uvvr_mean+ent_txt+tideres+cor_txt+HIMARSH*pca,
                         data=surv_dat,
                         family = binomial(link="logit"))


# store model selection factors in a table
mod_tab_int<-data.frame(Response=rep(c("Presence","Success"),length(mod_list_pres2)),
                           Model_Comparison=rep("High Marsh Quality",2*length(mod_list_pres2)),
                           Model_Name=rep(NA,2*length(mod_list_pres2)),
                           Function=rep(NA,2*length(mod_list_pres2)),
                           AIC=rep(NA,2*length(mod_list_pres2)),
                           dAIC=rep(NA,2*length(mod_list_pres2)))

mod_tab_int[1,"Model_Name"]<-"High Marsh * Raw Reflectance"
mod_tab_int[1,"Function"]<-deparse1(mod_list_pres2[[1]]$formula)
mod_tab_int[1,"AIC"]<-mod_list_pres2[[1]]$aic

mod_tab_int[2,"Model_Name"]<-"High Marsh * Raw Reflectance"
mod_tab_int[2,"Function"]<-deparse1(mod_list_surv2[[1]]$formula)
mod_tab_int[2,"AIC"]<-mod_list_surv2[[1]]$aic

mod_tab_int[3,"Model_Name"]<-"All Habitat with Interaction"
mod_tab_int[3,"Function"]<-deparse1(mod_list_pres2[[2]]$formula)
mod_tab_int[3,"AIC"]<-mod_list_pres2[[2]]$aic

mod_tab_int[4,"Model_Name"]<-"All Habitat with Interaction"
mod_tab_int[4,"Function"]<-deparse1(mod_list_surv2[[2]]$formula)
mod_tab_int[4,"AIC"]<-mod_list_surv2[[2]]$aic


mod_tab_int<-group_by(mod_tab_int,Response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(Response,dAIC)




# compare models
mod_tab<-rbind(mod_tab_habstr,mod_tab_int)%>%
  group_by(Response)%>%
  mutate(ΔAIC_within_Model_Comparison=round(dAIC,1),
         Overall_ΔAIC_across_Comparisons=round(AIC-min(AIC),1))%>%
  arrange(Response,Model_Comparison,dAIC)%>%
  ungroup()%>%
  dplyr::select(-Function,-dAIC)

write.csv(mod_tab,paste0(path_out,"Final_outputs/Model_Results/model_selection_table_",ab_type,".csv"), row.names = F)


# select the top ranked model based on AIC (any models within 2 delta AIC)
mod_pres<-mod_tab[mod_tab$Response=="Presence",]
form_pres<-mod_pres[mod_pres$AIC==min(mod_pres$AIC),]$fun

mod_surv<-mod_tab[mod_tab$Response=="Success",]
form_surv<-mod_surv[mod_surv$AIC==min(mod_surv$AIC),]$fun

# Use step function to find best model accounting for all combinations of variables
p.mod<-step(glm(as.formula(paste0("y~",paste(all_terms[which(!(all_terms%in%c("HIMARSH","pca")))],collapse = "+"),"+HIMARSH*pca")), pres_dat,family=binomial(link = "logit")))
# for pres, keeps all
s.mod<-step(glm(as.formula(paste0("y~",paste(all_terms[which(!(all_terms%in%c("HIMARSH","pca")))],collapse = "+"),"+HIMARSH*pca")), surv_dat,family=binomial(link = "logit")))
# for surv, removes ent_txt, uvvr_mean, ndvi - use cor and diff instead for these
#Use proportion of high marsh instead of high marsh at nest center - error in nest location and resolution uncertainty in veg layer



## 2. Predict to test data
#--------------------------------------

d.pres.glm<-list()
d.surv.glm<-list()

for (i in 1:k) {
  train <- pres_dat[pres_dat$group != i,]
  test <- pres_dat[pres_dat$group == i,]
  #mod.p <- glm(form_pres, data = train, family = binomial(link = "logit"))
  d.pres.glm[[i]] <- data.frame(id=test$id,
                                obs=test$y, 
                                pred=predict(mod.p,test, type="response"))
}

for (i in 1:k) {
  train <- surv_dat[surv_dat$group != i,]
  test <- surv_dat[surv_dat$group == i,]
  #mod.s <- glm(form_surv, data = train, family = binomial(link = "logit"))
  d.surv.glm[[i]] <- data.frame(id=test$id,
                                obs=test$y, 
                                pred=predict(mod.s,test, type="response"))
}


if(predict.surf==T){
## 3. Predict to Spatial Surface
#-------------------------------------------

# a) list predictor variables in the selected models
terms_pres<-unlist(strsplit(as.character(substring(form_pres,5)),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()%>%
  unique()

terms_pres<-unlist(strsplit(as.character(substring(form_surv,5)),split=" * ", fixed=TRUE))%>%
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
  names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","tideres","HIMARSH")

  # c) select just the layers that were used as predictor variables in the model
  mod_preds_p<-predictors[[names(predictors)%in%terms_pres]]
  mod_preds_s<-predictors[[names(predictors)%in%terms_surv]]
  
  # d) set areas outside marsh to NA (mask all predictors)
  mask<-predictors["Highmarsh"]
  mask[mask==0|mask==9|mask==7|mask==8]<-NA
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
  
  
  #Write predicted values to csv
  
  
  pres_out<-rast(list(preds_mask_p2,glm_predict_pres[[i]]))
  names(pres_out[[nlyr(pres_out)]])<-"predictions"
  mat_p[[i]]<-as.data.frame(terra::as.matrix(pres_out,wide=F))%>%
    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
  
  surv_out<-rast(list(preds_mask_s2,glm_predict_surv[[i]]))
  names(surv_out[[nlyr(surv_out)]])<-"predictions"
  mat_s[[i]]<-as.data.frame(terra::as.matrix(surv_out,wide=F))%>%
    mutate(id=paste0(seq(1,nrow(.),1),"z",i))
  
  
  
  
  #zone 1 is too big, process it in parts
  if(i==1){
    
    n_div<-4
    dat_zone_list<-list()
    row_start<-c()
    row_end<-c()
    
    
    pres_out<-rast(list(preds_mask_p2,glm_predict_pres[[i]]))
    names(pres_out[[nlyr(pres_out)]])<-"predictions"
    
    
    surv_out<-rast(list(preds_mask_s2,glm_predict_surv[[i]]))
    names(surv_out[[nlyr(surv_out)]])<-"predictions"
    
    #divide zone into 4
    for(j in 1:n_div){
      row_start_p[j]<-((round(nrow(pres_out)/n_div))*(j-1))+1
      row_end_p[j]<-(round(nrow(pres_out)/n_div))*j
      pres_out_sub<-pres_out[c(row_start_p[j]:row_end_p[j]),c(1:ncol(pres_out)),drop=F]
      
      row_start_s[j]<-((round(nrow(surv_out)/n_div))*(j-1))+1
      row_end_s[j]<-(round(nrow(surv_out)/n_div))*j
      surv_out_sub<-surv_out[c(row_start_s[j]:row_end_s[j]),c(1:ncol(surv_out)),drop=F]
      
      #coerce raster stack into a matrix with datasets (layers) as columns and cells as rows (wide = F will do this)
      mat_p_z1[[j]]<-as.data.frame(terra::as.matrix(pres_out_sub,wide=F))%>%
        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
      
      mat_s_z1[[j]]<-as.data.frame(terra::as.matrix(surv_out_sub,wide=F))%>%
        mutate(id=paste0(seq(1,nrow(.),1),"z",i))
    }
    
    
    #bind all the zone sections into one dataframe
    mat_p[[1]]<-do.call("rbind",mat_p_z1)
    mat_s[[1]]<-do.call("rbind",mat_s_z1)
  }
  
  
  #write each zone prediction dataset to matrix file
  write.csv(mat_p[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Placement/Z",i,"_predictionGLM_30m_",ab_type,"_placement.csv"),row.names = F)
  write.csv(mat_s[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Success/Z",i,"_predictionGLM_30m_",ab_type,"_success.csv"),row.names = F)
  
  #Write rasters
  writeRaster(glm_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_predsGLM_30m",ab_type,".tif"),overwrite=T)
  writeRaster(glm_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_predsGLM_30m",ab_type,".tif"),overwrite=T)
}
}