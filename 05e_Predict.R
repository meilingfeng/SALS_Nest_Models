
library(dismo)
library(tidyverse)
library(tidyterra)
library(lme4) #mixed effects models
library(terra)
#install.packages('terra', repos='https://rspatial.r-universe.dev') #development version

#**** NEED TO RESCALE VARABLES BEFORE SAMPLING AT NESTS

### Set up
# -------------------------------------------
source("05d_Model_Selection_Evaluation.R")


### Load Final models
#--------------------------------------------------
final_mod_pres
final_mod_surv

#final_mod_pres<-glm(y~ndvi+cor_txt+veg_code+HIMARSH+ent_txt+pca+uvvr_mean, 
#                    data=pres_train,
#                    family = binomial(link="logit"))

#final_mod_surv<-glm(y~ndvi+cor_txt+veg_code+HIMARSH+ent_txt+pca+uvvr_mean, 
#    data=surv_train,
#    family = binomial(link="logit"))

#for now just use unscaled lat and pca, remove lat- why is it messing up the predictions?


### Make Predictions ***modify to include separate terms for nest survival model - how do I incorporate site as random effect when extrapolating?
#----------------------------------
## 1. list predictor variables in the selected model
######
terms_pres<-unlist(strsplit(as.character(final_mod_pres$terms[[3]]),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()
#get rid of operators and interaction terms
terms_pres<-terms_pres[-1] #remove lat (9) for now

terms_surv<-unlist(strsplit(as.character(final_mod_surv$terms[[3]]),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()
#get rid of operators and interaction terms
terms_surv<-terms_surv[-1]

#empty list to hold prediction surfaces
glm_predict_pres<-list()
glm_predict_surv<-list()


mat_p<-list()
mat_s<-list()
mat_p_z1<-list()
mat_s_z1<-list()

## 2. loop predictions through each zone
######
i<-3
for (i in 1:length(file_list_all_zones)){
  
# a) load raster predictor layers for a particular zone
  #patch fix **** use zeroed NDVI
temp<-unlist(file_list_all_zones[[i]])
temp[4]<-sub("NDVI","zeroed_NDVI",temp[4])
temp[4]<-sub("Data","Outputs/Final_outputs",temp[4])
temp[4]<-sub("NAIP/","NAIP/Z",temp[4])
  #end patch fix
predictors<-rast(unlist(temp))

# b) name layers as their variables (rename veg_code as just Highmarsh since we're only using that one class for now)
names(predictors)<-c("Highmarsh","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","HIMARSH","LOMARSH")
#names(predictors)<-c("veg_code","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","HIMARSH","LOMARSH")

# c) select just the layers that were used as predictor variables in the model
mod_preds_p<-predictors[[names(predictors)%in%terms_pres]]
mod_preds_s<-predictors[[names(predictors)%in%terms_surv]]


# d) format predictor rasters

# get latitude 
if("latitude"%in%terms_pres){
lat<-rast(crs=crs(mod_preds_p[[1]]),extent=ext(mod_preds_p[[1]]),resolution=res(mod_preds_p[[1]]),vals=xyFromCell(terra::project(mod_preds_p[[1]],"EPSG:4617"),1:ncell(mod_preds_p[[1]]))[,2])
#lat<-rast(crs=crs(mod_preds[[1]]),extent=ext(mod_preds[[1]]),resolution=res(mod_preds[[1]]),vals=scale(xyFromCell(mod_preds[[1]],1:ncell(mod_preds[[1]]))[,2],center=F))
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

#if presence model terms include Highmarsh (or veg class)
if("Highmarsh"%in%terms_pres){
# make veg_class raster a factor that just indicates high marsh
t<-as.data.frame(mod_preds_p[["Highmarsh"]])%>%
  mutate(Highmarsh=as.factor(case_when(Highmarsh==1~ "1", #set high marsh to 1 (presence)
                                       Highmarsh%in%c(0,7,9) ~ NA_character_, #stream and upland set to NA
                                       (!(is.na(Highmarsh))&Highmarsh!=1)~ "0")))#set all other veg classes to 0 (absence of high marsh)))

r<-rast(crs=crs(mod_preds_p[["Highmarsh"]]),extent=ext(mod_preds_p[["Highmarsh"]]),resolution=res(mod_preds_p[["Highmarsh"]]),vals=t$Highmarsh)

# add to predictor list
mod_preds_p[["Highmarsh"]]<-r
}


#if survival model terms include Highmarsh (or veg class)
if("Highmarsh"%in%terms_surv){
  # make veg_class raster a factor that just indicates high marsh
  t<-as.data.frame(mod_preds_s[["Highmarsh"]])%>%
    mutate(Highmarsh=as.factor(case_when(Highmarsh==1~ "1", #set high marsh to 1 (presence)
                                         !(Highmarsh%in%c(1,0))~ "0",#set all other veg classes to 0 (absence of high marsh)
                                         Highmarsh==0 ~ NA_character_)))
  
  r<-rast(crs=crs(mod_preds_s[["Highmarsh"]]),extent=ext(mod_preds_s[["Highmarsh"]]),resolution=res(mod_preds_s[["Highmarsh"]]),vals=t$Highmarsh)
  
#add to predictor list
mod_preds_s[["Highmarsh"]]<-r

}



# e) set areas outside marsh to NA (mask all predictors)
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
  write.csv(mat_p[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Placement/Z",i,"_prediction_30m_",ab_type,"_placement.csv"),row.names = F)
  write.csv(mat_s[[i]],file=paste0(path_out,"/Final_outputs/Nest_Predictions/Success/Z",i,"_prediction_30m_",ab_type,"_success.csv"),row.names = F)

#Write rasters
writeRaster(glm_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_preds_30m",ab_type,".tif"),overwrite=T)
writeRaster(glm_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_preds_30m",ab_type,".tif"),overwrite=T)
}




### Model Evaluation 
#--------------------------------------
#look at mapped predictions
for(i in 1:length(glm_predict_pres)){
ggplot()+
    geom_spatraster(data=glm_predict_pres[[i]],aes(fill=lyr1))+
    scale_fill_viridis_c()+
    theme_classic()
ggsave(filename = paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_preds_30m",ab_type,".png"),
      # width=8,height = 8,
       dpi=300, units="in")

ggplot()+
  geom_spatraster(data=glm_predict_surv[[i]],aes(fill=lyr1))+
  scale_fill_viridis_c()+
  theme_classic()
ggsave(filename = paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_preds_30m",ab_type,".png"),
       # width=8,height = 8,
       dpi=300, units="in")
}


#Comparison to training data
#for GLM, can look at how much deviance is explained, whether there are patterns in residuals, whether there are points with high leverage

#does the model makes sense ecologically?
#do fitted functions (shape of modeled relationships) make sense?
#are there spatial patterns in model residuals?


