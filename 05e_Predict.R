
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
#terms_surv<-terms_surv[-1]

#empty list to hold prediction surfaces
glm_predict_pres<-list()
glm_predict_surv<-list()

i<-3

## 2. loop predictions through each zone
######
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
mask<-mod_preds_p["Highmarsh"]
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
#convert from raster stack to dataframe
#preds_df<-as.data.frame(preds_mask2, 
#                        xy=T, #inlcude cell coordinates
#                        cells=T, #include cell numbers
#                        na.rm=NA #remove cells with NA in all layers (outside marsh)
#                         )%>%
#  rename(long=x,lat=y)
#keep order in model function
#preds_mask2<-preds_mask2[[order(terms_pres)]]
#predict to prediction surface
#test<-predict.glm(final_mod_pres,newdata=preds_df,type="response")
glm_predict_pres[[i]]<-predict(preds_mask_p2,final_mod_pres, type="response")
glm_predict_surv[[i]]<-predict(preds_mask_s2,final_mod_surv, type="response")
}


#Write rasters
for(i in 1:length(glm_predict_pres)){
writeRaster(glm_predict_pres[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_preds_30m",ab_type,".tif"),overwrite=T)
writeRaster(glm_predict_surv[[i]],paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_preds_30m",ab_type,".tif"),overwrite=T)
}


### Model Evaluation 
#--------------------------------------
#look at mapped predictions
for(i in 1:length(glm_predict_pres)){
ggplot()+
    geom_spatraster(data=glm_predict_pres[[i]],aes(fill=lyr1))+
    scale_fill_hypso_c(direction=-1)+
    theme_classic()
ggsave(filename = paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_pres_preds_30m",ab_type,".png"),
      # width=8,height = 8,
       dpi=300, units="in")

ggplot()+
  geom_spatraster(data=glm_predict_surv[[i]],aes(fill=lyr1))+
  scale_fill_hypso_c(direction=-1)+
  theme_classic()
ggsave(filename = paste0(path_out,"Final_outputs/Nest_Predictions/Success/z",i,"_surv_preds_30m",ab_type,".png"),
       # width=8,height = 8,
       dpi=300, units="in")
}

test<-rast(paste0(dat_path,"Correll_Marsh_Zones/Zone4_DEM_30align.tif"))
ggplot()+
  geom_spatraster(data=test,aes(fill=Zone4_DEM_30align))+
  scale_fill_hypso_b(breaks=c(1:8),direction=-1)+
  theme_classic()
ggsave(filename = paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z4_veg_30m.png"),
       # width=8,height = 8,
       dpi=300, units="in")
#Comparison to training data
#for GLM, can look at how much deviance is explained, whether there are patterns in residuals, whether there are points with high leverage

#does the model makes sense ecologically?
#do fitted functions (shape of modeled relationships) make sense?
#are there spatial patterns in model residuals?


