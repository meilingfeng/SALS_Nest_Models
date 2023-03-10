
library(dismo)
library(tidyverse)
library(lme4) #mixed effects models
library(terra)

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. parameters
  # what prediction resolution?
reso<-30

  # use random background points or veg plots as absences?
ab_type<-"b"
#ab_type<-"v" 

  # pick which predictors and response to use for models
predictors <- "uvvr"
#predictors <- "no uvvr" #removes observations mising uvvr data
analysis<-"presence" #for nest presence
#analysis<-"survival" #for nest fates

  # veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))

## 3. Load observations with predictors
  # list of predictor files for each zone as a raster stack
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))

  # point predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_local.csv"))%>%
  filter(bp%in%c("p",ab_type))%>%
  mutate(presence=ifelse(bp=="p",1,0))%>%
  left_join(veg_codes,by="veg_class")

  # buffered predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"))%>%
  dplyr::select(id,HIMARSH,LOMARSH,POOL,PHRG,STRM,MUD,UPLND,TERRBRD)%>%
  right_join(dat,by="id")%>%
  mutate(veg_code=as.factor(veg_code))



## 4. Remove missing values from predictors
#data with all vegetation class and NAIP predictors
dat_comp<-dat[complete.cases(dat[ , c("veg_code","ent_txt")]),]
dat_comp[,c("HIMARSH","LOMARSH","PHRG","POOL","MUD","STRM","UPLND","TERRBRD")]<-dat_comp[,c("HIMARSH","LOMARSH","PHRG","POOL","MUD","STRM","UPLND","TERRBRD")]%>%
  replace(is.na(.),0)
#associated predictors
form<-y~ndvi+cor_txt+ent_txt+veg_code+HIMARSH+LOMARSH+precip+pca #+PHRG+POOL+MUD+STRM+UPLND+TERRBRD

#data with all predictors including UVVR
dat_comp_all<-dat_comp[complete.cases(dat_comp[ , c("uvvr_mean","uvvr_diff")]),]
#associated predictors
form_all<-y~ndvi+cor_txt+ent_txt+veg_class+HIMARSH+LOMARSH+PHRG+POOL+MUD+STRM+UPLND+TERRBRD+precip+pca+uvvr_mean+uvvr_diff

## 5. pick which predictors to use for models
if(predictors=="no uvvr"){
  mod_dat<-dat_comp[!is.na(dat_comp$fate),]
  mod_form<-form
}

if(predictors=="uvvr"){
  mod_dat<-dat_comp_all[!is.na(dat_comp_all$fate),]
  mod_form<-form_all
}



## 6. divide into testing and training data 80 training/20 testing
table(mod_dat$presence)
dat_split<-mod_dat %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
#check ratio of data
map_int(dat_split, nrow)
#check balance in each split dataset
table(dat_split[["train"]]$presence) #retains good balance
table(dat_split[["test"]]$presence)
#split into separate objects
train<-dat_split[["train"]]
test<-dat_split[["test"]]


## 7. pick which response to use for models
if(analysis=="survival"){
  y<-fate
}
if(predictors=="presence"){
  y<-presence
}



### Build Models
#--------------------------------------------------

## 1. GLMs with presence-background
m<-glm(mod_form, 
           data=train,
           family = binomial(link="logit"))

# How to choose how many variables to include?
# Stepwise simplification with Akaike Information Criterion
# penalized for complexity, looks for best fit for least complex model
aic_form<-step(m)[["formula"]] #this goes through each combination of variables to find the best fit


glm_final <- glm(formula = aic_form, data = dat_split[["train"]])
summary(glm_final)


#plogis() #converts log likelihood to probability of success




## 2. Bioclim with presence-only

bc <- bioclim(dat_comp[dat_comp$presence==1,c('veg_class', 'ndvi', "precip","HIMARSH","LOMARSH","PHRG","MUD","UPLND","TERRBRD","STRM")])

pairs(bc)






### Make Predictions
#----------------------------------
#load rasters
predictors<-rast(unlist(file_list_all_zones[[3]]))
#select variables in the model
aic_form
names(predictors)
names(predictors)<-c("veg_code","cor_txt","ent_txt","ndvi","pca","uvvr_diff","uvvr_mean","precip","HIMARSH","LOMARSH")
predictors<-predictors[[c(1,3:5,8:10)]]

#give each template cell a unique id
#ids<-list()
#for(i in 1:length(file_list_all_zones)){
#  r<-rast(file_list_all_zones[[i]][[1]])
#  r[r==0]<-NA
#  id<-mask(rast(extent=ext(r),crs=crs(r), resolution=res(r),vals="a"),r)
#  id[id=="a"]<-paste0(c(1:ncell(id[id=="a"])),"z",i)

#  ids[[i]]<-id}
#plot(id)
  
#need to make veg_class raster a factor
t<-as.data.frame(predictors[[1]])
t$veg_code<-as.factor(t$veg_code)
r<-rast(crs=crs(predictors[[1]]),extent=ext(predictors[[1]]),res=res(predictors[[1]]),vals=t$veg_code)
predictors[[1]]<-r

#set areas outside marsh to NA
mask<-predictors[[1]]
mask[mask=="0"]<-NA
predictors_mask<-mask(predictors,mask)

p_glm<-predict(predictors_mask,glm_final)

plot(p2)




### Model Evaluation 
#--------------------------------------
#Comparison to training data
#for GLM, can look at how much deviance is explained, whether there are patterns in residuals, whether there are points with high leverage

#does the model makes sense ecologically?
#do fitted functions (shape of modeled relationships) make sense?
#are there spatial patterns in model residuals?

#Correlation and AUC
p <- dat_comp_split[["train"]][dat_comp_split[["train"]]$presence==1,]
a <- dat_comp_split[["train"]][dat_comp_split[["train"]]$presence==0,]
evaluate(p,a,glm_final)
par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
       pch=c(21,24), col=c('red', 'blue'))
comb <- c(p,a)
group <- c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~group, col=c('blue', 'red'))