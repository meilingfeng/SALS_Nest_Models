library(tidyverse)
library(car)#vif
library(mgcv)
library(dismo)
library(sf)


### Set up
# -------------------------------------------
#source("05c_Independence_Interactions.R")
source("05b_candidate_models.R")  
  
## Model Evaluation
#---------------------------------------------------

#remove the full model that has collinearity
#mod_list_pres<-mod_list_pres[-1]
#mod_list_surv<-mod_list_surv[-1]


#Spatial sorting bias (SSB) - background points that cover a larger spatial extent have higher AUC
#SSB = 1 is no bias, close to zero is a lot of bias


pres<-filter(pres_dat,presence==1)%>%
  dplyr::select(longitude,latitude)
backgr<-filter(pres_dat,presence==0)%>%
  dplyr::select(longitude,latitude)
set.seed(0)
nr <- nrow(pres)
s <- sample(nr, 0.25 * nr)
p_train <- pres[-s, ]
p_test <- pres[s, ]
nr <- nrow(backgr)
set.seed(9)
s <- sample(nr, 0.25 * nr)
b_train <- backgr[-s, ]
b_test <- backgr[s, ]

sb <- ssb(p_test, b_test, p_train)
sb[,1] / sb[,2]

#There is a lot of bias. Veg points improve slightly but not much.

#i <- pwdSample(p_test, b_test, p_train, n=1, tr=0.1)
#p_test_pwd <- p_test[!is.na(i[,1]), ]
#b_test_pwd <- b_test[na.omit(as.vector(i)), ]
#sb_pwd <- ssb(p_test_pwd, b_test_pwd, p_train)
#sb[,1] / sb[,2]





# AUC - evalutate model performance using training data
for(i in 1:length(mod_list_pres_final)){
mod_tab[mod_tab$response=="Presence" & mod_tab$fun == deparse1(mod_list_pres_final[[i]]$formula), "AUC"]<- evaluate(filter(pres_test,y==1), filter(pres_test,y==0), mod_list_pres_final[[i]])@auc
}

for(i in 1:length(mod_list_surv_final)){
  mod_tab[mod_tab$response=="Success" & mod_tab$fun == deparse1(mod_list_surv_final[[i]]$formula), "AUC"]<- evaluate(filter(surv_test,y==1), filter(surv_test,y==0), mod_list_surv_final[[i]])@auc
}



#write.csv(mod_tab,paste0(path_out,"Final_outputs/candidate_model_selection_table_",ab_type,".csv"),row.names = F)


# select the top ranked model based on AIC and AUC
mod_pres<-mod_tab[mod_tab$response=="Presence",]
form_pres<-mod_pres[mod_pres$AIC<min(mod_pres$AIC)+2 & mod_pres$AUC > max(mod_pres$AUC)-2,]$fun
final_mod_pres<-glm(form_pres, data = pres_train, family = binomial(link = "logit"))

mod_surv<-mod_tab[mod_tab$response=="Success",]
form_surv<-mod_surv[mod_surv$AIC<min(mod_surv$AIC)+1 & mod_surv$AUC > max(mod_surv$AUC)-2,]$fun
final_mod_surv<-glm(form_surv, data = surv_train, family = binomial(link = "logit"))

final_mod_pres$fitted.values
#plot evaluation
e<-evaluate(filter(pres_test,y==1), filter(pres_test,y==0),final_mod_pres)
png(filename=paste0(path_out,"Final_outputs/pres_model_evaluation_",ab_type,".png"))
par(mfrow=c(1, 2))
density(e)
boxplot(e, col=c('blue', 'red'))
dev.off()


e<-evaluate(filter(surv_test,y==1), filter(surv_test,y==0),final_mod_surv)
png(paste0(path_out,"Final_outputs/surv_model_evaluation_",ab_type,".png"))
par(mfrow=c(1, 2))
density(e)
boxplot(e, col=c('blue', 'red'))
dev.off()
