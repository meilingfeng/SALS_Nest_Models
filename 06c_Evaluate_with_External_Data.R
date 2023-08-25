library(tidyverse)
library(dismo)
library(gbm)

## Set up 
#------------------------------------------------------------------
#Load data 
source("05a_Set_Model_Parameters.R")

# Create final models, don't do k-folds
predict.surf<-T
build<-F

# Run only with GLM and the better ML method
source("05c_boosted_regression_trees.R")
source("05b_GLM_Selection_Prediction.R")


model.names<-c("BRTs","GLM")



#Predict to validation datasets
#----------------------------------------------------
#GLMs
train <- pres_dat[pres_dat$group != i,]
test <- pres_dat[pres_dat$group == i,]
mod.p <- glm(form_pres, data = train, family = binomial(link = "logit"))
d.pres.glm[[i]] <- data.frame(id=test$id,
                              obs=test$y, 
                              pred=predict(mod.p,test, type="response"))

train <- surv_dat[surv_dat$group != i,]
test <- surv_dat[surv_dat$group == i,]
mod.s <- glm(form_surv, data = train, family = binomial(link = "logit"))
d.surv.glm[[i]] <- data.frame(id=test$id,
                              obs=test$y, 
                              pred=predict(mod.s,test, type="response"))


#BRTs


#combine each model's predictions into one df
d.pres<-cbind(d.pres.brt[[i]],d.pres.mxt[[i]][,3],d.pres.glm[[i]][,3])%>%
    rename("Maxent"="d.pres.mxt[[i]][, 3]","GLM"="d.pres.glm[[i]][, 3]","BRTs"="pred")
  
d.surv<-cbind(d.surv.brt[[i]],d.surv.mxt[[i]][,3],d.surv.glm[[i]][,3])%>%
    rename("Maxent"="d.surv.mxt[[i]][, 3]","GLM"="d.surv.glm[[i]][, 3]","BRTs"="pred")



## Get Nest presence and nest survival Discrimination and Calibration metrics for each modeling method
#-----------------------------------------------------------------------------------------------------------
#Use AUC and TSS since they respond better to balanced data
#Use Kappa over Brier for calibration since it also responds better to balanced data
#Thinning to balance improves model performance, especially if discrimination or binary calibration is the goal
#thinning just the majority also helps retain all presence records for parametric methods.
#objects to hold model performance of each fold
thr.pres<-list()
thr.surv<-list()
thr.p<-c()
thr.s<-c()

accu.pres<-list()
accu.surv<-list()

# calculate metrics for each fold
for(i in 1:k){
  # get optimal threshold for Maximizing Kappa
  thr.p<-optimal.thresholds(d.pres[[i]],opt.methods = c("MaxKappa","MaxSens+Spec"))
  thr.pres[[i]]<-c(thr.p$BRTs,thr.p$Maxent,thr.p$GLM)
  
  thr.s<-optimal.thresholds(d.surv[[i]],opt.methods = c("MaxKappa","MaxSens+Spec"))
  thr.surv[[i]]<-c(thr.s$BRTs,thr.s$Maxent,thr.s$GLM)
  
  
  # accuracy metric tables
  accu.p<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]][c(2,4,6)])
  accu.p$Kappa<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]][c(1,3,5)])$Kappa
  accu.p$Kappa.sd<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]][c(1,3,5)])$Kappa.sd
  accu.p[, -c(1, 2)] <- signif(accu.p[, -c(1, 2)], digits = 3)
  accu.pres[[i]]<-accu.p
  
  accu.s<-presence.absence.accuracy(d.surv[[i]],threshold=thr.surv[[i]][c(2,4,6)])
  accu.s$Kappa<-presence.absence.accuracy(d.surv[[i]],threshold = thr.surv[[i]][c(1,3,5)])$Kappa
  accu.s$Kappa.sd<-presence.absence.accuracy(d.surv[[i]],threshold = thr.surv[[i]][c(1,3,5)])$Kappa.sd
  accu.s[, -c(1, 2)] <- signif(accu.s[, -c(1, 2)], digits = 3)
  accu.surv[[i]]<-accu.s
}


# Then calculate the mean and SD of each metric for each modeling method
auc.p<-list() #AUC is a continuous measure of discrimination (0.8-1.0 is really good)
kappa.p<-list() #Kappa is a binary measure of calibration (above 0.4 is good)
tss.p<-list() # True SKill Statistic is a binary measure of discrimination
thres.p<-list() #use a threshold  that maximizes Kappa (Accuracy)

auc.s<-list()
kappa.s<-list()
tss.s<-list()
thres.s<-list()

for(i in 1:k){
  auc.p[[i]]<-accu.pres[[i]][,"AUC"]
  kappa.p[[i]]<-accu.pres[[i]][,"Kappa"]
  tss.p[[i]]<-accu.pres[[i]][,"sensitivity"]+accu.pres[[i]][,"specificity"]-1
  thres.p[[i]]<-accu.pres[[i]][,"threshold"]
  
  auc.s[[i]]<-accu.surv[[i]][,"AUC"]
  kappa.s[[i]]<-accu.surv[[i]][,"Kappa"]
  tss.s[[i]]<-accu.surv[[i]][,"sensitivity"]+accu.surv[[i]][1,"specificity"]-1
  thres.s[[i]]<-accu.surv[[i]][,"threshold"]
}


eval.tab<-data.frame(model=model.names,
                     auc.p=unlist(pmap(auc.p,mean)),
                     auc.p.sd=rep(NA,length(model.names)),
                     auc.s=unlist(pmap(auc.s,mean)),
                     auc.s.sd=rep(NA,length(model.names)),
                     tss.p=unlist(pmap(tss.p,mean)),
                     tss.p.sd=rep(NA,length(model.names)),
                     tss.s=unlist(pmap(tss.s,mean)),
                     tss.s.sd=rep(NA,length(model.names)),
                     kappa.p=unlist(pmap(kappa.p,mean)),
                     kappa.p.sd=rep(NA,length(model.names)),
                     kappa.s=unlist(pmap(kappa.s,mean)),
                     kappa.s.sd=rep(NA,length(model.names)))

for (i in 1:length(model.names)){ #for each model
  temp.auc.p<-c()
  temp.auc.s<-c()
  temp.tss.p<-c()
  temp.tss.s<-c()
  temp.kappa.p<-c()
  temp.kappa.s<-c()
  for(j in 1:k){#get all the folds
    temp.auc.p[j]<-auc.p[[j]][[i]]
    temp.auc.s[j]<-auc.s[[j]][[i]]
    temp.tss.p[j]<-tss.p[[j]][[i]]
    temp.tss.s[j]<-tss.s[[j]][[i]]
    temp.kappa.p[j]<-kappa.p[[j]][[i]]
    temp.kappa.s[j]<-kappa.s[[j]][[i]]
  }
  eval.tab[i, "auc.p.sd"]<-round(sd(temp.auc.p,na.rm=T),3)
  eval.tab[i, "auc.s.sd"]<-round(sd(temp.auc.s,na.rm=T),3)
  eval.tab[i, "tss.p.sd"]<-round(sd(temp.tss.p,na.rm=T),3)
  eval.tab[i, "tss.s.sd"]<-round(sd(temp.tss.s,na.rm=T),3)
  eval.tab[i, "kappa.p.sd"]<-round(sd(temp.kappa.p,na.rm=T),3)
  eval.tab[i, "kappa.s.sd"]<-round(sd(temp.kappa.s,na.rm=T),3) # alternatively use SE ", 95% CI ",round(mean(auc_p)+1.96*(sd(auc_p)/length(auc_p)),3),"-",round(mean(auc_p)-1.96*(sd(auc_p)/length(auc_p)),3))
  
}

write.csv(eval.tab,paste0(path_out,"Final_outputs/Model_Results/model_evaluation_table_",ab_type,".csv"), row.names = F)

