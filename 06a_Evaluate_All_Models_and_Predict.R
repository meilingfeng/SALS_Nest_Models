library(tidyverse)
library(dismo)
library(gbm)
library(ggpubr)

## Set up 
#------------------------------------------------------------------

#if running all species, set sals_only to FALSE.
sals_only<-T

#number of folds. Make sure this matches 05a_Set_Model_Parameters.R
k<-5

# Don't predict to spatial surface yet
build<-T
predict.surf<-F

if(sals_only==T){
# Compare 3 modeling methods when looking at just SALS. If running all species, just do the BRTs to assess performance
source("05c_boosted_regression_trees.R")
source("05d_Maxent.R")
source("05b_GLM_ModelComparison_Prediction.R")
  
model.names<-c("BRTs","Maxent","GLM")

}else{
source("05c_boosted_regression_trees.R")  

model.names<-c("BRTs")

}


if(sals_only==T){

#combine each model's predictions into one df per fold
d.pres<-list() #list of dfs for each fold for predicting nest presence
d.surv<-list() #list of dfs for each fold for predicting nest survival

for(i in 1:k){
  set.seed(123)
  #since training and testing datasets are of different lengths depending on the method, randomly sample the same amount of predictions
n_rows<-min(nrow(d.pres.brt[[i]]),nrow(d.pres.mxt[[i]]),nrow(d.pres.glm[[i]]))
d.pres[[i]]<-cbind(d.pres.brt[[i]][sample(nrow(d.pres.brt[[i]]),n_rows,replace = F),],
                   d.pres.mxt[[i]][sample(nrow(d.pres.mxt[[i]]),n_rows,replace = F),3],
                   d.pres.glm[[i]][sample(nrow(d.pres.glm[[i]]),n_rows,replace = F),3])%>%
  rename("Maxent"="d.pres.mxt[[i]][sample(nrow(d.pres.mxt[[i]]), n_rows, replace = F), ","GLM"="d.pres.glm[[i]][sample(nrow(d.pres.glm[[i]]), n_rows, replace = F), ","BRTs"="pred")

n_rows<-min(nrow(d.surv.brt[[i]]),nrow(d.surv.mxt[[i]]),nrow(d.surv.glm[[i]]))
d.surv[[i]]<-cbind(d.surv.brt[[i]][sample(nrow(d.surv.brt[[i]]),n_rows,replace = F),],
                   d.surv.mxt[[i]][sample(nrow(d.surv.mxt[[i]]),n_rows,replace = F),3],
                   d.surv.glm[[i]][sample(nrow(d.surv.glm[[i]]),n_rows,replace = F),3])%>%
  rename("Maxent"="d.surv.mxt[[i]][sample(nrow(d.surv.mxt[[i]]), n_rows, replace = F), ","GLM"="d.surv.glm[[i]][sample(nrow(d.surv.glm[[i]]), n_rows, replace = F), ","BRTs"="pred")
}

}else{
d.pres<-lapply(d.pres.brt,function(x)
  rename(x,"BRTs"="pred"))
}

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

# calculate performance metrics for each fold
for(i in 1:k){
# get optimal threshold for Maximizing percent correctly classified (while this can have issues in unbalanced prevalence (over predicting 0 for rare species), our data is balanced)
thr.p<-optimal.thresholds(d.pres[[i]],opt.methods = c("MaxPCC"))

if(sals_only==T){
thr.pres[[i]]<-c(thr.p$BRTs,thr.p$Maxent,thr.p$GLM)
}else{
thr.pres[[i]]<-thr.p$BRTs
}

thr.s<-optimal.thresholds(d.surv[[i]],opt.methods = c("MaxPCC"))
thr.surv[[i]]<-c(thr.s$BRTs,thr.s$Maxent,thr.s$GLM)
if(sals_only==T){
thr.surv[[i]]<-c(thr.s$BRTs,thr.s$Maxent,thr.s$GLM)
}else{
thr.surv[[i]]<-thr.s$BRTs
}

# accuracy metric tables
if(sals_only==T){
accu.p<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]])
accu.p[, -c(1, 2)] <- signif(accu.p[, -c(1, 2)], digits = 3)
accu.pres[[i]]<-accu.p

accu.s<-presence.absence.accuracy(d.surv[[i]],threshold=thr.surv[[i]])

accu.s[, -c(1, 2)] <- signif(accu.s[, -c(1, 2)], digits = 3)
accu.surv[[i]]<-accu.s
}else{
  accu.p<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]])
  accu.p[, -c(1, 2)] <- signif(accu.p[, -c(1, 2)], digits = 3)#round numbers
  accu.pres[[i]]<-accu.p
  
  accu.s<-presence.absence.accuracy(d.surv[[i]],threshold=thr.surv[[i]])
  accu.s[, -c(1, 2)] <- signif(accu.s[, -c(1, 2)], digits = 3)
  accu.surv[[i]]<-accu.s  
}
}


# Then calculate the mean and SD of each metric for each modeling method
auc.p<-list() #AUC is a continuous measure of discrimination (0.8-1.0 is really good)
kappa.p<-list() #Kappa is a binary measure of calibration (above 0.4 is good)
tss.p<-list() # True SKill Statistic is a binary measure of discrimination
sen.p<-list() #sensitivity (omission error)
spec.p<-list() #specificity (commission error)
pcc.p<-list() # percent correctly classified (overall accuracy)
thres.p<-list() #use a threshold  that maximizes PCC (overall Accuracy)

auc.s<-list()
kappa.s<-list()
tss.s<-list()
sen.s<-list() #sensitivity (omission error)
spec.s<-list() #specificity (commission error)
pcc.s<-list() # percent correctly classified (overall accuracy)
thres.s<-list() #use a threshold  that maximizes PCC (overall Accuracy)

for(i in 1:k){
auc.p[[i]]<-accu.pres[[i]][,"AUC"]
kappa.p[[i]]<-accu.pres[[i]][,"Kappa"]
tss.p[[i]]<-accu.pres[[i]][,"sensitivity"]+accu.pres[[i]][,"specificity"]-1
sen.p[[i]]<-accu.pres[[i]][,"sensitivity"]
spec.p[[i]]<-accu.pres[[i]][,"specificity"]
pcc.p[[i]]<-accu.pres[[i]][,"PCC"]
thres.p[[i]]<-accu.pres[[i]][,"threshold"]

auc.s[[i]]<-accu.surv[[i]][,"AUC"]
kappa.s[[i]]<-accu.surv[[i]][,"Kappa"]
tss.s[[i]]<-accu.surv[[i]][,"sensitivity"]+accu.surv[[i]][1,"specificity"]-1
sen.s[[i]]<-accu.surv[[i]][,"sensitivity"]
spec.s[[i]]<-accu.surv[[i]][,"specificity"]
pcc.s[[i]]<-accu.surv[[i]][,"PCC"]
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
           sen.p=unlist(pmap(sen.p,mean)),
           sen.p.sd=rep(NA,length(model.names)),
           sen.s=unlist(pmap(sen.s,mean)),
           sen.s.sd=rep(NA,length(model.names)),
           spec.p=unlist(pmap(spec.p,mean)),
           spec.p.sd=rep(NA,length(model.names)),
           spec.s=unlist(pmap(spec.s,mean)),
           spec.s.sd=rep(NA,length(model.names)),
           pcc.p=unlist(pmap(pcc.p,mean)),
           pcc.p.sd=rep(NA,length(model.names)),
           pcc.s=unlist(pmap(pcc.s,mean)),
           pcc.s.sd=rep(NA,length(model.names)),
           kappa.p=unlist(pmap(kappa.p,mean)),
           kappa.p.sd=rep(NA,length(model.names)),
           kappa.s=unlist(pmap(kappa.s,mean)),
           kappa.s.sd=rep(NA,length(model.names)))

for (i in 1:length(model.names)){ #for each model
  temp.auc.p<-c()
  temp.auc.s<-c()
  temp.tss.p<-c()
  temp.tss.s<-c()
  temp.sen.p<-c()
  temp.sen.s<-c()
  temp.spec.p<-c()
  temp.spec.s<-c()
  temp.pcc.p<-c()
  temp.pcc.s<-c()
  temp.kappa.p<-c()
  temp.kappa.s<-c()
  for(j in 1:k){#get all the folds
    temp.auc.p[j]<-auc.p[[j]][[i]]
    temp.auc.s[j]<-auc.s[[j]][[i]]
    temp.tss.p[j]<-tss.p[[j]][[i]]
    temp.tss.s[j]<-tss.s[[j]][[i]]
    temp.sen.p[j]<-sen.p[[j]][[i]]
    temp.sen.s[j]<-sen.s[[j]][[i]]
    temp.spec.p[j]<-spec.p[[j]][[i]]
    temp.spec.s[j]<-spec.s[[j]][[i]]
    temp.pcc.p[j]<-pcc.p[[j]][[i]]
    temp.pcc.s[j]<-pcc.s[[j]][[i]]
    temp.kappa.p[j]<-kappa.p[[j]][[i]]
    temp.kappa.s[j]<-kappa.s[[j]][[i]]
  }
  eval.tab[i, "auc.p.sd"]<-round(sd(temp.auc.p,na.rm=T),3)
  eval.tab[i, "auc.s.sd"]<-round(sd(temp.auc.s,na.rm=T),3)
  eval.tab[i, "tss.p.sd"]<-round(sd(temp.tss.p,na.rm=T),3)
  eval.tab[i, "tss.s.sd"]<-round(sd(temp.tss.s,na.rm=T),3)
  eval.tab[i, "sen.p.sd"]<-round(sd(temp.sen.p,na.rm=T),3)
  eval.tab[i, "sen.s.sd"]<-round(sd(temp.sen.s,na.rm=T),3)
  eval.tab[i, "spec.p.sd"]<-round(sd(temp.spec.p,na.rm=T),3)
  eval.tab[i, "spec.s.sd"]<-round(sd(temp.spec.s,na.rm=T),3)
  eval.tab[i, "pcc.p.sd"]<-round(sd(temp.pcc.p,na.rm=T),3)
  eval.tab[i, "pcc.s.sd"]<-round(sd(temp.pcc.s,na.rm=T),3)
  eval.tab[i, "kappa.p.sd"]<-round(sd(temp.kappa.p,na.rm=T),3)
  eval.tab[i, "kappa.s.sd"]<-round(sd(temp.kappa.s,na.rm=T),3) # alternatively use SE ", 95% CI ",round(mean(auc_p)+1.96*(sd(auc_p)/length(auc_p)),3),"-",round(mean(auc_p)-1.96*(sd(auc_p)/length(auc_p)),3))

}

write.csv(eval.tab,paste0(path_out,"Final_outputs/Model_Results/model_evaluation_table_2_16_25.csv"), row.names = F)


# BRT has the best performance across all metrics. 
# In general worse predictions for survival. 
# Both survival and presence have better discrimination than calibration.
# high AUC ... all these metrics are influenced by prevalence. Increasing the area of background can increase the scores. Argument for veg points? try extrapolating predictor rasters to points outside their range.


## Visual of Discrimination evaluation (important if making presence absence maps) - histogram of presences to absences 
#---------------------------------------------------------------------------------------------------------------------------------
## nest presence
  # first compare across models
disc_p_brt<-ggplot(d.pres[[1]],aes(BRTs,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="",
       title=paste0("BRTs (AUC = ",eval.tab[1,"auc.p"],", PCC = ",eval.tab[1,"pcc.p"],", TSS = ",eval.tab[1,"tss.p"],")"))
disc_p_mxt<-ggplot(d.pres[[1]],aes(Maxent,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="Density of Observations",
       title=paste0("Maxent (AUC = ",eval.tab[2,"auc.p"],", PCC = ",eval.tab[2,"pcc.p"],", TSS = ",eval.tab[2,"tss.p"],")"))
disc_p_glm<-ggplot(d.pres[[1]],aes(GLM,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Site",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="",
       title=paste0("GLM (AUC = ",eval.tab[3,"auc.p"],", PCC = ",eval.tab[3,"pcc.p"],", TSS = ",eval.tab[3,"tss.p"],")"))
(disc_p_brt/disc_p_mxt/disc_p_glm)+
  plot_layout(guides = "collect")+
  plot_annotation(title=paste0("N=",nrow(pres_dat[pres_dat$group!=1&pres_dat$y==1,]),", Prevalance = ",round(nrow(pres_dat[pres_dat$group!=1&pres_dat$y==1,])/nrow(pres_dat[pres_dat$group!=1,]),2)))


ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_plots_all_mods_pres.jpeg"),
       width=8,height=10,dpi=300,units = "in")


  # then for the best modeling method, compare across regions
pres_regions<-d.pres[[1]]%>%
  left_join(pres_dat[,c("id","region")],by="id")%>%
  mutate(y=ifelse(obs==0,"Absent","Present"),
         region=case_when(region==1~"Maine",
                          region==3~"Connecticut",
                          region==4~"New York/Long Island",
                          region==5~ "New Jersey"),
         region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")))%>%
  filter(!is.na(region))

disc_p_regions<-ggplot(pres_regions,aes(BRTs,color=as.factor(y),fill=as.factor(y)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Success",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations")+
  facet_wrap(~region,scales = "free")
disc_p_regions+plot_annotation(title="BRT Nest Site Discrimination by Region")

#ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_region_plots_pres_brt.jpeg"),
#       width=8,height=8,dpi=300,units = "in")



## nest survival
disc_s_brt<-ggplot(d.surv[[1]],aes(BRTs,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="",
       title=paste0("BRTs (AUC = ",eval.tab[1,"auc.s"],", PCC = ",eval.tab[1,"pcc.s"],", TSS = ",eval.tab[1,"tss.s"],")"))
disc_s_mxt<-ggplot(d.surv[[1]],aes(Maxent,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="Density of Observations",
       title=paste0("Maxent (AUC = ",eval.tab[2,"auc.s"],", PCC = ",eval.tab[2,"pcc.s"],", TSS = ",eval.tab[2,"tss.s"],")"))
disc_s_glm<-ggplot(d.surv[[1]],aes(GLM,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Success",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="",
       title=paste0("GLM (AUC = ",eval.tab[3,"auc.s"],", PCC = ",eval.tab[3,"pcc.s"],", TSS = ",eval.tab[3,"tss.s"],")"))
(disc_s_brt/disc_s_mxt/disc_s_glm)+
  plot_layout(guides = "collect")+
  plot_annotation(title=paste0("N=",nrow(surv_dat[surv_dat$group!=1&surv_dat$y==1,]),", Prevalance = ",round(nrow(surv_dat[surv_dat$group!=1&surv_dat$y==1,])/nrow(surv_dat[surv_dat$group!=1,]),2)))

ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_plots_all_mods_surv.jpeg"),
       width=8,height=10,dpi=300,units = "in")


  # then for the best model, compare across regions
surv_regions<-d.surv[[5]]%>%
  left_join(pres_dat[,c("id","region")],by="id")%>%
  mutate(y=ifelse(obs==0,"Failed","Fledged"),
         region=case_when(region==1~"Maine",
                          region==3~"Connecticut",
                          region==4~"New York/Long Island",
                          region==5~ "New Jersey"),
         region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")))



disc_s_regions<-ggplot(surv_regions,aes(BRTs,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Success",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations")+
  facet_wrap(~region,scales = "free")
disc_s_regions+plot_annotation(title="BRT Nest Survival Discrimination by Region")

#ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_region_plots_surv_brt.jpeg"),
#       width=8,height=8,dpi=300,units = "in")


#AUC plots
jpeg(paste0(path_out,"Final_outputs/Model_Results/AUC_plots.jpeg"))
par(oma = c(0, 0, 0, 0), mfrow = c(1, 2), cex = 0.7, cex.lab = 1.5)

auc.roc.plot(d.pres[[5]], which.model=1,color = T, legend.cex = 0.4, mark=5,main="")
mtext("Presence", side = 3, line = 2.5, cex = 1.6)
mtext(paste0("(", round(nrow(pres_dat[pres_dat$group!=5&pres_dat$y==1,])/nrow(pres_dat[pres_dat$group!=5,]),2),
             " Prevalence )"), side = 3, line = 0.5, cex = 1.3)
auc.roc.plot(d.surv[[5]], which.model=1,color = T, legend.cex = 0.4, mark=5,main="")
mtext("Survival", side = 3, line = 2.5, cex = 1.6)
mtext(paste0("(", round(nrow(surv_dat[surv_dat$group!=5&surv_dat$y==1,])/nrow(surv_dat[surv_dat$group!=5,]),2),
             " Prevalence )"), side = 3, line = 0.5, cex = 1.3)

dev.off()




## Visual of Calibration evaluation (important if making probability maps)
#--------------------------------------------------------------------------------------------------------------
# groups locations into bins based on their predicted probabilities
# then calculates the ratio of observations with presence vs total observations in that bin
# gives confidence interval for each bin and N of observations in each bin
jpeg(paste0(path_out,"Final_outputs/Model_Results/calibration_plots_",ab_type,".jpeg"))
par(oma = c(0, 5, 0, 0), mar = c(3, 3, 3, 1), mfrow = c(2,3), cex = 0.7, cex.lab = 1.4, mgp = c(2, 0.5, 0))

  for (mod in 1:3) {
    calibration.plot(d.pres[[5]], which.model = mod, color = mod+1, main = model.names[mod], xlab = "", ylab = "")
    if (mod == 1) {
     mtext("Presence", side = 2, line = 6.3, cex = 1.6)
     mtext(paste("(", round(nrow(pres_dat[pres_dat$group!=5&pres_dat$y==1,])/nrow(pres_dat[pres_dat$group!=5,]),2),
                     "Prevalence )"), side = 2, line = 4.3, cex = 1.3)
    }
  }
for (mod in 1:3) {
  calibration.plot(d.surv[[5]], which.model = mod, color = mod+1, main = model.names[mod], xlab = "", ylab = "")
  if (mod == 1) {
    mtext("Survival", side = 2, line = 6.3, cex = 1.6)
    mtext(paste("(", round(nrow(surv_dat[surv_dat$group!=5&surv_dat$y==1,])/nrow(surv_dat[surv_dat$group!=5,]),2),
                "Prevalence )"), side = 2, line = 4.3, cex = 1.3)
  }
}


mtext("Predicted Probability of Occurrence", side = 1, line = -1,
       cex = 1.4, outer = TRUE)
mtext("Observed Occurrence as Proportion of Sites Surveyed",
          side = 2, line = -1, cex = 1.4, outer = TRUE)

dev.off()



## Variable Importance (for the best modeling method)
#-----------------------------------------------------------------------------------

var_p<-left_join(contrib_pres[[1]],contrib_pres[[2]],by="var")%>%
                       left_join(contrib_pres[[3]],by="var")%>%
                       left_join(contrib_pres[[4]], by="var")%>%
                       left_join(contrib_pres[[5]],by="var")%>%
  pivot_longer(cols = -1,names_to = "fold",values_to = "importance")


p5<-ggerrorplot(var_p, x = "var", y = "importance", 
                desc_stat = "mean_sd", color = "black",
                add = "jitter", add.params = list(color = "darkgray"))+
  coord_flip()+
  labs(y="",x="")#,title="Probability of Nest Presence")


var_s<-left_join(contrib_surv[[1]],contrib_surv[[2]],by="var")%>%
  left_join(contrib_surv[[3]],by="var")%>%
  left_join(contrib_surv[[4]], by="var")%>%
  left_join(contrib_surv[[5]],by="var")%>%
  pivot_longer(cols = -1,names_to = "fold",values_to = "importance")

p6<-ggerrorplot(var_s, x = "var", y = "importance", 
                desc_stat = "mean_sd", color = "black",
                add = "jitter", add.params = list(color = "darkgray"))+
  coord_flip()+
  labs(y="Variable Importance",x="")#,title="Probability of Nesting Success")

(p5/p6)+plot_annotation(tag_levels = "A")

ggsave(paste0(path_out,"Final_outputs/Model_Results/var_importance_pres_surv_brt.png"),
       width=8,height=8,dpi=300,units = "in")




# Predict
# Predict to spatial surface, don't build models
predict.surf<-T
build<-F
# Run only the best method
source("05c_boosted_regression_trees.R")

