library(tidyverse)
library(dismo)
library(gbm)

## Set up 
#------------------------------------------------------------------
#Load data 
source("05a_Set_Model_Parameters.R")

# Don't predict to spatial surface yet
predict.surf<-F
# Run only the test data predictions using 3 modeling methods
source("05c_boosted_regression_trees.R")
source("05d_Maxent.R")
source("05b_GLM_Selection_Prediction.R")

model.names<-c("BRTs","Maxent","GLM")


#combine each model's predictions into one df per fold
d.pres<-list() #list of dfs for each fold for predicting nest presence
d.surv<-list() #list of dfs for each fold for predicting nest survival
for(i in 1:k){
d.pres[[i]]<-cbind(d.pres.brt[[i]],d.pres.mxt[[i]][,3],d.pres.glm[[i]][,3])%>%
  rename("Maxent"="d.pres.mxt[[i]][, 3]","GLM"="d.pres.glm[[i]][, 3]","BRTs"="pred")

d.surv[[i]]<-cbind(d.surv.brt[[i]],d.surv.mxt[[i]][,3],d.surv.glm[[i]][,3])%>%
  rename("Maxent"="d.surv.mxt[[i]][, 3]","GLM"="d.surv.glm[[i]][, 3]","BRTs"="pred")
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

accu.pres<-list()
accu.surv<-list()

# calculate metrics for each fold
for(i in 1:k){
# get optimal threshold for Maximizing Kappa
thr.p<-optimal.thresholds(d.pres[[i]],opt.methods = "MaxKappa")
thr.pres[[i]]<-c(thr.p[1,]$BRTs,thr.p[1,]$Maxent,thr.p[1,]$GLM)

thr.s<-optimal.thresholds(d.surv[[i]],opt.methods = "MaxKappa")
thr.surv[[i]]<-c(thr.s[1,]$BRTs,thr.s[1,]$Maxent,thr.s[1,]$GLM)


# accuracy metric tables
accu.p<-presence.absence.accuracy(d.pres[[i]],threshold = thr.pres[[i]])
accu.p[, -c(1, 2)] <- signif(accu.p[, -c(1, 2)], digits = 3)
accu.pres[[i]]<-accu.p

accu.s<-presence.absence.accuracy(d.surv[[i]],threshold=thr.surv[[i]])
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
# BRT has the best performance across all metrics. 
# In general worse predictions for survival. 
# Both survival and presence have better discrimination than calibration.
# high AUC ... all these metrics are influenced by prevalence. Increasing the area of background can increase the scores. Argument for veg points? try extrapolating predictor rasters to points outside their range.


## Visual of Discrimination evaluation (important if making presence absence maps) - histogram of presences to absences 
#---------------------------------------------------------------------------------------------------------------------------------
## nest presence
  # first compare across models
disc_p_brt<-ggplot(test_p,aes(BRTs,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="",
       title=paste0("BRTs (AUC = ",eval.tab[1,"auc.p"],", MaxKappa = ",eval.tab[1,"kappa.p"],", TSS = ",eval.tab[1,"tss.p"],")"))
disc_p_mxt<-ggplot(test_p,aes(Maxent,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="Density of Observations",
       title=paste0("Maxent (AUC = ",eval.tab[2,"auc.p"],", MaxKappa = ",eval.tab[2,"kappa.p"],", TSS = ",eval.tab[2,"tss.p"],")"))
disc_p_glm<-ggplot(test_p,aes(GLM,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Site",
       fill="Observed Nest Site",
       color="Observed Nest Site",
       y="",
       title=paste0("GLM (AUC = ",eval.tab[3,"auc.p"],", MaxKappa = ",eval.tab[3,"kappa.p"],", TSS = ",eval.tab[3,"tss.p"],")"))
(disc_p_brt/disc_p_mxt/disc_p_glm)+
  plot_layout(guides = "collect")+
  plot_annotation(title=paste0("N=",nrow(pres_dat[pres_dat$group!=5&pres_dat$y==1,]),", Prevalance = ",round(nrow(pres_dat[pres_dat$group!=5&pres_dat$y==1,])/nrow(pres_dat[pres_dat$group!=5,]),2)))


ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_plots_all_mods_pres_",ab_type,".jpeg"),
       width=8,height=10,dpi=300,units = "in")


  # then for the best model, compare across regions
pres_regions<-d.pres[[5]]%>%
  left_join(pres_dat[,c("id","region")],by="id")%>%
  mutate(y=ifelse(obs==0,"Absent","Present"),
         region=case_when(region==1~"Maine",
                          region==3~"Connecticut",
                          region==4~"New York/Long Island",
                          region==5~ "New Jersey"),
         region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")))

disc_p_regions<-ggplot(pres_regions,aes(BRTs,color=as.factor(y),fill=as.factor(y)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Success",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations")+
  facet_wrap(~region,scales = "free")
disc_p_regions+plot_annotation(title="BRT Nest Site Discrimination by Region")

ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_region_plots_pres_",ab_type,"_brt.jpeg"),
       width=8,height=8,dpi=300,units = "in")



## nest survival
disc_s_brt<-ggplot(test_s,aes(BRTs,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="",
       title=paste0("BRTs (AUC = ",eval.tab[1,"auc.s"],", MaxKappa = ",eval.tab[1,"kappa.s"],", TSS = ",eval.tab[1,"tss.s"],")"))
disc_s_mxt<-ggplot(test_s,aes(Maxent,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="Density of Observations",
       title=paste0("Maxent (AUC = ",eval.tab[2,"auc.s"],", MaxKappa = ",eval.tab[2,"kappa.s"],", TSS = ",eval.tab[2,"tss.s"],")"))
disc_s_glm<-ggplot(test_s,aes(GLM,color=as.factor(obs),fill=as.factor(obs)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Probability of Nest Success",
       fill="Observed Nest Success",
       color="Observed Nest Success",
       y="",
       title=paste0("GLM (AUC = ",eval.tab[3,"auc.s"],", MaxKappa = ",eval.tab[3,"kappa.s"],", TSS = ",eval.tab[3,"tss.s"],")"))
(disc_s_brt/disc_s_mxt/disc_s_glm)+
  plot_layout(guides = "collect")+
  plot_annotation(title=paste0("N=",nrow(surv_dat[surv_dat$group!=5&surv_dat$y==1,]),", Prevalance = ",round(nrow(surv_dat[surv_dat$group!=5&surv_dat$y==1,])/nrow(surv_dat[surv_dat$group!=5,]),2)))

ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_plots_all_mods_surv_",ab_type,".jpeg"),
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

ggsave(paste0(path_out,"Final_outputs/Model_Results/discrim_region_plots_surv_",ab_type,"_brt.jpeg"),
       width=8,height=8,dpi=300,units = "in")


#AUC plots
par(oma = c(0, 0, 0, 0), mfrow = c(1, 2), cex = 0.7, cex.lab = 1.5)

auc.roc.plot(d.pres[[5]], which.model=1,color = T, legend.cex = 0.4, mark=5,main="")
mtext("Presence", side = 3, line = 2.5, cex = 1.6)
mtext(paste0("(", round(nrow(pres_dat[pres_dat$group!=5&pres_dat$y==1,])/nrow(pres_dat[pres_dat$group!=5,]),2),
             " Prevalence )"), side = 3, line = 0.5, cex = 1.3)
auc.roc.plot(d.surv[[5]], which.model=1,color = T, legend.cex = 0.4, mark=5,main="")
mtext("Survival", side = 3, line = 2.5, cex = 1.6)
mtext(paste0("(", round(nrow(surv_dat[surv_dat$group!=5&surv_dat$y==1,])/nrow(surv_dat[surv_dat$group!=5,]),2),
             " Prevalence )"), side = 3, line = 0.5, cex = 1.3)






## Visual of Calibration evaluation (important if making probability maps)
#--------------------------------------------------------------------------------------------------------------
# groups locations into bins based on their predicted probabilities
# then calculates the ratio of observations with presence vs total observations in that bin
# gives confidence interval for each bin and N of observations in each bin
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
  labs(y="",x="",title="Nest Presence BRTs")


var_s<-left_join(contrib_surv[[1]],contrib_surv[[2]],by="var")%>%
  left_join(contrib_surv[[3]],by="var")%>%
  left_join(contrib_surv[[4]], by="var")%>%
  left_join(contrib_surv[[5]],by="var")%>%
  pivot_longer(cols = -1,names_to = "fold",values_to = "importance")

p6<-ggerrorplot(var_s, x = "var", y = "importance", 
                desc_stat = "mean_sd", color = "black",
                add = "jitter", add.params = list(color = "darkgray"))+
  coord_flip()+
  labs(y="Variable Importance",x="",title="Nest Survival BRTs")

(p5/p6)+plot_annotation(tag_levels = "A")

ggsave(paste0(path_out,"Final_outputs/Model_Results/var_importance_pres_surv_",ab_type,"_brt.png"),
       width=8,height=8,dpi=300,units = "in")
