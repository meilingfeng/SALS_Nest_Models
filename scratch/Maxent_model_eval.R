
me_p <- maxent(dplyr::select(pres_dat_train,all_of(all_terms)), # predictors (to fit categorical variables, indicate them using ",factors='variablename'")
               dplyr::select(pres_dat_train,y)) # occurrences

me_s <- maxent(dplyr::select(surv_dat_train,all_of(all_terms)), 
               dplyr::select(surv_dat_train,y))



# get fitted values


ptest_p <- pres_dat_test[pres_dat_test$y==1,]%>%dplyr::select(all_of(terms))#nest presence
atest_p <- pres_dat_test[pres_dat_test$y==0,]%>%dplyr::select(all_of(terms))#nest absence

ptest_s <- surv_dat_test[surv_dat_test$y==1,]%>%dplyr::select(all_of(terms))#nest success
atest_s <- surv_dat_test[surv_dat_test$y==0,]%>%dplyr::select(all_of(terms))#nest fail


#presence
testp_p <- data.frame(pred=predict(me_p, ptest_p),
                      y=rep("Present",nrow(ptest_p)))
testa_p <- data.frame(pred=predict(me_p, atest_p),
                      y=rep("Absent",nrow(atest_p)))
test_p<-rbind(testp_p,testa_p)

#survival
testp_s <- data.frame(preds=predict(me_s, ptest_s),
                      y=rep("Fledged",nrow(ptest_s)))
testa_s <- data.frame(preds=predict(me_s, atest_s),
                      y=rep("Failed",nrow(atest_s)))
test_s<-rbind(testp_s,testa_s)

#presence
me_acc_p[[i]]<-ggplot(test_p,aes(preds,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Predicted Probability of Nest Site",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations",title=paste0("Maxent (N=",nrow(testp_p),", ",nrow(testa_p)," (Present, Absent)"))

#ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_pres_",ab_type,"_mxnt.png"),
#      width=8,height=10,dpi=300,units = "in")

#survival
me_acc_s[[i]]<-ggplot(test_s,aes(preds,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Predicted Probability of Nest Success",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations",title=paste0("Maxent (N=",nrow(testp_s),", ",nrow(testa_s)," (Present, Absent)"))

#ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_surv_",ab_type,"_mxnt.png"),
#      width=4,height=6,dpi=300,units = "in")





# get model performance metrics (AUC)
auc_p[[i]] <- evaluate(me_p, p=ptest_p, a=atest_p)@auc
auc_s[[i]] <- evaluate(me_s, p=ptest_s, a=atest_s)@auc

#threshold(e_p)

#plot(e_p, 'ROC')
}

auc_p<-unlist(auc_p)
auc_p<-paste0("AUC = ",round(mean(auc_p),3),", 95% CI ",round(mean(auc_p)+1.96*(sd(auc_p)/length(auc_p)),3),"-",round(mean(auc_p)-1.96*(sd(auc_p)/length(auc_p)),3))
auc_p

auc_s<-unlist(auc_s)
auc_s<-paste0("AUC = ",round(mean(auc_s),3),", 95% CI ",round(mean(auc_s)+1.96*(sd(auc_s)/length(auc_p)),3),"-",round(mean(auc_s)-1.96*(sd(auc_s)/length(auc_s)),3))
auc_s




me_acc_p[[1]]
me_acc_p[[2]]
me_acc_p[[3]]
me_acc_p[[4]]
me_acc_p[[5]]

me_acc_s[[1]]
me_acc_s[[2]]
me_acc_s[[3]]
me_acc_s[[4]]
me_acc_s[[5]]



me_p <- maxent(dplyr::select(pres_dat,all_of(terms)), dplyr::select(pres_dat,y)) 

pres_dat2<-pres_dat%>%mutate(
  preds=predict(me_p, pres_dat),
  y=ifelse(y==0,"Absent","Present"),
  region=case_when(region==1~"Maine",
                   region==3~"Connecticut",
                   region==4~"New York/Long Island",
                   region==5~ "New Jersey"),
  region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")))

me_s <- maxent(dplyr::select(surv_dat,all_of(terms)), dplyr::select(surv_dat,y)) 

surv_dat2<-surv_dat%>%mutate(
  preds=predict(me_s, surv_dat),
  y=ifelse(y==0,"Failed","Fledged"),
  region=case_when(region==1~"Maine",
                   region==3~"Connecticut",
                   region==4~"New York/Long Island",
                   region==5~ "New Jersey"),
  region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")))

#presence
p1<-ggplot(pres_dat2,aes(preds,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations")
p2<-ggplot(pres_dat2,aes(preds,color=as.factor(y),fill=as.factor(y)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Predicted Probability of Nest Success",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations")+
  facet_wrap(~region,scales = "free")
(p1/p2)+plot_layout(guides = "collect")+plot_annotation(title=auc_p)
ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_pres_",ab_type,"_mxnt.png"),
       width=8,height=10,dpi=300,units = "in")


#survival
p3<-ggplot(surv_dat2,aes(preds,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations")
p4<-ggplot(surv_dat2,aes(preds,color=y,fill=y))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Predicted Probability of Nest Success",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations")+
  facet_wrap(~region,scales = "free")
(p3/p4)+plot_layout(guides = "collect")+plot_annotation(title=auc_s)

ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_surv_",ab_type,"_mxnt.png"),
       width=8,height=10,dpi=300,units = "in")

var_p<-as.data.frame(t(data.frame(a=var_imp_p[[1]],b=var_imp_p[[2]], c=var_imp_p[[3]], d=var_imp_p[[4]], e= var_imp_p[[5]])))%>%
  pivot_longer(cols = everything(),names_to = "variable",values_to = "importance")


p5<-ggerrorplot(var_p, x = "variable", y = "importance", 
                desc_stat = "mean_sd", color = "black",
                add = "jitter", add.params = list(color = "darkgray"))+
  coord_flip()+
  labs(y="",x="")


var_s<-as.data.frame(t(data.frame(a=var_imp_s[[1]],b=var_imp_s[[2]], c=var_imp_s[[3]], d=var_imp_s[[4]], e= var_imp_s[[5]])))%>%
  pivot_longer(cols = everything(),names_to = "variable",values_to = "importance")


p6<-ggerrorplot(var_s, x = "variable", y = "importance", 
                desc_stat = "mean_sd", color = "black",
                add = "jitter", add.params = list(color = "darkgray"))+
  coord_flip()+
  labs(y="Variable Contribution (%)",x="")

(p5/p6)+plot_annotation(tag_levels = "A")
ggsave(paste0(path_out,"Final_outputs/Model_Results/var_importance_pres_surv_",ab_type,"_mxnt.png"),
       width=8,height=8,dpi=300,units = "in")
