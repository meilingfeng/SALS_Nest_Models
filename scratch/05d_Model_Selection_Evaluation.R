library(tidyverse)
library(car)#vif
library(mgcv)
library(dismo)
library(sf)
library(patchwork)
library(ggforce)
library(officer)
library(flextable)
library(data.table)

### Set up
# -------------------------------------------
source("05b_Regression_Model_Selection.R")  
  
## Model Evaluation
#---------------------------------------------------

#Spatial sorting bias (SSB) - background points that cover a larger spatial extent have higher AUC
#SSB = 1 is no bias, close to zero is a lot of bias


pres<-filter(pres_dat,y==1)%>%
  dplyr::select(longitude,latitude,group)
backgr<-filter(pres_dat,y==0)%>%
  dplyr::select(longitude,latitude,group)

p_train <- pres[pres$group!=1, ]
p_test <- pres[pres$group==1, ]

b_train <- backgr[backgr$group!=1, ]
b_test <- backgr[backgr$group==1, ]

sb <- ssb(p_test, b_test, p_train)
sb[,1] / sb[,2]

#0.4 not great but not terrible

i <- pwdSample(p_test, b_test, p_train, n=1, tr=0.1)
p_test_pwd <- p_test[!is.na(i[,1]), ]
b_test_pwd <- b_test[na.omit(as.vector(i)), ]
sb_pwd <- ssb(p_test_pwd, b_test_pwd, p_train)
sb[,1] / sb[,2]




# select the top ranked model based on AIC (any models within 2 delta AIC)
mod_pres<-mod_tab[mod_tab$response=="Presence",]
form_pres<-mod_pres[mod_pres$AIC==min(mod_pres$AIC),]$fun
final_mod_pres<-glm(form_pres, data = pres_dat, family = binomial(link = "logit"))

mod_surv<-mod_tab[mod_tab$response=="Success",]
form_surv<-mod_surv[mod_surv$AIC==min(mod_surv$AIC),]$fun
final_mod_surv<-glm(form_surv, data = surv_dat, family = binomial(link = "logit"))

e_p <- list()
for (i in 1:k) {
  train <- pres_dat[pres_dat$group != i,]
  test <- pres_dat[pres_dat$group == i,]
  mod <- glm(form_pres, data = train, family = binomial(link = "logit"))
  e_p[[i]] <- evaluate(p=test[test$y==1,], a=test[test$y==0,], mod)
}
auc_p <- sapply(e_p, function(x){x@auc})

auc_p<-paste0("AUC = ",round(mean(auc_p),3),", 95% CI ",round(mean(auc_p)+1.96*(sd(auc_p)/length(auc_p)),3),"-",round(mean(auc_p)-1.96*(sd(auc_p)/length(auc_p)),3))

e_s <- list()
for (i in 1:k) {
  train <- surv_dat[surv_dat$group != i,]
  test <- surv_dat[surv_dat$group == i,]
  mod <- glm(form_surv, data = train, family = binomial(link = "logit"))
  e_s[[i]] <- evaluate(p=test[test$y==1,], a=test[test$y==0,], mod)
}
auc_s <- sapply(e_s, function(x){x@auc})

auc_s<-paste0("AUC = ",round(mean(auc_s),3),", 95% CI ",round(mean(auc_s)+1.96*(sd(auc_s)/length(auc_s)),3),"-",round(mean(auc_s)-1.96*(sd(auc_s)/length(auc_s)),3))




### Model selection table
# variable order for plotting and tables
hypoth<-c('Additional Habitat Characteristics', 'High Marsh Quality','Unclassified Habitat Characteristics', 'Geographic Variation')
mod_names<-c("High Marsh","Vegetation Health (NDVI)", "Marsh Resilience (UVVR)", "Habitat Dissimilarity (Entropy)", "High Marsh + NDVI" ,"High Marsh + UVVR","High Marsh + Entropy","Full Additional Habitat",
             "High Marsh * Unclassified Habitat",
             "Unclassified Habitat","High Marsh + Unclassified Habitat", "Additional Habitat + Unclassified Habitat","Latitude", "High Marsh + Latitude","Additional Habitat + Latitude","Unclassified Habitat + Latitude", "Additional and Unclassified Habitat + Latitude")


tab<-mod_tab%>%
  mutate(class=factor(class,levels=hypoth),
         name=factor(name,levels=mod_names))
tab<-tab%>%
  arrange(response,class,name)%>%
  dplyr::select(-max_VIF,-VIF_var)

tab_p<-filter(tab,response=="Presence")%>%dplyr::select(-response)
tab_s<-filter(tab,response=="Success")%>%dplyr::select(-response)

tab_p_flx<-flextable(tab_p)%>%
  ##specify grouped column header widths and names
  ##how many columns are in each grouped header
  #add_header_row(colwidths = c(1,1,3,3,3,3), 
                 ##Names of each grouped header
                 #values = c("","","BC Priorities",
                  #          "BCH Priorities",
                  #          "BC Priorities with IPLCs",
                  #          "BCH Priorities with IPLCs"))%>%
  #specify sub-column header names (rename variables)
  set_header_labels(class="Hypothesis",
                    name="Model Name", 
                    fun="Model Function",
                    AIC="AIC",
                    dAIC="Hypothesis ΔAIC",
                    overall_dAIC="Overall ΔAIC") %>%
  
  ##Set all missing values to "---"
  #colformat_num(j=1:length(tab), na_str = "---")%>%
  
  #column header font format
  fontsize(size=11, part='header') %>%
  #table values font format
  fontsize(size=10, part='body') %>%
  #cell alignment throughout the whole table
  flextable::align(align="center",part = "all")%>%
  #column widths
  flextable::width(j=c(1,2,3),3,unit = "in")%>%
  flextable::width(j=c(4,5,6),1,unit = "in")%>%
  #Merge the rows by hypothesis 
  merge_v(j=1)%>%
  #add horizontal border lines for each row starting
  border_inner_h(border=fp_border(color='#CAD4D1', style='solid', width=0.5)) %>%
  #add thick horizontal border lines dividing each region section
  #border(j=1:length(tab_p),i=c(8,11,12), 
  #       border.bottom =fp_border(color='#666666', style='solid', width=1.5))%>%
  #Make the background for the top model darker to stand out
  bg(i=which(tab_p$fun==form_pres),bg="#FAE29C")%>%
  #Make the background for the top additional habitat model darker to stand out
  #border(j=1:length(tab_p),i=which(tab_p$class=="Additional Habitat Characteristics" & tab_p$dAIC==0),
  #       border.top=fp_border(color="#FAE29C",style="solid",width=2),
  #       border.bottom=fp_border(color="#FAE29C",style="solid",width=2))%>%
  #Make all other regional sections and the header section white background
  bg(j=1:length(tab_p),i=c(1:(which(tab_p$fun==form_pres)-1),(which(tab_p$fun==form_pres)+1):nrow(tab_p)),bg="white")%>%
  bg(bg="white",part = "header")



if(!exists(paste0(path_out,"Final_outputs/Model_Results/pres_model_selection_table_",ab_type,".png"))){
  save_as_image(tab_p_flx,paste0(path_out,"Final_outputs/Model_Results/pres_model_selection_table_",ab_type,".png"))
}


tab_s_flx<-flextable(tab_s)%>%
  ##specify grouped column header widths and names
  ##how many columns are in each grouped header
  #add_header_row(colwidths = c(1,1,3,3,3,3), 
  ##Names of each grouped header
  #values = c("","","BC Priorities",
  #          "BCH Priorities",
  #          "BC Priorities with IPLCs",
  #          "BCH Priorities with IPLCs"))%>%
  #specify sub-column header names (rename variables)
  set_header_labels(class="Hypothesis",
                    name="Model Name", 
                    fun="Model Function",
                    AIC="AIC",
                    dAIC="Hypothesis ΔAIC",
                    overall_dAIC="Overall ΔAIC") %>%
  
  ##Set all missing values to "---"
  #colformat_num(j=1:length(tab), na_str = "---")%>%
  
  #column header font format
  fontsize(size=11, part='header') %>%
  #table values font format
  fontsize(size=10, part='body') %>%
  #cell alignment throughout the whole table
  flextable::align(align="center",part = "all")%>%
  #column widths
  flextable::width(j=c(1,2,3),3,unit = "in")%>%
  flextable::width(j=c(4,5,6),1,unit = "in")%>%
  #Merge the rows by hypothesis 
  merge_v(j=1)%>%
  #add horizontal border lines for each row starting
  border_inner_h(border=fp_border(color='#CAD4D1', style='solid', width=0.5)) %>%
  #add thick horizontal border lines dividing each region section
  #border(j=1:length(tab_s),i=c(8,11,12), 
  #       border.bottom =fp_border(color='#666666', style='solid', width=1.5))%>%
  #Make the background for the top model darker to stand out
  bg(i=which(tab_s$fun==form_surv),bg="#FAE29C")%>%
  #Make the background for the top additional habitat model darker to stand out
  #border(j=1:length(tab_s),i=which(tab_s$class=="Additional Habitat Characteristics" & tab_s$dAIC==0),
  #       border.top=fp_border(color="#FAE29C",style="solid",width=2),
  #       border.bottom=fp_border(color="#FAE29C",style="solid",width=2))%>%
#Make all other regional sections and the header section white background
bg(j=1:length(tab_s),i=c(1:(which(tab_s$fun==form_surv)-1),(which(tab_s$fun==form_surv)+1):nrow(tab_s)),bg="white")%>%
bg(bg="white",part = "header")



if(!exists(paste0(path_out,"Final_outputs/Model_Results/surv_model_selection_table_",ab_type,".png"))){
  save_as_image(tab_s_flx,paste0(path_out,"Final_outputs/Model_Results/surv_model_selection_table_",ab_type,".png"))
}


## plot model selection

p1<-mod_tab%>%
  filter(response=="Success")%>%
  mutate(overall_dAIC_plot=case_when(overall_dAIC==0~0.2,overall_dAIC!=0 ~overall_dAIC),
         class=factor(class, levels=hypoth),
         name=factor(name, levels = rev(mod_names)))%>%
  dplyr::select(name,overall_dAIC,overall_dAIC_plot, class)%>%
  ggplot(aes(y=name,x=overall_dAIC_plot, fill=overall_dAIC_plot)) +
  theme_classic(base_size = 12) +
  geom_col(width = 0.9, position = position_dodge(0.3),)+
  geom_text(aes(label = overall_dAIC), hjust = -0.5)+
  scale_fill_viridis_c(direction=-1)+
  labs(y="Models",x="ΔAIC",fill="ΔAIC")+
  facet_col(~class, scales = "free_y",space="free")



p2<-mod_tab%>%
  filter(response=="Presence")%>%
  mutate(overall_dAIC_plot=case_when(overall_dAIC==0~10,overall_dAIC!=0~overall_dAIC),
         class=factor(class, levels=hypoth),
         name=factor(name, levels = rev(mod_names)))%>%
  dplyr::select(name,overall_dAIC,overall_dAIC_plot, class)%>%
  ggplot(aes(y=name,x=overall_dAIC_plot, fill=overall_dAIC_plot)) +
  theme_classic(base_size = 12) +
  geom_col(width = 0.9, position = position_dodge(0.3),)+
  geom_text(aes(label = overall_dAIC), hjust = -0.5)+
  scale_fill_viridis_c(direction=-1)+
  labs(y="Models",x="ΔAIC",fill="ΔAIC")+
  facet_col(~class, scales = "free_y",space="free")


p1
ggsave(filename=paste0(path_out,"Final_outputs/Model_Results/Candidate_model_AIC_success_",ab_type,".png"), width = 10, height = 6, dpi = "retina")


p2
ggsave(filename=paste0(path_out,"Final_outputs/Model_Results/Candidate_model_AIC_placement_",ab_type,".png"), width = 16, height = 6, dpi = "retina")





### Plot predicted (fitted values) against original 0,1s to see how they correlate

#format data for plotting
pres_dat2<-pres_dat%>%mutate(preds=round(final_mod_pres$fitted.values,2),
                                 region=case_when(region==1~"Maine",
                                                  region==3~"Connecticut",
                                                  region==4~"New York/Long Island",
                                                  region==5~ "New Jersey"),
                                 region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")),
                             y=ifelse(y==0,"Absent","Present"))%>%
  filter(!is.na(region))

surv_dat2<-surv_dat%>%mutate(preds=round(final_mod_surv$fitted.values,2),
                                 region=case_when(region==1~"Maine",
                                                  region==3~"Connecticut",
                                                  region==4~"New York/Long Island",
                                                  region==5~ "New Jersey"),
                                 region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey")),
                             y=ifelse(y==0,"Failed","Fledged"))%>%
  filter(!is.na(region))


# scatter plots
#ggplot(surv_train2, aes(x=preds,y=y,group=as.factor(region)))+
#  geom_jitter(height=0.05,aes(color=as.factor(region)))+
#  geom_smooth(method="glm",method.args=list(family="binomial"),aes(colour=as.factor(region)))

#ac_s<-ggplot(surv_train2%>%mutate(region=case_when(region==1~"Maine",
#                                            region==3~"Connecticut",
#                                            region==4~"New York/Long Island",
#                                            region==5~ "New Jersey"),
 #                          region= factor(region,levels=c("Maine","Connecticut","New York/Long Island", "New Jersey"))), 
#       aes(x=preds,y=y))+
#  geom_jitter(height=0.05,alpha=0.5)+
#  geom_smooth(method="glm",method.args=list(family="binomial"))+
#  scale_y_continuous(labels=c("Failed","Fledged"), breaks=c(0,1))+
#  labs(x="Fitted Values",y="Observed Nest Success")+
#  theme_bw(base_size=12)+
#  facet_wrap(~region,scales="free_x",nrow=1)

#ac_p<-ggplot(pres_train2, 
#       aes(x=preds,y=y))+
#  geom_jitter(height=0.05,alpha=0.5)+
#  geom_smooth(method="glm",method.args=list(family="binomial"))+
#  scale_y_continuous(labels=c("Absent","Present"), breaks=c(0,1))+
#  labs(x="",y="Observed Nest Presence")+
#  theme_bw(base_size=12)+
#  facet_wrap(~region,scales="free_x",nrow=1)

#(ac_p/ac_s)
#ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_surv_pres_",ab_type,".png"),
#       width=18,height=7,dpi=300,units = "in")



#presence
p4<-ggplot(pres_dat2,aes(preds,color=as.factor(y),fill=as.factor(y)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations",title="All Regions")

p5<-ggplot(pres_dat2,aes(preds,color=as.factor(y),fill=as.factor(y)))+
  geom_density(alpha=0.5,adjust=1)+
  scale_fill_viridis_d(end = 0.8)+
  scale_color_viridis_d(end=0.75)+
  theme_bw(base_size=12)+
  scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
  labs(x="Predicted Probability of Nest Site",fill="Observed Nest Site",color="Observed Nest Site",y="Density of Observations")+
  facet_wrap(~region,scales = "free")

(p4/p5)+plot_layout(guides = "collect")+plot_annotation(title=auc_p)

ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_pres_",ab_type,".png"),
       width=8,height=10,dpi=300,units = "in")

  #survival
p7<-  ggplot(surv_dat2,aes(preds,color=as.factor(y),fill=as.factor(y)))+
     geom_density(alpha=0.5,adjust=1)+
     scale_fill_viridis_d(end = 0.8)+
     scale_color_viridis_d(end=0.75)+
     theme_bw(base_size=12)+
     scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
     labs(x="",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations",title="All Regions")
  
p8<-ggplot(surv_dat2,aes(preds,color=as.factor(y),fill=as.factor(y)))+
     geom_density(alpha=0.5,adjust=1)+
     scale_fill_viridis_d(end = 0.8)+
     scale_color_viridis_d(end=0.75)+
     theme_bw(base_size=12)+
     scale_x_continuous(limits = c(-0.2,1.2),breaks=c(0,.5,1))+
     labs(x="Predicted Probability of Nest Success)",fill="Observed Nest Success",color="Observed Nest Success",y="Density of Observations")+
     facet_wrap(~region,scales = "free")
(p7/p8)+plot_layout(guides = "collect")+plot_annotation(title=auc_s)



ggsave(paste0(path_out,"Final_outputs/Model_Results/observed_fitted_density_surv_",ab_type,".png"),
       width=8,height=10,dpi=300,units = "in")

