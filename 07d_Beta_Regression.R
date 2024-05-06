library(glmmTMB)
library(StepBeta)
library(ggstats)
library(usdm)
library(glmmTMB)
library(buildmer)
library(broom.mixed)
library(ggplot2)
library(dplyr)
## 1. Data Preparation for following model building:
#-------------------------------------------------------
speciesnames<-c("SALS","SESP","CLRA","WILL","NESP","HYBR")
for (s in 1:length(speciesnames)){
  rapid_dat10.model=read.csv(paste0(path_out,"Final_outputs/",speciesnames[s],"_Nest_Predictions_at_Rapid_Survey_Points.csv"))
  
  #Convert 0s and 1s of nest presence and nestling survival in order to build beta regression models later.
  #the equation is from: https://psycnet.apa.org/doiLanding?doi=10.1037%2F1082-989X.11.1.54
  #N=length(rapid_dat10$presence)
  #rapid_dat10.model$presence=(rapid_dat10.model$presence*(N-1)+0.5)/N
  #rapid_dat10.model$survival=(rapid_dat10.model$survival*(N-1)+0.5)/N
  #centering Latitude
  rapid_dat10.model[,c(5,19:32)]=scale(rapid_dat10.model[,c(5,19:32)])
  #Check collinearity problems of predictor variables
  rapid_dat10.vif<-as.data.frame(rapid_dat10.model[,c(19:32)])
  vifcor(rapid_dat10.vif,th=0.9)
  vifstep(rapid_dat10.vif,th=10)
  #No collineraity problems

  
## 2. Perform stepped beta regression modeling 
#--------------------------------------------------------
#  2.1 for Nesting Survival
  #summary the global model using all predictor variables (percentage)

  summary(fullmod_survival_pct<-glmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
        patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
        high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
        upland_pct+trees_pct+water_pct+Lat+I(Lat^2),data=rapid_dat10.model,family=ordbeta(link = "logit")))
  #use stepmod function to simplify the global model
  stepmod_survival_pct<-buildglmmTMB(survival~alt_tall_pct+distichlis_pct+gerardii_pct+
                                   patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                                   high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                                   upland_pct+trees_pct+water_pct+Lat+I(Lat^2),data=rapid_dat10.model,family=ordbeta(link = "logit"))
  #check the output (saltmarsh_border_pct is not significant)
  summary(stepmod_survival_pct)
  #record the simplified model
#  simplified_surv_pct_Lat2<-betareg(survival ~ alt_tall_pct+ distichlis_pct+gerardii_pct +patens_pct + alt_short_pct+ phrag_pct +low_marsh_pct+ high_marsh_pct + brackish_border_pct  + upland_pct + water_pct  + Lat+I(Lat^2), data = rapid_dat10.model )
  #record the simplified model with Lat^2 removed for further comparison
  simplified_surv_pct<-stepmod_survival_pct@model
  
#  2.2 for Nesting Presence
  
  ##summary the global model using all predictor variables (percentage)
  summary(fullmod_pres_pct<-glmmTMB(
    presence~alt_tall_pct+distichlis_pct+gerardii_pct+
      patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
      high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
      trees_pct+water_pct+upland_pct+Lat+I(Lat^2),data=rapid_dat10.model,family=ordbeta(link = "logit")))
  #use stepmod function to simplify the global model
  stepmod_pres_pct<-buildglmmTMB(presence~alt_tall_pct+distichlis_pct+gerardii_pct+
                 patens_pct+alt_short_pct+phrag_pct+low_marsh_pct+
                 high_marsh_pct+brackish_border_pct+saltmarsh_border_pct+
                 upland_pct+trees_pct+water_pct+Lat+I(Lat^2),data=rapid_dat10.model,family=ordbeta(link = "logit"))
  #check the output (saltmarsh_border_pct is not significant)
  summary(stepmod_pres_pct)
  #record the simplified model
#  simplified_pres_pct_Lat2<-glmmTMB(formula =  presence ~  alt_tall_pct + gerardii_pct+patens_pct + alt_short_pct +phrag_pct  + high_marsh_pct +trees_pct+Lat +I(Lat^2),data = rapid_dat10.model)
  #record the simplified model with Lat^2 removed for further comparison
  simplified_pres_pct<-stepmod_pres_pct@model
  
# 2.3 coefficient plots
  ggcoef_model(fullmod_survival_pct)
  ggsave(paste0(path_out,"Plots/",speciesnames[s],"_full_surv_coef.png"),width=12,height=6,units = "in")
  ggcoef_model(simplified_surv_pct)
  ggsave(paste0(path_out,"Plots/",speciesnames[s],"_simp_surv_coef.png"),width=12,height=6,units = "in")
  
 # ggcoef_model(simplified_surv_pct_Lat2)
  ggcoef_model(fullmod_pres_pct)
  ggsave(paste0(path_out,"Plots/",speciesnames[s],"_full_pres_coef.png"),width=12,height=6,units = "in")
  
  ggcoef_model(simplified_pres_pct)
  ggsave(paste0(path_out,"Plots/",speciesnames[s],"_simp_pres_coef.png"),width=12,height=6,units = "in")
  
 # ggcoef_model(simplified_pres_pct_Lat2)  
  
  
## 3. Compare different beta regression models
#-------------------------------------------------------
#  3.1 for Nesting Survival
  ## Mod List
  ## 1. Select model 
  #------------------------------------------
  ## Start a list of potential models
  mod_list_surv1<-list()
  # How much variation does each variable explain relative to each other?
  #alt_tall
  mod_list_surv1[[1]]<-glmmTMB(survival~alt_tall_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[2]]<-glmmTMB(survival~alt_tall_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #distichlis
  mod_list_surv1[[3]]<-glmmTMB(survival~distichlis_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[4]]<-glmmTMB(survival~distichlis_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #gerardii
  mod_list_surv1[[5]]<-glmmTMB(survival~gerardii_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[6]]<-glmmTMB(survival~gerardii_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #patens
  mod_list_surv1[[7]]<-glmmTMB(survival~patens_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[8]]<-glmmTMB(survival~patens_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #alt_short
  mod_list_surv1[[9]]<-glmmTMB(survival~alt_short_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[10]]<-glmmTMB(survival~alt_short_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #phrag
  mod_list_surv1[[11]]<-glmmTMB(survival~phrag_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[12]]<-glmmTMB(survival~phrag_pct+Lat+I(Lat^2),
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #HIGHMARSH
  mod_list_surv1[[13]]<-glmmTMB(survival~high_marsh_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[14]]<-glmmTMB(survival~high_marsh_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Barckish Border
  mod_list_surv1[[15]]<-glmmTMB(survival~brackish_border_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[16]]<-glmmTMB(survival~brackish_border_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Non-saltmarsh
  
  mod_list_surv1[[17]]<-glmmTMB(survival~upland_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[18]]<-glmmTMB(survival~upland_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[19]]<-glmmTMB(survival~water_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[20]]<-glmmTMB(survival~water_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[21]]<-glmmTMB(survival~low_marsh_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[22]]<-glmmTMB(survival~low_marsh_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[23]]<-glmmTMB(survival~saltmarsh_border_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[24]]<-glmmTMB(survival~saltmarsh_border_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[25]]<-glmmTMB(survival~trees_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[26]]<-glmmTMB(survival~trees_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_surv1[[27]]<- simplified_surv_pct
  mod_list_surv1[[28]]<- fullmod_survival_pct
 # mod_list_surv1[[27]]<- simplified_surv_pct_Lat2

  # store model selection factors in a table
  mod_tab_habstr<-data.frame(Response=rep("survival",length(mod_list_surv1)),
                             Model_Comparison=rep("Veg Species Occurrence",length(mod_list_surv1)),
                             Model_Name=rep(NA,length(mod_list_surv1)),
                             Function=rep(NA,length(mod_list_surv1)),
                             AIC=rep(NA,length(mod_list_surv1)),
                             dAIC=rep(NA,length(mod_list_surv1)))
  
  mod_tab_habstr[1,"Model_Name"]<-"alt tall pres/abs"
  mod_tab_habstr[1,"Function"]<-deparse1(mod_list_surv1[[1]]$formula)
  mod_tab_habstr[1,"AIC"]<-AIC(mod_list_surv1[[1]])
  
  mod_tab_habstr[2,"Model_Name"]<-"alt tall pct"
  mod_tab_habstr[2,"Function"]<-deparse1(mod_list_surv1[[2]]$formula)
  mod_tab_habstr[2,"AIC"]<-AIC(mod_list_surv1[[2]])
  
  mod_tab_habstr[3,"Model_Name"]<-"distichlis pres/abs"
  mod_tab_habstr[3,"Function"]<-deparse1(mod_list_surv1[[3]]$formula)
  mod_tab_habstr[3,"AIC"]<-AIC(mod_list_surv1[[3]])
  
  mod_tab_habstr[4,"Model_Name"]<-"distichlis pct"
  mod_tab_habstr[4,"Function"]<-deparse1(mod_list_surv1[[4]]$formula)
  mod_tab_habstr[4,"AIC"]<-AIC(mod_list_surv1[[4]])
  
  mod_tab_habstr[5,"Model_Name"]<-"gerardii pres/abs"
  mod_tab_habstr[5,"Function"]<-deparse1(mod_list_surv1[[5]]$formula)
  mod_tab_habstr[5,"AIC"]<-AIC(mod_list_surv1[[5]])
  
  mod_tab_habstr[6,"Model_Name"]<-"gerardii pct"
  mod_tab_habstr[6,"Function"]<-deparse1(mod_list_surv1[[6]]$formula)
  mod_tab_habstr[6,"AIC"]<-AIC(mod_list_surv1[[6]])
  
  mod_tab_habstr[7,"Model_Name"]<-"patens pres/abs"
  mod_tab_habstr[7,"Function"]<-deparse1(mod_list_surv1[[7]]$formula)
  mod_tab_habstr[7,"AIC"]<-AIC(mod_list_surv1[[7]])
  
  mod_tab_habstr[8,"Model_Name"]<-"patens pct"
  mod_tab_habstr[8,"Function"]<-deparse1(mod_list_surv1[[8]]$formula)
  mod_tab_habstr[8,"AIC"]<-AIC(mod_list_surv1[[8]])
  
  mod_tab_habstr[9,"Model_Name"]<-"alt short pres/abs"
  mod_tab_habstr[9,"Function"]<-deparse1(mod_list_surv1[[9]]$formula)
  mod_tab_habstr[9,"AIC"]<-AIC(mod_list_surv1[[9]])
  
  mod_tab_habstr[10,"Model_Name"]<-"alt short pct"
  mod_tab_habstr[10,"Function"]<-deparse1(mod_list_surv1[[10]]$formula)
  mod_tab_habstr[10,"AIC"]<-AIC(mod_list_surv1[[10]])
  
  mod_tab_habstr[11,"Model_Name"]<-"phrag pres/abs"
  mod_tab_habstr[11,"Function"]<-deparse1(mod_list_surv1[[11]]$formula)
  mod_tab_habstr[11,"AIC"]<-AIC(mod_list_surv1[[11]])
  
  mod_tab_habstr[12,"Model_Name"]<-"phrag pct"
  mod_tab_habstr[12,"Function"]<-deparse1(mod_list_surv1[[12]]$formula)
  mod_tab_habstr[12,"AIC"]<-AIC(mod_list_surv1[[12]])
  
  mod_tab_habstr[13,"Model_Name"]<-"high marsh pres/abs"
  mod_tab_habstr[13,"Function"]<-deparse1(mod_list_surv1[[13]]$formula)
  mod_tab_habstr[13,"AIC"]<-AIC(mod_list_surv1[[13]])
  
  mod_tab_habstr[14,"Model_Name"]<-"high marsh pct"
  mod_tab_habstr[14,"Function"]<-deparse1(mod_list_surv1[[14]]$formula)
  mod_tab_habstr[14,"AIC"]<-AIC(mod_list_surv1[[14]])
  
  mod_tab_habstr[15,"Model_Name"]<-"brackish border pres/abs"
  mod_tab_habstr[15,"Function"]<-deparse1(mod_list_surv1[[15]]$formula)
  mod_tab_habstr[15,"AIC"]<-AIC(mod_list_surv1[[15]])
  
  mod_tab_habstr[16,"Model_Name"]<-"brackish border pct"
  mod_tab_habstr[16,"Function"]<-deparse1(mod_list_surv1[[16]]$formula)
  mod_tab_habstr[16,"AIC"]<-AIC(mod_list_surv1[[16]])
  
  mod_tab_habstr[17,"Model_Name"]<-"upland pres/abs"
  mod_tab_habstr[17,"Function"]<-deparse1(mod_list_surv1[[17]]$formula)
  mod_tab_habstr[17,"AIC"]<-AIC(mod_list_surv1[[17]])
  
  mod_tab_habstr[18,"Model_Name"]<-"upland pct"
  mod_tab_habstr[18,"Function"]<-deparse1(mod_list_surv1[[18]]$formula)
  mod_tab_habstr[18,"AIC"]<-AIC(mod_list_surv1[[18]])
  
  mod_tab_habstr[19,"Model_Name"]<-"water pres/abs"
  mod_tab_habstr[19,"Function"]<-deparse1(mod_list_surv1[[19]]$formula)
  mod_tab_habstr[19,"AIC"]<-AIC(mod_list_surv1[[19]])
  
  mod_tab_habstr[20,"Model_Name"]<-"water pct"
  mod_tab_habstr[20,"Function"]<-deparse1(mod_list_surv1[[20]]$formula)
  mod_tab_habstr[20,"AIC"]<-AIC(mod_list_surv1[[20]])
  
  mod_tab_habstr[21,"Model_Name"]<-"low marsh pres/abs"
  mod_tab_habstr[21,"Function"]<-deparse1(mod_list_surv1[[21]]$formula)
  mod_tab_habstr[21,"AIC"]<-AIC(mod_list_surv1[[21]])
  
  mod_tab_habstr[22,"Model_Name"]<-"low marsh pct"
  mod_tab_habstr[22,"Function"]<-deparse1(mod_list_surv1[[22]]$formula)
  mod_tab_habstr[22,"AIC"]<-AIC(mod_list_surv1[[22]])
  
  mod_tab_habstr[23,"Model_Name"]<-"saltmarsh border pres/abs"
  mod_tab_habstr[23,"Function"]<-deparse1(mod_list_surv1[[22]]$formula)
  mod_tab_habstr[23,"AIC"]<-AIC(mod_list_surv1[[23]])
  
  mod_tab_habstr[24,"Model_Name"]<-"saltmarsh border pct"
  mod_tab_habstr[24,"Function"]<-deparse1(mod_list_surv1[[24]]$formula)
  mod_tab_habstr[24,"AIC"]<-AIC(mod_list_surv1[[24]])
  
  mod_tab_habstr[25,"Model_Name"]<-"trees pres/abs"
  mod_tab_habstr[25,"Function"]<-deparse1(mod_list_surv1[[25]]$formula)
  mod_tab_habstr[25,"AIC"]<-AIC(mod_list_surv1[[25]])
  
  mod_tab_habstr[26,"Model_Name"]<-"trees pct"
  mod_tab_habstr[26,"Function"]<-deparse1(mod_list_surv1[[26]]$formula)
  mod_tab_habstr[26,"AIC"]<-AIC(mod_list_surv1[[26]])
  
  mod_tab_habstr[27,"Model_Name"]<-"simplified model with Lat and Lat^2"
  mod_tab_habstr[27,"Function"]<-deparse1(mod_list_surv1[[27]]$formula)
  mod_tab_habstr[27,"AIC"]<-AIC(mod_list_surv1[[27]])
  
  mod_tab_habstr[28,"Model_Name"]<-"global model"
  mod_tab_habstr[28,"Function"]<-deparse1(mod_list_surv1[[28]]$formula)
  mod_tab_habstr[28,"AIC"]<-AIC(mod_list_surv1[[28]])
  
# mod_tab_habstr[29,"Model_Name"]<-"simplified model with Lat and Lat^2"
#  mod_tab_habstr[29,"Function"]<-deparse1(mod_list_surv1[[29]]$formula)
#  mod_tab_habstr[29,"AIC"]<-AIC(mod_list_surv1[[29]])
  mod_tab_habstr<-group_by(mod_tab_habstr,Response)%>%mutate(dAIC=AIC-min(AIC))%>%
    ungroup()%>%
    arrange(Response,Model_Comparison,dAIC)
  
  
  # compare models
  mod_tab<-mod_tab_habstr%>%
    group_by(Response)%>%
    mutate(ΔAIC_within_Model_Comparison=round(dAIC,1),
           Overall_ΔAIC_across_Comparisons=round(AIC-min(AIC),1))%>%
    arrange(Response,Model_Comparison,dAIC)%>%
    ungroup()%>%
    dplyr::select(-dAIC) 
  mod_tab
  
  write.csv(mod_tab,paste0(path_out,"Final_outputs/Model_Results/",speciesnames[s],"_rapid_veg_survival_model_selection_table",".csv"), row.names = F)  
#  3.2 for Nesting Presence
  ## Mod List
  # Select model 
  #------------------------------------------
  ## Start a list of potential models
  mod_list_pres1<-list()
  # How much variation does each variable explain relative to each other?
  #alt_tall
  mod_list_pres1[[1]]<-glmmTMB(presence~alt_tall_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[2]]<-glmmTMB(presence~alt_tall_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #distichlis
  mod_list_pres1[[3]]<-glmmTMB(presence~distichlis_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[4]]<-glmmTMB(presence~distichlis_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #gerardii
  mod_list_pres1[[5]]<-glmmTMB(presence~gerardii_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[6]]<-glmmTMB(presence~gerardii_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #patens
  mod_list_pres1[[7]]<-glmmTMB(presence~patens_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[8]]<-glmmTMB(presence~patens_pct+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #alt_short
  mod_list_pres1[[9]]<-glmmTMB(presence~alt_short_pres+Lat+I(Lat^2), 
                               data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[10]]<-glmmTMB(presence~alt_short_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #phrag
  mod_list_pres1[[11]]<-glmmTMB(presence~phrag_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[12]]<-glmmTMB(presence~phrag_pct+Lat+I(Lat^2),
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #HIGHMARSH
  mod_list_pres1[[13]]<-glmmTMB(presence~high_marsh_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[14]]<-glmmTMB(presence~high_marsh_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Barckish Border
  mod_list_pres1[[15]]<-glmmTMB(presence~brackish_border_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[16]]<-glmmTMB(presence~brackish_border_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Upland
  
  mod_list_pres1[[17]]<-glmmTMB(presence~upland_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[18]]<-glmmTMB(presence~upland_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Trees
  mod_list_pres1[[19]]<-glmmTMB(presence~trees_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[20]]<-glmmTMB(presence~trees_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #Water
  mod_list_pres1[[21]]<-glmmTMB(presence~water_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[22]]<-glmmTMB(presence~water_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #LowMarsh
  mod_list_pres1[[23]]<-glmmTMB(presence~low_marsh_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[24]]<-glmmTMB(presence~low_marsh_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  #SaltmarshBorder
  mod_list_pres1[[25]]<-glmmTMB(presence~saltmarsh_border_pres+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  mod_list_pres1[[26]]<-glmmTMB(presence~saltmarsh_border_pct+Lat+I(Lat^2), 
                                data=rapid_dat10.model,family=ordbeta(link = "logit")
  )
  
  mod_list_pres1[[27]]<- simplified_pres_pct
  mod_list_pres1[[28]]<- fullmod_pres_pct
  #mod_list_pres1[[29]]<- simplified_pres_pct_Lat2
  
  # store model selection factors in a table
  mod_tab_habstr<-data.frame(Response=rep("Presence",length(mod_list_pres1)),
                             Model_Comparison=rep("Veg Species Occurrence",length(mod_list_pres1)),
                             Model_Name=rep(NA,length(mod_list_pres1)),
                             Function=rep(NA,length(mod_list_pres1)),
                             AIC=rep(NA,length(mod_list_pres1)),
                             dAIC=rep(NA,length(mod_list_pres1)))
  
  mod_tab_habstr[1,"Model_Name"]<-"alt tall pres/abs"
  mod_tab_habstr[1,"Function"]<-deparse1(mod_list_pres1[[1]]$formula)
  mod_tab_habstr[1,"AIC"]<-AIC(mod_list_pres1[[1]])
  
  mod_tab_habstr[2,"Model_Name"]<-"alt tall pct"
  mod_tab_habstr[2,"Function"]<-deparse1(mod_list_pres1[[2]]$formula)
  mod_tab_habstr[2,"AIC"]<-AIC(mod_list_pres1[[2]])
  
  mod_tab_habstr[3,"Model_Name"]<-"distichlis pres/abs"
  mod_tab_habstr[3,"Function"]<-deparse1(mod_list_pres1[[3]]$formula)
  mod_tab_habstr[3,"AIC"]<-AIC(mod_list_pres1[[3]])
  
  mod_tab_habstr[4,"Model_Name"]<-"distichlis pct"
  mod_tab_habstr[4,"Function"]<-deparse1(mod_list_pres1[[4]]$formula)
  mod_tab_habstr[4,"AIC"]<-AIC(mod_list_pres1[[4]])
  
  mod_tab_habstr[5,"Model_Name"]<-"gerardii pres/abs"
  mod_tab_habstr[5,"Function"]<-deparse1(mod_list_pres1[[5]]$formula)
  mod_tab_habstr[5,"AIC"]<-AIC(mod_list_pres1[[5]])
  
  mod_tab_habstr[6,"Model_Name"]<-"gerardii pct"
  mod_tab_habstr[6,"Function"]<-deparse1(mod_list_pres1[[6]]$formula)
  mod_tab_habstr[6,"AIC"]<-AIC(mod_list_pres1[[6]])
  
  mod_tab_habstr[7,"Model_Name"]<-"patens pres/abs"
  mod_tab_habstr[7,"Function"]<-deparse1(mod_list_pres1[[7]]$formula)
  mod_tab_habstr[7,"AIC"]<-AIC(mod_list_pres1[[7]])
  
  mod_tab_habstr[8,"Model_Name"]<-"patens pct"
  mod_tab_habstr[8,"Function"]<-deparse1(mod_list_pres1[[8]]$formula)
  mod_tab_habstr[8,"AIC"]<-AIC(mod_list_pres1[[8]])
  
  mod_tab_habstr[9,"Model_Name"]<-"alt short pres/abs"
  mod_tab_habstr[9,"Function"]<-deparse1(mod_list_pres1[[9]]$formula)
  mod_tab_habstr[9,"AIC"]<-AIC(mod_list_pres1[[9]])
  
  mod_tab_habstr[10,"Model_Name"]<-"alt short pct"
  mod_tab_habstr[10,"Function"]<-deparse1(mod_list_pres1[[10]]$formula)
  mod_tab_habstr[10,"AIC"]<-AIC(mod_list_pres1[[10]])
  
  mod_tab_habstr[11,"Model_Name"]<-"phrag pres/abs"
  mod_tab_habstr[11,"Function"]<-deparse1(mod_list_pres1[[11]]$formula)
  mod_tab_habstr[11,"AIC"]<-AIC(mod_list_pres1[[11]])
  
  mod_tab_habstr[12,"Model_Name"]<-"phrag pct"
  mod_tab_habstr[12,"Function"]<-deparse1(mod_list_pres1[[12]]$formula)
  mod_tab_habstr[12,"AIC"]<-AIC(mod_list_pres1[[12]])
  
  mod_tab_habstr[13,"Model_Name"]<-"high marsh pres/abs"
  mod_tab_habstr[13,"Function"]<-deparse1(mod_list_pres1[[13]]$formula)
  mod_tab_habstr[13,"AIC"]<-AIC(mod_list_pres1[[13]])
  
  mod_tab_habstr[14,"Model_Name"]<-"high marsh pct"
  mod_tab_habstr[14,"Function"]<-deparse1(mod_list_pres1[[14]]$formula)
  mod_tab_habstr[14,"AIC"]<-AIC(mod_list_pres1[[14]])
  
  mod_tab_habstr[15,"Model_Name"]<-"brackish border pres/abs"
  mod_tab_habstr[15,"Function"]<-deparse1(mod_list_pres1[[15]]$formula)
  mod_tab_habstr[15,"AIC"]<-AIC(mod_list_pres1[[15]])
  
  mod_tab_habstr[16,"Model_Name"]<-"brackish border pct"
  mod_tab_habstr[16,"Function"]<-deparse1(mod_list_pres1[[16]]$formula)
  mod_tab_habstr[16,"AIC"]<-AIC(mod_list_pres1[[16]])
  
  mod_tab_habstr[17,"Model_Name"]<-"upland pres/abs"
  mod_tab_habstr[17,"Function"]<-deparse1(mod_list_pres1[[17]]$formula)
  mod_tab_habstr[17,"AIC"]<-AIC(mod_list_pres1[[17]])
  
  mod_tab_habstr[18,"Model_Name"]<-"upland pct"
  mod_tab_habstr[18,"Function"]<-deparse1(mod_list_pres1[[18]]$formula)
  mod_tab_habstr[18,"AIC"]<-AIC(mod_list_pres1[[18]])
  
  mod_tab_habstr[19,"Model_Name"]<-"trees pres/abs"
  mod_tab_habstr[19,"Function"]<-deparse1(mod_list_pres1[[19]]$formula)
  mod_tab_habstr[19,"AIC"]<-AIC(mod_list_pres1[[19]])
  
  mod_tab_habstr[20,"Model_Name"]<-"trees pct"
  mod_tab_habstr[20,"Function"]<-deparse1(mod_list_pres1[[20]]$formula)
  mod_tab_habstr[20,"AIC"]<-AIC(mod_list_pres1[[20]])
  
  mod_tab_habstr[21,"Model_Name"]<-"water pres/abs"
  mod_tab_habstr[21,"Function"]<-deparse1(mod_list_pres1[[21]]$formula)
  mod_tab_habstr[21,"AIC"]<-AIC(mod_list_pres1[[21]])
  
  mod_tab_habstr[22,"Model_Name"]<-"water pct"
  mod_tab_habstr[22,"Function"]<-deparse1(mod_list_pres1[[22]]$formula)
  mod_tab_habstr[22,"AIC"]<-AIC(mod_list_pres1[[22]])
  
  mod_tab_habstr[23,"Model_Name"]<-"low marsh pres/abs"
  mod_tab_habstr[23,"Function"]<-deparse1(mod_list_pres1[[23]]$formula)
  mod_tab_habstr[23,"AIC"]<-AIC(mod_list_pres1[[23]])
  
  mod_tab_habstr[24,"Model_Name"]<-"low marsh pct"
  mod_tab_habstr[24,"Function"]<-deparse1(mod_list_pres1[[24]]$formula)
  mod_tab_habstr[24,"AIC"]<-AIC(mod_list_pres1[[24]])
  
  mod_tab_habstr[25,"Model_Name"]<-"saltmarsh border pres/abs"
  mod_tab_habstr[25,"Function"]<-deparse1(mod_list_pres1[[25]]$formula)
  mod_tab_habstr[25,"AIC"]<-AIC(mod_list_pres1[[25]])
  
  mod_tab_habstr[26,"Model_Name"]<-"saltmarsh border pct"
  mod_tab_habstr[26,"Function"]<-deparse1(mod_list_pres1[[26]]$formula)
  mod_tab_habstr[26,"AIC"]<-AIC(mod_list_pres1[[26]])
  
  mod_tab_habstr[27,"Model_Name"]<-"simplified model"
  mod_tab_habstr[27,"Function"]<-deparse1(mod_list_pres1[[27]]$formula)
  mod_tab_habstr[27,"AIC"]<-AIC(mod_list_pres1[[27]])
  
  mod_tab_habstr[28,"Model_Name"]<-"global model"
  mod_tab_habstr[28,"Function"]<-deparse1(mod_list_pres1[[28]]$formula)
  mod_tab_habstr[28,"AIC"]<-AIC(mod_list_pres1[[28]])
  
  #mod_tab_habstr[29,"Model_Name"]<-"simplified model with Lat and Lat^2"
  #mod_tab_habstr[29,"Function"]<-deparse1(mod_list_pres1[[29]]$formula)
  #mod_tab_habstr[29,"AIC"]<-AIC(mod_list_pres1[[29]])
  mod_tab_habstr<-group_by(mod_tab_habstr,Response)%>%mutate(dAIC=AIC-min(AIC))%>%
    ungroup()%>%
    arrange(Response,Model_Comparison,dAIC)
  
  
  # compare models
  mod_tab<-mod_tab_habstr%>%
    group_by(Response)%>%
    mutate(ΔAIC_within_Model_Comparison=round(dAIC,1),
           Overall_ΔAIC_across_Comparisons=round(AIC-min(AIC),1))%>%
    arrange(Response,Model_Comparison,dAIC)%>%
    ungroup()%>%
    dplyr::select(-dAIC) 
  mod_tab
  
  write.csv(mod_tab,paste0(path_out,"Final_outputs/Model_Results/",speciesnames[s],"_rapid_veg_presence_model_selection_table",".csv"), row.names = F)

  simp_surv_coef<-as.data.frame(fixef(simplified_surv_pct)$cond)
  colnames(simp_surv_coef)<-"coef"
  simp_surv_coef$plogis<-apply(simp_surv_coef,1,FUN = function(x){plogis(x[1]+simp_surv_coef$coef[1])-plogis(simp_surv_coef$coef[1])})
  write.csv(simp_surv_coef,paste0(path_out,"Final_outputs/Model_Results/",speciesnames[s],"_simp_surv_coef_table",".csv"), row.names = F)
  
  simp_pres_coef<-as.data.frame(fixef(simplified_pres_pct)$cond)
  colnames(simp_pres_coef)<-"coef"
  simp_pres_coef$plogis<-apply(simp_pres_coef,1,FUN = function(x){plogis(x[1]+simp_pres_coef$coef[1])-plogis(simp_pres_coef$coef[1])})
  write.csv(simp_pres_coef,paste0(path_out,"Final_outputs/Model_Results/",speciesnames[s],"_simp_pres_coef_table",".csv"), row.names = F)
  
    }
