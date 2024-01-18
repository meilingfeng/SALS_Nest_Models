# store model selection factors in a table
mod_tab_habstr<-data.frame(Response=rep(c("Presence","Success"),length(mod_list_pres1)),
                           Model_Name=rep(NA,2*length(mod_list_pres1)),
                           Function=rep(NA,2*length(mod_list_pres1)),
                           AIC=rep(NA,2*length(mod_list_pres1)),
                           dAIC=rep(NA,2*length(mod_list_pres1)))

mod_tab_habstr[1,"Model_Name"]<-"High Marsh Habitat"
mod_tab_habstr[1,"Function"]<-deparse1(mod_list_pres1[[1]]$formula)
mod_tab_habstr[1,"AIC"]<-mod_list_pres1[[1]]$aic

mod_tab_habstr[2,"Model_Name"]<-"High Marsh Habitat"
mod_tab_habstr[2,"Function"]<-deparse1(mod_list_surv1[[1]]$formula)
mod_tab_habstr[2,"AIC"]<-mod_list_surv1[[1]]$aic

mod_tab_habstr[3,"Model_Name"]<-"Additional Habitat Features"
mod_tab_habstr[3,"Function"]<-deparse1(mod_list_pres1[[2]]$formula)
mod_tab_habstr[3,"AIC"]<-mod_list_pres1[[2]]$aic

mod_tab_habstr[4,"Model_Name"]<-"Additional Habitat Features"
mod_tab_habstr[4,"Function"]<-deparse1(mod_list_surv1[[2]]$formula)
mod_tab_habstr[4,"AIC"]<-mod_list_surv1[[2]]$aic

mod_tab_habstr[5,"Model_Name"]<-"High Marsh + Additional Habitat Features"
mod_tab_habstr[5,"Function"]<-deparse1(mod_list_pres1[[3]]$formula)
mod_tab_habstr[5,"AIC"]<-mod_list_pres1[[3]]$aic

mod_tab_habstr[6,"Model_Name"]<-"High Marsh + Additional Habitat Features"
mod_tab_habstr[6,"Function"]<-deparse1(mod_list_surv1[[3]]$formula)
mod_tab_habstr[6,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_habstr[7,"Model_Name"]<-"Unclassified Habitat Features"
mod_tab_habstr[7,"Function"]<-deparse1(mod_list_pres1[[4]]$formula)
mod_tab_habstr[7,"AIC"]<-mod_list_pres1[[4]]$aic

mod_tab_habstr[8,"Model_Name"]<-"Unclassified Habitat Features"
mod_tab_habstr[8,"Function"]<-deparse1(mod_list_surv1[[4]]$formula)
mod_tab_habstr[8,"AIC"]<-mod_list_surv1[[4]]$aic

mod_tab_habstr[9,"Model_Name"]<-"Elevation"
mod_tab_habstr[9,"Function"]<-deparse1(mod_list_pres1[[5]]$formula)
mod_tab_habstr[9,"AIC"]<-mod_list_pres1[[5]]$aic

mod_tab_habstr[10,"Model_Name"]<-"Elevation"
mod_tab_habstr[10,"Function"]<-deparse1(mod_list_surv1[[5]]$formula)
mod_tab_habstr[10,"AIC"]<-mod_list_surv1[[5]]$aic

mod_tab_habstr[11,"Model_Name"]<-"Vegetation Features"
mod_tab_habstr[11,"Function"]<-deparse1(mod_list_pres1[[6]]$formula)
mod_tab_habstr[11,"AIC"]<-mod_list_pres1[[6]]$aic

mod_tab_habstr[12,"Model_Name"]<-"Vegetation Features"
mod_tab_habstr[12,"Function"]<-deparse1(mod_list_surv1[[6]]$formula)
mod_tab_habstr[12,"AIC"]<-mod_list_surv1[[6]]$aic

mod_tab_habstr[13,"Model_Name"]<-"High Marsh Quality: Reflectance"
mod_tab_habstr[13,"Function"]<-deparse1(mod_list_pres1[[7]]$formula)
mod_tab_habstr[13,"AIC"]<-mod_list_pres1[[7]]$aic

mod_tab_habstr[14,"Model_Name"]<-"High Marsh Quality: Reflectance"
mod_tab_habstr[14,"Function"]<-deparse1(mod_list_surv1[[7]]$formula)
mod_tab_habstr[14,"AIC"]<-mod_list_surv1[[7]]$aic

mod_tab_habstr[15,"Model_Name"]<-"High Marsh Quality: NDVI"
mod_tab_habstr[15,"Function"]<-deparse1(mod_list_pres1[[8]]$formula)
mod_tab_habstr[15,"AIC"]<-mod_list_pres1[[8]]$aic

mod_tab_habstr[16,"Model_Name"]<-"High Marsh Quality: NDVI"
mod_tab_habstr[16,"Function"]<-deparse1(mod_list_surv1[[8]]$formula)
mod_tab_habstr[16,"AIC"]<-mod_list_surv1[[8]]$aic


mod_tab<-group_by(mod_tab_habstr,Response)%>%mutate(ΔAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(Response,ΔAIC)