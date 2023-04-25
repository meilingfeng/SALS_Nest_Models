library(tidyverse)
library(car)#vif

### Set up
# -------------------------------------------
source("05a_Set_Model_Parameters.R")




## Start a list of potential models

#Use proportion of high marsh over high marsh at nest - account for error in nest location and resolution uncertainty in veg layer

# 1. Habitat characteristics not captured by High Marsh

# compare a model using just proportion of high marsh to models just using NDVI, Entropy, UVVR
# compare models using additive combinations of each: highmarsh+NDVI. Highmarsh+UVVR. Highmarsh+PC, 
# Compare model using combinations of non-highmarsh: NDVI+UVVR, NDVI+PCA, PCA+UVVR, PCA+UVVR+NDVI
mod_list_pres1<-list()
mod_list_surv1<-list()



# How much variation does each variable explain relative to each other?
# NDVI
mod_list_pres1[[1]]<-glm(y~ndvi, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[1]]<-glm(y~ndvi, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# Entropy
mod_list_pres1[[2]]<-glm(y~ent_txt, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[2]]<-glm(y~ent_txt, 
                        data=surv_dat,
                        family = binomial(link="logit"))
# UVVR
mod_list_pres1[[3]]<-glm(y~uvvr_mean, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[3]]<-glm(y~uvvr_mean, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# HIGH MARSH
mod_list_pres1[[4]]<-glm(y~HIMARSH, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv1[[4]]<-glm(y~HIMARSH, 
                        data=surv_dat,
                        family = binomial(link="logit"))

#Does adding NDVI or UVVR or Entropy improve high marsh?
# HIGH MARSH +NDVI
mod_list_pres1[[5]]<-glm(y~HIMARSH+ndvi, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[5]]<-glm(y~HIMARSH+ndvi, 
                         data=surv_dat,
                         family = binomial(link="logit"))
# HIGH MARSH + UVVR
mod_list_pres1[[6]]<-glm(y~HIMARSH+uvvr_mean, 
                          data=pres_dat,
                          family = binomial(link="logit"))
mod_list_surv1[[6]]<-glm(y~HIMARSH+uvvr_mean, 
                          data=surv_dat,
                          family = binomial(link="logit"))
# HIGH MARSH + Entropy
mod_list_pres1[[7]]<-glm(y~HIMARSH+ent_txt, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[7]]<-glm(y~HIMARSH+ent_txt, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# ALL
mod_list_pres1[[8]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv1[[8]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH,
                         data=surv_dat,
                         family = binomial(link="logit"))
vif_pres1<-vif(mod_list_pres1[[8]])
vif_surv1<-vif(mod_list_surv1[[8]])


# store model selection factors in a table
mod_tab_habstr<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres1)),
                           class=rep("Additional Habitat Characteristics",2*length(mod_list_pres1)),
                           name=rep(NA,2*length(mod_list_pres1)),
                           fun=rep(NA,2*length(mod_list_pres1)),
                           max_VIF=rep(NA,2*length(mod_list_pres1)),
                           VIF_var=rep(NA,2*length(mod_list_pres1)),
                           AIC=rep(NA,2*length(mod_list_pres1)),
                           dAIC=rep(NA,2*length(mod_list_pres1)))

mod_tab_habstr[1,"name"]<-"Vegetation Health (NDVI)"
mod_tab_habstr[1,"fun"]<-deparse1(mod_list_pres1[[1]]$formula)
mod_tab_habstr[1,"AIC"]<-mod_list_pres1[[1]]$aic

mod_tab_habstr[2,"name"]<-"Vegetation Health (NDVI)"
mod_tab_habstr[2,"fun"]<-deparse1(mod_list_surv1[[1]]$formula)
mod_tab_habstr[2,"AIC"]<-mod_list_surv1[[1]]$aic

mod_tab_habstr[3,"name"]<-"Habitat Dissimilarity (Entropy)"
mod_tab_habstr[3,"fun"]<-deparse1(mod_list_pres1[[2]]$formula)
mod_tab_habstr[3,"AIC"]<-mod_list_pres1[[2]]$aic

mod_tab_habstr[4,"name"]<-"Habitat Dissimilarity (Entropy)"
mod_tab_habstr[4,"fun"]<-deparse1(mod_list_surv1[[2]]$formula)
mod_tab_habstr[4,"AIC"]<-mod_list_surv1[[2]]$aic

mod_tab_habstr[5,"name"]<-"Marsh Resilience (UVVR)"
mod_tab_habstr[5,"fun"]<-deparse1(mod_list_pres1[[3]]$formula)
mod_tab_habstr[5,"AIC"]<-mod_list_pres1[[3]]$aic

mod_tab_habstr[6,"name"]<-"Marsh Resilience (UVVR)"
mod_tab_habstr[6,"fun"]<-deparse1(mod_list_surv1[[3]]$formula)
mod_tab_habstr[6,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_habstr[7,"name"]<-"High Marsh"
mod_tab_habstr[7,"fun"]<-deparse1(mod_list_pres1[[4]]$formula)
mod_tab_habstr[7,"AIC"]<-mod_list_pres1[[4]]$aic

mod_tab_habstr[8,"name"]<-"High Marsh"
mod_tab_habstr[8,"fun"]<-deparse1(mod_list_surv1[[4]]$formula)
mod_tab_habstr[8,"AIC"]<-mod_list_surv1[[4]]$aic

mod_tab_habstr[9,"name"]<-"High Marsh + NDVI"
mod_tab_habstr[9,"fun"]<-deparse1(mod_list_pres1[[5]]$formula)
mod_tab_habstr[9,"AIC"]<-mod_list_pres1[[5]]$aic

mod_tab_habstr[10,"name"]<-"High Marsh + NDVI"
mod_tab_habstr[10,"fun"]<-deparse1(mod_list_surv1[[5]]$formula)
mod_tab_habstr[10,"AIC"]<-mod_list_surv1[[5]]$aic

mod_tab_habstr[11,"name"]<-"High Marsh + UVVR"
mod_tab_habstr[11,"fun"]<-deparse1(mod_list_pres1[[6]]$formula)
mod_tab_habstr[11,"AIC"]<-mod_list_pres1[[6]]$aic

mod_tab_habstr[12,"name"]<-"High Marsh + UVVR"
mod_tab_habstr[12,"fun"]<-deparse1(mod_list_surv1[[6]]$formula)
mod_tab_habstr[12,"AIC"]<-mod_list_surv1[[6]]$aic


mod_tab_habstr[13,"name"]<-"High Marsh + Entropy"
mod_tab_habstr[13,"fun"]<-deparse1(mod_list_pres1[[7]]$formula)
mod_tab_habstr[13,"AIC"]<-mod_list_pres1[[7]]$aic

mod_tab_habstr[14,"name"]<-"High Marsh + Entropy"
mod_tab_habstr[14,"fun"]<-deparse1(mod_list_surv1[[7]]$formula)
mod_tab_habstr[14,"AIC"]<-mod_list_surv1[[7]]$aic

mod_tab_habstr[15,"name"]<-"Full Additional Habitat"
mod_tab_habstr[15,"fun"]<-deparse1(mod_list_pres1[[8]]$formula)
mod_tab_habstr[15,"max_VIF"]<-max(vif_pres1)
mod_tab_habstr[15,"VIF_var"]<-names(vif_pres1[which(vif_pres1==max(vif_pres1))])
mod_tab_habstr[15,"AIC"]<-mod_list_pres1[[8]]$aic

mod_tab_habstr[16,"name"]<-"Full Additional Habitat"
mod_tab_habstr[16,"fun"]<-deparse1(mod_list_surv1[[8]]$formula)
mod_tab_habstr[16,"max_VIF"]<-max(vif_surv1)
mod_tab_habstr[16,"VIF_var"]<-names(vif_surv1[which(vif_surv1==max(vif_surv1))])
mod_tab_habstr[16,"AIC"]<-mod_list_surv1[[8]]$aic


mod_tab_habstr<-group_by(mod_tab_habstr,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,class,dAIC)





# 2. Is there any additional information in reflectance that has not been classified/or explained? - Does adding PCA improve UVVR, NDVI, and High marsh?

mod_list_pres2<-list()
mod_list_surv2<-list()

# PCA
mod_list_pres2[[1]]<-glm(y~pca, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv2[[1]]<-glm(y~pca, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# Highmarsh + PCA
mod_list_pres2[[2]]<-glm(y~HIMARSH+pca, 
                         data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv2[[2]]<-glm(y~HIMARSH+pca, 
                         data=surv_dat,
                         family = binomial(link="logit"))

# ALL
#mod_list_pres2[[3]]<-glm(as.formula(paste0(mod_tab_habstr[mod_tab_habstr$response=="Presence" & mod_tab_habstr$dAIC==0, 4],"+pca")), 
#                         data=pres_dat,
#                         family = binomial(link="logit"))
#mod_list_surv2[[3]]<-glm(as.formula(paste0(mod_tab_habstr[mod_tab_habstr$response=="Success" & mod_tab_habstr$dAIC==0, 4],"+pca")),
#                         data=surv_dat,
#                         family = binomial(link="logit"))


#Additional Habitat with Unclassified habitat
mod_list_pres2[[3]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+pca,data=pres_dat,
                         family = binomial(link="logit"))
mod_list_surv2[[3]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+pca,
                         data=surv_dat,
                         family = binomial(link="logit"))
vif_pres2<-vif(mod_list_pres2[[3]])
vif_surv2<-vif(mod_list_surv2[[3]])

# store model selection factors in a table
mod_tab_unclass<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres2)),
                           class=rep("Unclassified Habitat Characteristics",2*length(mod_list_pres2)),
                           name=rep(NA,2*length(mod_list_pres2)),
                           fun=rep(NA,2*length(mod_list_pres2)),
                           max_VIF=rep(NA,2*length(mod_list_pres2)),
                           VIF_var=rep(NA,2*length(mod_list_pres2)),
                           AIC=rep(NA,2*length(mod_list_pres2)),
                           dAIC=rep(NA,2*length(mod_list_pres2)))


mod_tab_unclass[1,"name"]<-"Unclassified Habitat"
mod_tab_unclass[1,"fun"]<-deparse1(mod_list_pres2[[1]]$formula)
mod_tab_unclass[1,"AIC"]<-mod_list_pres2[[1]]$aic

mod_tab_unclass[2,"name"]<-"Unclassified Habitat"
mod_tab_unclass[2,"fun"]<-deparse1(mod_list_surv2[[1]]$formula)
mod_tab_unclass[2,"AIC"]<-mod_list_surv2[[1]]$aic

mod_tab_unclass[3,"name"]<-"High Marsh + Unclassified Habitat"
mod_tab_unclass[3,"fun"]<-deparse1(mod_list_pres2[[2]]$formula)
mod_tab_unclass[3,"AIC"]<-mod_list_pres2[[2]]$aic

mod_tab_unclass[4,"name"]<-"High Marsh + Unclassified Habitat"
mod_tab_unclass[4,"fun"]<-deparse1(mod_list_surv2[[2]]$formula)
mod_tab_unclass[4,"AIC"]<-mod_list_surv2[[2]]$aic

mod_tab_unclass[5,"name"]<-"Additional Habitat + Unclassified Habitat"
mod_tab_unclass[5,"fun"]<-deparse1(mod_list_pres2[[3]]$formula)
mod_tab_unclass[5,"max_VIF"]<-max(vif_pres2)[1]
mod_tab_unclass[5,"VIF_var"]<-names(vif_pres2[which(vif_pres2==max(vif_pres2))])[1]
mod_tab_unclass[5,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_unclass[6,"name"]<-"Additional Habitat + Unclassified Habitat"
mod_tab_unclass[6,"fun"]<-deparse1(mod_list_surv2[[3]]$formula)
mod_tab_unclass[6,"max_VIF"]<-max(vif_surv2)[1]
mod_tab_unclass[6,"VIF_var"]<-names(vif_surv2[which(vif_surv2==max(vif_surv2))])[1]
mod_tab_unclass[6,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_unclass<-group_by(mod_tab_unclass,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,class,dAIC)





# 3. Quality within high marsh (interactions) - sparrow presence does not equal nesting presence
# within high marsh habitat, certain habitat qualities are more suitable for nesting
  # eg healthier vegetation (NDVI) grows taller and denser
  # compare highmarsh only, ndvi only, highmarsh + ndvi to high marsh*ndvi
mod_list_pres3<-list()
mod_list_surv3<-list()

#do the interaction models do better than the additive models?


# high marsh and PC reflect
mod_list_pres3[[1]]<-glm(y~HIMARSH*pca, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv3[[1]]<-glm(y~HIMARSH*pca, 
                        data=surv_dat,
                        family = binomial(link="logit"))


# store model selection factors in a table
mod_tab_int<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres3)),
                           class=rep("High Marsh Quality",2*length(mod_list_pres3)),
                           name=rep(NA,2*length(mod_list_pres3)),
                           fun=rep(NA,2*length(mod_list_pres3)),
                           max_VIF=rep(NA,2*length(mod_list_pres3)),
                           VIF_var=rep(NA,2*length(mod_list_pres3)),
                           AIC=rep(NA,2*length(mod_list_pres3)),
                           dAIC=rep(NA,2*length(mod_list_pres3)))

mod_tab_int[1,"name"]<-"High Marsh * Unclassified Habitat"
mod_tab_int[1,"fun"]<-deparse1(mod_list_pres3[[1]]$formula)
mod_tab_int[1,"AIC"]<-mod_list_pres3[[1]]$aic

mod_tab_int[2,"name"]<-"High Marsh * Unclassified Habitat"
mod_tab_int[2,"fun"]<-deparse1(mod_list_surv3[[1]]$formula)
mod_tab_int[2,"AIC"]<-mod_list_surv3[[1]]$aic


mod_tab_int<-group_by(mod_tab_int,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)





################################
# compare models
mod_tab<-rbind(mod_tab_habstr,mod_tab_int,mod_tab_unclass)%>%
  group_by(response)%>%
  mutate(overall_dAIC=AIC-min(AIC))%>%
  arrange(response,overall_dAIC)%>%
  ungroup()

# PCA more important for nest placement, high marsh presence more important for success
# Nest in high marsh (local habitat) more important than surrounding habitat for both
# Correlation more important than entropy, but not as important as proportion high marsh

#################################




# 4. Do these habitat relationships change across the breeding range? (eg more success at center of range?)
# Latitude explains variation in these relationships throughout the breeding range
# add latitude on its own and to the top models
mod_list_pres4<-list()
mod_list_surv4<-list()

# compare best models of habitat with latitude



# Just latitude
mod_list_pres4[[1]]<-glm(y~latitude, 
                        data=pres_dat,
                        family = binomial(link="logit"))
mod_list_surv4[[1]]<-glm(y~latitude, 
                        data=surv_dat,
                        family = binomial(link="logit"))

# latitude with high marsh
mod_list_pres4[[2]]<-glm(y~HIMARSH+latitude,
                         data=pres_dat,
                         family = binomial(link="logit"))

mod_list_surv4[[2]]<-glm(y~HIMARSH+latitude,
                         data=surv_dat,
                         family = binomial(link="logit"))

# latitude with additional habitat
mod_list_pres4[[3]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+latitude,
                         data=pres_dat,
                         family = binomial(link="logit"))

mod_list_surv4[[3]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+latitude,
                         data=surv_dat,
                         family = binomial(link="logit"))

# latitude unclassified habitat
mod_list_pres4[[4]]<-glm(y~HIMARSH+pca+latitude,
                         data=pres_dat,
                         family = binomial(link="logit"))

mod_list_surv4[[4]]<-glm(y~HIMARSH+pca+latitude,
                         data=surv_dat,
                         family = binomial(link="logit"))



# latitude with additional and unclassified habitat
mod_list_pres4[[5]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+pca+latitude,
                         data=pres_dat,
                         family = binomial(link="logit"))

mod_list_surv4[[5]]<-glm(y~ndvi+uvvr_mean+ent_txt+HIMARSH+pca+latitude,
                         data=surv_dat,
                         family = binomial(link="logit"))






# store model selection factors in a table
mod_tab_lat<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres4)),
                        class=rep("Geographic Variation",2*length(mod_list_pres4)),
                        name=rep(NA,2*length(mod_list_pres4)),
                        fun=rep(NA,2*length(mod_list_pres4)),
                        max_VIF=rep(NA,2*length(mod_list_pres4)),
                        VIF_var=rep(NA,2*length(mod_list_pres4)),
                        AIC=rep(NA,2*length(mod_list_pres4)),
                        dAIC=rep(NA,2*length(mod_list_pres4)))



mod_tab_lat[1,"name"]<-"Latitude"
mod_tab_lat[1,"fun"]<-deparse1(mod_list_pres4[[1]]$formula)
mod_tab_lat[1,"AIC"]<-mod_list_pres4[[1]]$aic

mod_tab_lat[2,"name"]<-"Latitude"
mod_tab_lat[2,"fun"]<-deparse1(mod_list_surv4[[1]]$formula)
mod_tab_lat[2,"AIC"]<-mod_list_surv4[[1]]$aic


mod_tab_lat[3,"name"]<-"High Marsh + Latitude"
mod_tab_lat[3,"fun"]<-deparse1(mod_list_pres4[[2]]$formula)
mod_tab_lat[3,"AIC"]<-mod_list_pres4[[2]]$aic

mod_tab_lat[4,"name"]<-"High Marsh + Latitude"
mod_tab_lat[4,"fun"]<-deparse1(mod_list_surv4[[2]]$formula)
mod_tab_lat[4,"AIC"]<-mod_list_surv4[[2]]$aic


mod_tab_lat[5,"name"]<-"Additional Habitat + Latitude"
mod_tab_lat[5,"fun"]<-deparse1(mod_list_pres4[[3]]$formula)
mod_tab_lat[5,"AIC"]<-mod_list_pres4[[3]]$aic

mod_tab_lat[6,"name"]<-"Additional Habitat + Latitude"
mod_tab_lat[6,"fun"]<-deparse1(mod_list_surv4[[3]]$formula)
mod_tab_lat[6,"AIC"]<-mod_list_surv4[[3]]$aic


mod_tab_lat[7,"name"]<-"Unclassified Habitat + Latitude"
mod_tab_lat[7,"fun"]<-deparse1(mod_list_pres4[[4]]$formula)
mod_tab_lat[7,"AIC"]<-mod_list_pres4[[4]]$aic

mod_tab_lat[8,"name"]<-"Unclassified Habitat + Latitude"
mod_tab_lat[8,"fun"]<-deparse1(mod_list_surv4[[4]]$formula)
mod_tab_lat[8,"AIC"]<-mod_list_surv4[[4]]$aic


mod_tab_lat[9,"name"]<-"Additional and Unclassified Habitat + Latitude"
mod_tab_lat[9,"fun"]<-deparse1(mod_list_pres4[[5]]$formula)
mod_tab_lat[9,"AIC"]<-mod_list_pres4[[5]]$aic

mod_tab_lat[10,"name"]<-"Additional and Unclassified Habitat + Latitude"
mod_tab_lat[10,"fun"]<-deparse1(mod_list_surv4[[5]]$formula)
mod_tab_lat[10,"AIC"]<-mod_list_surv4[[5]]$aic



mod_tab_lat<-group_by(mod_tab_lat,response)%>%
  mutate(dAIC=AIC-min(AIC),
         overall_dAIC=NA)%>%
  ungroup()%>%
  arrange(response,dAIC)




################################
# combine all tables
mod_tab<-rbind(mod_tab,mod_tab_lat)%>%
  group_by(response)%>%
  mutate(overall_dAIC=round(AIC-min(AIC)),
         dAIC=round(dAIC),
         AIC=round(AIC),
         max_VIF=round(max_VIF,2))%>%
  arrange(response,overall_dAIC)%>%
  ungroup()


mod_list_pres_final<-c(mod_list_pres1,mod_list_pres2,mod_list_pres3,mod_list_pres4)
mod_list_surv_final<-c(mod_list_surv1,mod_list_surv2,mod_list_surv3,mod_list_surv4)






# PCA more important for nest placement, high marsh presence more important for success
# Nest in high marsh (local habitat) more important than surrounding habitat for both
# Correlation more important than entropy, but not as important as proportion high marsh

#################################