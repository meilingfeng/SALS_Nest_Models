library(tidyverse)
library(car)#vif

### Set up
# -------------------------------------------
source("05a_Set_Model_Parameters.R")




## Start a list of potential models


# 1. Habitat Structure: Does NDVI and spectral data explain additional information about important habitat characteristics beyond structure inherent in plant species form and species composition?
# Compare NDVI, PCA, and High Marsh vegetation classification
mod_list_pres1<-list()
mod_list_surv1<-list()
# ALL
mod_list_pres1[[1]]<-glm(y~ndvi+pca+Highmarsh, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv1[[1]]<-glm(y~ndvi+pca+Highmarsh, 
                        data=surv_train,
                        family = binomial(link="logit"))
vif_pres<-vif(mod_list_pres1[[1]])
vif_surv<-vif(mod_list_surv1[[1]])


# NDVI
mod_list_pres1[[2]]<-glm(y~ndvi, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv1[[2]]<-glm(y~ndvi, 
                        data=surv_train,
                        family = binomial(link="logit"))

# PCA
mod_list_pres1[[3]]<-glm(y~pca, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv1[[3]]<-glm(y~pca, 
                        data=surv_train,
                        family = binomial(link="logit"))

# HIGH MARSH
mod_list_pres1[[4]]<-glm(y~Highmarsh, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv1[[4]]<-glm(y~Highmarsh, 
                        data=surv_train,
                        family = binomial(link="logit"))

# store model selection factors in a table
mod_tab_habstr<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres1)),
                           class=rep("Habitat Structure",2*length(mod_list_pres1)),
                           name=rep(NA,2*length(mod_list_pres1)),
                           fun=rep(NA,2*length(mod_list_pres1)),
                           max_VIF=rep(NA,2*length(mod_list_pres1)),
                           VIF_var=rep(NA,2*length(mod_list_pres1)),
                           AIC=rep(NA,2*length(mod_list_pres1)),
                           dAIC=rep(NA,2*length(mod_list_pres1)),
                           AUC=rep(NA,2*length(mod_list_pres1)))

mod_tab_habstr[1,"name"]<-"All"
mod_tab_habstr[1,"fun"]<-deparse1(mod_list_pres1[[1]]$formula)
mod_tab_habstr[1,"max_VIF"]<-max(vif_pres)
mod_tab_habstr[1,"VIF_var"]<-names(vif_pres[which(vif_pres==max(vif_pres))])
mod_tab_habstr[1,"AIC"]<-mod_list_pres1[[1]]$aic

mod_tab_habstr[2,"name"]<-"All"
mod_tab_habstr[2,"fun"]<-deparse1(mod_list_surv1[[1]]$formula)
mod_tab_habstr[2,"max_VIF"]<-max(vif_surv)
mod_tab_habstr[2,"VIF_var"]<-names(vif_surv[which(vif_surv==max(vif_surv))])
mod_tab_habstr[2,"AIC"]<-mod_list_surv1[[1]]$aic

mod_tab_habstr[3,"name"]<-"NDVI"
mod_tab_habstr[3,"fun"]<-deparse1(mod_list_pres1[[2]]$formula)
mod_tab_habstr[3,"AIC"]<-mod_list_pres1[[2]]$aic

mod_tab_habstr[4,"name"]<-"NDVI"
mod_tab_habstr[4,"fun"]<-deparse1(mod_list_surv1[[2]]$formula)
mod_tab_habstr[4,"AIC"]<-mod_list_surv1[[2]]$aic

mod_tab_habstr[5,"name"]<-"Reflectance PC"
mod_tab_habstr[5,"fun"]<-deparse1(mod_list_pres1[[3]]$formula)
mod_tab_habstr[5,"AIC"]<-mod_list_pres1[[3]]$aic

mod_tab_habstr[6,"name"]<-"Reflectance PC"
mod_tab_habstr[6,"fun"]<-deparse1(mod_list_surv1[[3]]$formula)
mod_tab_habstr[6,"AIC"]<-mod_list_surv1[[3]]$aic

mod_tab_habstr[7,"name"]<-"High Marsh"
mod_tab_habstr[7,"fun"]<-deparse1(mod_list_pres1[[4]]$formula)
mod_tab_habstr[7,"AIC"]<-mod_list_pres1[[4]]$aic

mod_tab_habstr[8,"name"]<-"High Marsh"
mod_tab_habstr[8,"fun"]<-deparse1(mod_list_surv1[[4]]$formula)
mod_tab_habstr[8,"AIC"]<-mod_list_surv1[[4]]$aic


mod_tab_habstr<-group_by(mod_tab_habstr,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)




# 2. Scale: Does habitat matter more at the nest location or surrounding habitat? 
# Compare high marsh at nest to high marsh prevalence in surrounding buffer
mod_list_pres2<-list()
mod_list_surv2<-list()

#ALL
mod_list_pres2[[1]]<-glm(y~HIMARSH+Highmarsh, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv2[[1]]<-glm(y~HIMARSH+Highmarsh, 
                        data=surv_train,
                        family = binomial(link="logit"))
vif_pres<-vif(mod_list_pres2[[1]])
vif_surv<-vif(mod_list_surv2[[1]])


# Proportion high marsh
mod_list_pres2[[2]]<-glm(y~HIMARSH, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv2[[2]]<-glm(y~HIMARSH, 
                        data=surv_train,
                        family = binomial(link="logit"))

# Nest in high marsh
mod_list_pres2[[3]]<-glm(y~Highmarsh, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv2[[3]]<-glm(y~Highmarsh, 
                        data=surv_train,
                        family = binomial(link="logit"))





# store model selection factors in a table
mod_tab_scale<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres2)),
                           class=rep("Scale",2*length(mod_list_pres2)),
                           name=rep(NA,2*length(mod_list_pres2)),
                           fun=rep(NA,2*length(mod_list_pres2)),
                           max_VIF=rep(NA,2*length(mod_list_pres2)),
                           VIF_var=rep(NA,2*length(mod_list_pres2)),
                           AIC=rep(NA,2*length(mod_list_pres2)),
                           dAIC=rep(NA,2*length(mod_list_pres2)),
                           AUC=rep(NA,2*length(mod_list_pres2)))

mod_tab_scale[1,"name"]<-"All"
mod_tab_scale[1,"fun"]<-deparse1(mod_list_pres2[[1]]$formula)
mod_tab_scale[1,"max_VIF"]<-max(vif_pres)
mod_tab_scale[1,"VIF_var"]<-names(vif_pres[which(vif_pres==max(vif_pres))])[1]
mod_tab_scale[1,"AIC"]<-mod_list_pres2[[1]]$aic

mod_tab_scale[2,"name"]<-"All"
mod_tab_scale[2,"fun"]<-deparse1(mod_list_surv2[[1]]$formula)
mod_tab_scale[2,"max_VIF"]<-max(vif_surv)
mod_tab_scale[2,"VIF_var"]<-names(vif_surv[which(vif_surv==max(vif_surv))])[1]
mod_tab_scale[2,"AIC"]<-mod_list_surv2[[1]]$aic

mod_tab_scale[3,"name"]<-"Proportion High Marsh"
mod_tab_scale[3,"fun"]<-deparse1(mod_list_pres2[[2]]$formula)
mod_tab_scale[3,"AIC"]<-mod_list_pres2[[2]]$aic

mod_tab_scale[4,"name"]<-"Proportion High Marsh"
mod_tab_scale[4,"fun"]<-deparse1(mod_list_surv2[[2]]$formula)
mod_tab_scale[4,"AIC"]<-mod_list_surv2[[2]]$aic

mod_tab_scale[5,"name"]<-"Nest in High Marsh"
mod_tab_scale[5,"fun"]<-deparse1(mod_list_pres2[[3]]$formula)
mod_tab_scale[5,"AIC"]<-mod_list_pres2[[3]]$aic

mod_tab_scale[6,"name"]<-"Nest in High Marsh"
mod_tab_scale[6,"fun"]<-deparse1(mod_list_surv2[[3]]$formula)
mod_tab_scale[6,"AIC"]<-mod_list_surv2[[3]]$aic


mod_tab_scale<-group_by(mod_tab_scale,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)




# 3. Spatial Structure: Does the spatial composition of the surrounding habitat matter? Edges provide taller vegetation and elevation?
# Compare entropy and correlation with proportion high marsh
mod_list_pres3<-list()
mod_list_surv3<-list()

# ALL
mod_list_pres3[[1]]<-glm(y~ent_txt+cor_txt+HIMARSH, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv3[[1]]<-glm(y~ent_txt+cor_txt+HIMARSH, 
                        data=surv_train,
                        family = binomial(link="logit"))
vif_pres<-vif(mod_list_pres3[[1]])
vif_surv<-vif(mod_list_surv3[[1]])


# Entropy -dissimilarity
mod_list_pres3[[2]]<-glm(y~ent_txt, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv3[[2]]<-glm(y~ent_txt, 
                        data=surv_train,
                        family = binomial(link="logit"))


# Correlation - similarity
mod_list_pres3[[3]]<-glm(y~cor_txt, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv3[[3]]<-glm(y~cor_txt, 
                        data=surv_train,
                        family = binomial(link="logit"))


# store model selection factors in a table
mod_tab_spatstr<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres3)),
                           class=rep("Spatial Structure",2*length(mod_list_pres3)),
                           name=rep(NA,2*length(mod_list_pres3)),
                           fun=rep(NA,2*length(mod_list_pres3)),
                           max_VIF=rep(NA,2*length(mod_list_pres3)),
                           VIF_var=rep(NA,2*length(mod_list_pres3)),
                           AIC=rep(NA,2*length(mod_list_pres3)),
                           dAIC=rep(NA,2*length(mod_list_pres3)),
                           AUC=rep(NA,2*length(mod_list_pres3)))

mod_tab_spatstr[1,"name"]<-"All"
mod_tab_spatstr[1,"fun"]<-deparse1(mod_list_pres3[[1]]$formula)
mod_tab_spatstr[1,"max_VIF"]<-max(vif_pres)
mod_tab_spatstr[1,"VIF_var"]<-names(vif_pres[which(vif_pres==max(vif_pres))])
mod_tab_spatstr[1,"AIC"]<-mod_list_pres3[[1]]$aic

mod_tab_spatstr[2,"name"]<-"All"
mod_tab_spatstr[2,"fun"]<-deparse1(mod_list_surv3[[1]]$formula)
mod_tab_spatstr[2,"max_VIF"]<-max(vif_surv)
mod_tab_spatstr[2,"VIF_var"]<-names(vif_surv[which(vif_surv==max(vif_surv))])
mod_tab_spatstr[2,"AIC"]<-mod_list_surv3[[1]]$aic

mod_tab_spatstr[3,"name"]<-"Entropy"
mod_tab_spatstr[3,"fun"]<-deparse1(mod_list_pres3[[2]]$formula)
mod_tab_spatstr[3,"AIC"]<-mod_list_pres3[[2]]$aic

mod_tab_spatstr[4,"name"]<-"Entropy"
mod_tab_spatstr[4,"fun"]<-deparse1(mod_list_surv3[[2]]$formula)
mod_tab_spatstr[4,"AIC"]<-mod_list_surv3[[2]]$aic

mod_tab_spatstr[5,"name"]<-"Correlation"
mod_tab_spatstr[5,"fun"]<-deparse1(mod_list_pres3[[3]]$formula)
mod_tab_spatstr[5,"AIC"]<-mod_list_pres3[[3]]$aic

mod_tab_spatstr[6,"name"]<-"Correlation"
mod_tab_spatstr[6,"fun"]<-deparse1(mod_list_surv3[[3]]$formula)
mod_tab_spatstr[6,"AIC"]<-mod_list_surv3[[3]]$aic


mod_tab_spatstr<-group_by(mod_tab_spatstr,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)



################################
# compare models
mod_tab<-rbind(mod_tab_habstr,mod_tab_scale,mod_tab_spatstr)

# PCA more important for nest placement, high marsh presence more important for success
# Nest in high marsh (local habitat) more important than surrounding habitat for both
# Correlation more important than entropy, but not as important as proportion high marsh

#################################


# 4. Do multiple components influence nest placement and survival?
mod_list_pres4<-list()
mod_list_surv4<-list()

# Compare top habitat structure and top scale together
mod_list_pres4[[1]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Scale" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4],"+",
       sub(".*~", "", mod_tab[mod_tab$class=="Habitat Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
       data=pres_train,
       family = binomial(link="logit"))
mod_list_surv4[[1]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Scale" & mod_tab$response=="Success" & mod_tab$dAIC==0, 4],"+",
                               sub(".*~", "", mod_tab[mod_tab$class=="Habitat Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
                        data=pres_train,
                        family = binomial(link="logit"))
vif_pres1<-vif(mod_list_pres4[[1]])
vif_surv1<-vif(mod_list_surv4[[1]])


# Compare top habitat structure and top spatial structure together
mod_list_pres4[[2]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Spatial Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4],"+",
                               sub(".*~", "", mod_tab[mod_tab$class=="Habitat Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv4[[2]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Spatial Structure" & mod_tab$response=="Success" & mod_tab$dAIC==0, 4],"+",
                               sub(".*~", "", mod_tab[mod_tab$class=="Habitat Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
                        data=pres_train,
                        family = binomial(link="logit"))
vif_pres2<-vif(mod_list_pres4[[2]])
vif_surv2<-vif(mod_list_surv4[[2]])


# Compare top Scale and top spatial structure together
mod_list_pres4[[3]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Spatial Structure" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4],"+",
                               sub(".*~", "", mod_tab[mod_tab$class=="Scale" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv4[[3]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Spatial Structure" & mod_tab$response=="Success" & mod_tab$dAIC==0, 4],"+",
                               sub(".*~", "", mod_tab[mod_tab$class=="Scale" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4]))), 
                        data=pres_train,
                        family = binomial(link="logit"))
vif_pres3<-vif(mod_list_pres4[[3]])
vif_surv3<-vif(mod_list_surv4[[3]])



# store model selection factors in a table
mod_tab_combo<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres4)),
                        class=rep("Combined Classes",2*length(mod_list_pres4)),
                        name=rep(NA,2*length(mod_list_pres4)),
                        fun=rep(NA,2*length(mod_list_pres4)),
                        max_VIF=rep(NA,2*length(mod_list_pres4)),
                        VIF_var=rep(NA,2*length(mod_list_pres4)),
                        AIC=rep(NA,2*length(mod_list_pres4)),
                        dAIC=rep(NA,2*length(mod_list_pres4)),
                        AUC=rep(NA,2*length(mod_list_pres4)))

mod_tab_combo[1,"name"]<-"Habitat Structure + Scale"
mod_tab_combo[1,"fun"]<-deparse1(mod_list_pres4[[1]]$formula)
mod_tab_combo[1,"max_VIF"]<-max(vif_pres1)
mod_tab_combo[1,"VIF_var"]<-names(vif_pres1[which(vif_pres1==max(vif_pres1))])[1]
mod_tab_combo[1,"AIC"]<-mod_list_pres4[[1]]$aic

mod_tab_combo[2,"name"]<-"Habitat Structure + Scale"
mod_tab_combo[2,"fun"]<-deparse1(mod_list_surv4[[1]]$formula)
mod_tab_combo[2,"max_VIF"]<-max(vif_surv1)
mod_tab_combo[2,"VIF_var"]<-names(vif_surv1[which(vif_surv1==max(vif_surv1))])[1]
mod_tab_combo[2,"AIC"]<-mod_list_surv4[[1]]$aic

mod_tab_combo[3,"name"]<-"Habitat + Spatial Structure"
mod_tab_combo[3,"fun"]<-deparse1(mod_list_pres4[[2]]$formula)
mod_tab_combo[3,"max_VIF"]<-max(vif_pres2)
mod_tab_combo[3,"VIF_var"]<-names(vif_pres2[which(vif_pres2==max(vif_pres2))])[1]
mod_tab_combo[3,"AIC"]<-mod_list_pres4[[2]]$aic

mod_tab_combo[4,"name"]<-"Habitat + Spatial Structure"
mod_tab_combo[4,"fun"]<-deparse1(mod_list_surv4[[2]]$formula)
mod_tab_combo[4,"max_VIF"]<-max(vif_surv2)
mod_tab_combo[4,"VIF_var"]<-names(vif_surv2[which(vif_surv2==max(vif_surv2))])[1]
mod_tab_combo[4,"AIC"]<-mod_list_surv4[[2]]$aic

mod_tab_combo[5,"name"]<-"Spatial Structure + Scale"
mod_tab_combo[5,"fun"]<-deparse1(mod_list_pres4[[3]]$formula)
mod_tab_combo[5,"max_VIF"]<-max(vif_pres3)
mod_tab_combo[5,"VIF_var"]<-names(vif_pres3[which(vif_pres3==max(vif_pres3))])[1]
mod_tab_combo[5,"AIC"]<-mod_list_pres4[[3]]$aic

mod_tab_combo[6,"name"]<-"Spatial Structure + Scale"
mod_tab_combo[6,"fun"]<-deparse1(mod_list_surv4[[3]]$formula)
mod_tab_combo[6,"max_VIF"]<-max(vif_surv3)
mod_tab_combo[6,"VIF_var"]<-names(vif_surv3[which(vif_surv3==max(vif_surv3))])[1]
mod_tab_combo[6,"AIC"]<-mod_list_surv4[[3]]$aic




mod_tab_combo<-group_by(mod_tab_combo,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)



################################
# compare models
mod_tab<-rbind(mod_tab,mod_tab_combo)

#################################



# 5. Does marsh resilience impact habitat structure and spatial structure? Can birds detect resilience? 

mod_list_pres5<-list()
mod_list_surv5<-list()

# compare best models of habitat structure with uvvr

# UVVR with Habitat Structure
mod_list_pres5[[1]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Combined Classes" & mod_tab$response=="Presence" & mod_tab$dAIC==0, 4],"+uvvr_mean")), 
                        data=pres_train,
                        family = binomial(link="logit"))

mod_list_surv5[[1]]<-glm(as.formula(paste0(mod_tab[mod_tab$class=="Combined Classes" & mod_tab$response=="Success" & mod_tab$dAIC==0, 4],"+uvvr_mean")), 
                        data=surv_train,
                        family = binomial(link="logit"))

vif_pres1<-vif(mod_list_pres5[[1]])
vif_surv1<-vif(mod_list_surv5[[1]])




# Just UVVR
mod_list_pres5[[2]]<-glm(y~uvvr_mean, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv5[[2]]<-glm(y~uvvr_mean, 
                        data=surv_train,
                        family = binomial(link="logit"))





# store model selection factors in a table
mod_tab_res<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres5)),
                        class=rep("Marsh Resilience",2*length(mod_list_pres5)),
                        name=rep(NA,2*length(mod_list_pres5)),
                        fun=rep(NA,2*length(mod_list_pres5)),
                        max_VIF=rep(NA,2*length(mod_list_pres5)),
                        VIF_var=rep(NA,2*length(mod_list_pres5)),
                        AIC=rep(NA,2*length(mod_list_pres5)),
                        dAIC=rep(NA,2*length(mod_list_pres5)),
                        AUC=rep(NA,2*length(mod_list_pres5)))

mod_tab_res[1,"name"]<-"Habitat + Resilience"
mod_tab_res[1,"fun"]<-deparse1(mod_list_pres5[[1]]$formula)
mod_tab_res[1,"max_VIF"]<-max(vif_pres1)
mod_tab_res[1,"VIF_var"]<-names(vif_pres1[which(vif_pres1==max(vif_pres1))])[1]
mod_tab_res[1,"AIC"]<-mod_list_pres5[[1]]$aic

mod_tab_res[2,"name"]<-"Habitat + Resilience"
mod_tab_res[2,"fun"]<-deparse1(mod_list_surv5[[1]]$formula)
mod_tab_res[2,"max_VIF"]<-max(vif_surv1)
mod_tab_res[2,"VIF_var"]<-names(vif_surv1[which(vif_surv1==max(vif_surv1))])[1]
mod_tab_res[2,"AIC"]<-mod_list_surv5[[1]]$aic

mod_tab_res[3,"name"]<-"Marsh Resilience"
mod_tab_res[3,"fun"]<-deparse1(mod_list_pres5[[2]]$formula)
mod_tab_res[3,"AIC"]<-mod_list_pres5[[2]]$aic

mod_tab_res[4,"name"]<-"Marsh Resilience"
mod_tab_res[4,"fun"]<-deparse1(mod_list_surv5[[2]]$formula)
mod_tab_res[4,"AIC"]<-mod_list_surv5[[2]]$aic




mod_tab_res<-group_by(mod_tab_res,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)





# 6. Do these habitat relationships change across the breeding range? (eg more success at center of range?)

mod_list_pres6<-list()
mod_list_surv6<-list()

# compare best models of habitat with latitude

# latitude with best habitat model
mod_list_pres6[[1]]<-glm(as.formula(paste0(mod_tab_res[mod_tab_res$class=="Marsh Resilience" & mod_tab_res$response=="Presence" & mod_tab_res$dAIC==0, 4],"+latitude")), 
                        data=pres_train,
                        family = binomial(link="logit"))

mod_list_surv6[[1]]<-glm(as.formula(paste0(mod_tab_res[mod_tab_res$class=="Marsh Resilience" & mod_tab_res$response=="Success" & mod_tab_res$dAIC==0, 4],"+latitude")), 
                        data=surv_train,
                        family = binomial(link="logit"))

vif_pres1<-vif(mod_list_pres6[[1]])
vif_surv1<-vif(mod_list_surv6[[1]])




# Just latitude
mod_list_pres6[[2]]<-glm(y~latitude, 
                        data=pres_train,
                        family = binomial(link="logit"))
mod_list_surv6[[2]]<-glm(y~latitude, 
                        data=surv_train,
                        family = binomial(link="logit"))





# store model selection factors in a table
mod_tab_lat<-data.frame(response=rep(c("Presence","Success"),length(mod_list_pres6)),
                        class=rep("Latitude",2*length(mod_list_pres6)),
                        name=rep(NA,2*length(mod_list_pres6)),
                        fun=rep(NA,2*length(mod_list_pres6)),
                        max_VIF=rep(NA,2*length(mod_list_pres6)),
                        VIF_var=rep(NA,2*length(mod_list_pres6)),
                        AIC=rep(NA,2*length(mod_list_pres6)),
                        dAIC=rep(NA,2*length(mod_list_pres6)),
                        AUC=rep(NA,2*length(mod_list_pres6)))

mod_tab_lat[1,"name"]<-"Habitat + Latitude"
mod_tab_lat[1,"fun"]<-deparse1(mod_list_pres6[[1]]$formula)
mod_tab_lat[1,"max_VIF"]<-max(vif_pres1)
mod_tab_lat[1,"VIF_var"]<-names(vif_pres1[which(vif_pres1==max(vif_pres1))])[1]
mod_tab_lat[1,"AIC"]<-mod_list_pres6[[1]]$aic

mod_tab_lat[2,"name"]<-"Habitat + Latitude"
mod_tab_lat[2,"fun"]<-deparse1(mod_list_surv6[[1]]$formula)
mod_tab_lat[2,"max_VIF"]<-max(vif_surv1)
mod_tab_lat[2,"VIF_var"]<-names(vif_surv1[which(vif_surv1==max(vif_surv1))])[1]
mod_tab_lat[2,"AIC"]<-mod_list_surv6[[1]]$aic

mod_tab_lat[3,"name"]<-"Latitude"
mod_tab_lat[3,"fun"]<-deparse1(mod_list_pres6[[2]]$formula)
mod_tab_lat[3,"AIC"]<-mod_list_pres6[[2]]$aic

mod_tab_lat[4,"name"]<-"Latitude"
mod_tab_lat[4,"fun"]<-deparse1(mod_list_surv6[[2]]$formula)
mod_tab_lat[4,"AIC"]<-mod_list_surv6[[2]]$aic




mod_tab_lat<-group_by(mod_tab_lat,response)%>%mutate(dAIC=AIC-min(AIC))%>%
  ungroup()%>%
  arrange(response,dAIC)




################################
# combine all tables
mod_tab<-rbind(mod_tab,mod_tab_res,mod_tab_lat)


mod_list_pres_final<-c(mod_list_pres1,mod_list_pres2,mod_list_pres3,mod_list_pres4,mod_list_pres5,mod_list_pres6)
mod_list_surv_final<-c(mod_list_surv1,mod_list_surv2,mod_list_surv3,mod_list_surv4,mod_list_surv5,mod_list_surv6)


# PCA more important for nest placement, high marsh presence more important for success
# Nest in high marsh (local habitat) more important than surrounding habitat for both
# Correlation more important than entropy, but not as important as proportion high marsh

#################################