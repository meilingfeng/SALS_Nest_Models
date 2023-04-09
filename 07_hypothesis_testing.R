library(tidyverse)
library(sf)
library(car)
library(viridis)
library("factoextra")
library(patchwork)
#https://rspatial.org/sdm/1_sdm_introduction.html

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

# point predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_local.csv"))%>%
  filter(bp %in% c("p","v"))%>%
  mutate(presence=ifelse(bp=="p",1,0))

# buffered predictors
#presence
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"))%>%
  dplyr::select(id,HIMARSH,LOMARSH,POOL,PHRG,STRM,MUD,TERRBRD)%>% #remove UPLND, keep phrag since directly management related
  right_join(dat,by="id")%>%
  mutate(Highmarsh=as.factor(ifelse(veg_class=="HIMARSH",1,0)))
#survival
dat_s<-dat%>% 
  filter(!(is.na(fate)))%>%
  mutate(fate=as.factor(fate))

# data from modeling
source("05d_Model_Selection_Evaluation.R")

#is there a difference in means between nest/non-nest or survive/non-survive?

# Compare the means of two groups
#----------------------------------------------------------------------------------------
#Format data for plotting
dat_s<- surv_train%>%
  pivot_longer(c("ent_txt","uvvr_mean","ndvi","pca","HIMARSH","ent_txt"),names_to = "variable",values_to = "values")
dat_p<- pres_train%>%
  pivot_longer(c("ent_txt","uvvr_mean","ndvi","pca","HIMARSH","ent_txt"),names_to = "variable",values_to = "values")

# 1. Plot differences
s<-ggplot(dat_s, aes(x = as.factor(fate), y = values)) + 
  geom_jitter(position = position_jitter(0.05),
              alpha=0.1,
              aes(color=as.factor(fate))) +
  geom_boxplot(
    alpha=0,
    width = .12, 
    ## remove outliers
    outlier.color = NA ## `outlier.shape = NA` works as well
  ) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = 1, 
    ## adjust height
    width = 0.7, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0, 
    point_colour="red",
    aes(fill=as.factor(fate))
  ) +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA))+
  theme_bw()+
  scale_fill_viridis(discrete=T,option="B",end=0.7)+#"magma" (or "A")"inferno" (or "B")"plasma" (or "C")"viridis" (or "D")"cividis" (or "E")"rocket" (or "F")"mako" (or "G")"turbo" (or "H")
  scale_color_viridis(discrete=T,option="B",end=0.7)+
  labs(x="Nest Fate",y="Variable")+
  facet_wrap(~variable,scales = "free",nrow = 1)+
  scale_x_discrete(labels=c("Failed","Fledged"))+
  theme(legend.position = "none")


#visually the only thing that looks different is higher entropy, pca, ndvi in successful nests



p<-ggplot(dat_p, aes(x = as.factor(presence), y = values)) + 
  geom_jitter(position = position_jitter(0.05),
              alpha=0.1,
              aes(color=as.factor(presence))) +
  geom_boxplot(
    alpha=0,
    width = .12, 
    ## remove outliers
    outlier.color = NA ## `outlier.shape = NA` works as well
  ) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = 1, 
    ## adjust height
    width = 0.7, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0, 
    point_colour="red",
    aes(fill=as.factor(presence))
  ) +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA))+
  theme_bw()+
  scale_fill_viridis(discrete=T,option="G",end=0.7)+#"magma" (or "A")"inferno" (or "B")"plasma" (or "C")"viridis" (or "D")"cividis" (or "E")"rocket" (or "F")"mako" (or "G")"turbo" (or "H")
  scale_color_viridis(discrete=T,option="G",end=0.7)+
  labs(x="Nest Site Present?",y="Variable")+
  facet_wrap(~variable,scales = "free",nrow=1)+
  scale_x_discrete(labels=c("Absent","Present"))+
  theme(legend.position = "none")



s/p
ggsave(paste0(path_out,"Final_outputs/habitat_comparison_surv_pres_b.png"),
       width=12,height=7,dpi=300,units = "in")





# 3. two-sample t-test -absent minus present difference
# The two-sample t-test is used to compare the means of two groups.-absent minus present difference
# a function takes the response ~ predictor, the categorical variable defining groups is the predictor
# need to specify that you assume the variances of the 2 groups are equal - this tells it to do a classic 2 sample t test (can lead to type 1 error)

#remove categorical predictors
terms_surv_adj<-terms_surv[-3]
terms_pres_adj<-terms_pres[-6]

#test output tables
output_p<-tibble(response = rep("Presence",length(terms_pres_adj)),
       variable = terms_pres_adj,
       n = rep(nrow(pres_train),length(terms_pres_adj)),
       Levenes_p=rep(NA,length(response)),
       mean_dif = rep(NA,length(response)),
       mean_dif_CI = rep(NA_character_,length(response)),
       t_p = rep(NA,length(response))
       )
output_s<-tibble(response = rep("Success",length(terms_surv_adj)),
                 variable = terms_surv_adj,
                 n = rep(nrow(surv_train),length(terms_surv_adj)),
                 Levenes_p=rep(NA,length(response)),
                 mean_dif = rep(NA,length(response)),
                 mean_dif_CI = rep(NA_character_,length(response)),
                 t_p = rep(NA,length(response))
)


#presence model variable t-tests
  for(j in 1:length(terms_pres_adj)){
pres_train$x<- unlist(as.vector(dplyr::select(pres_train,contains(terms_pres_adj[j]))))
# Compare the variances between groups-Levene's test (car)
output_p[j,"Levenes_p"]<-round(leveneTest(data = pres_train, x ~ as.factor(y), center = mean)$`Pr(>F)`[1],3)# sig p value (Pr>F) indicates you reject the null hypothesis that the groups have the same variances 

if(output_p[j,"Levenes_p"]>0.05){
test<-t.test(x ~ y, data = pres_train, var.equal = TRUE)# gives the CI about the mean difference between groups (fail, fledge) as well as the estimated mean in each group
}else{
test<-t.test(x ~ y, data = pres_train, var.equal = FALSE) #welch's test
}

output_p[j,"mean_dif"]<-round(test$estimate[1]-test$estimate[2],3)
output_p[j,"mean_dif_CI"]<-paste0(round(test$conf.int[1],3)," : ",round(test$conf.int[2],3))
output_p[j,"t_p"]<-round(test$p.value,3)
  }

#success model variable t-tests
for(j in 1:length(terms_surv_adj)){
  surv_train$x<- unlist(as.vector(dplyr::select(surv_train,contains(terms_surv_adj[j]))))
  # Compare the variances between groups-Levene's test (car)
  output_s[j,"Levenes_p"]<-round(leveneTest(data = surv_train, x ~ as.factor(y), center = mean)$`Pr(>F)`[1],3)# sig p value (Pr>F) indicates you reject the null hypothesis that the groups have the same variances 
  
  if(output_p[j,"Levenes_p"]>0.05){
    test<-t.test(x ~ y, data = surv_train, var.equal = TRUE)# gives the CI about the mean difference between groups (fail, fledge) as well as the estimated mean in each group
  }else{
    test<-t.test(x ~ y, data = surv_train, var.equal = FALSE) #welch's test
  }
  
  output_s[j,"mean_dif"]<-round(test$estimate[1]-test$estimate[2],3)
  output_s[j,"mean_dif_CI"]<-paste0(round(test$conf.int[1],3)," : ",round(test$conf.int[2],3))
  output_s[j,"t_p"]<-round(test$p.value,3)
}


write.csv(rbind(output_p,output_s),paste0(path_out,"Final_outputs/t_test_results",ab_type,".csv"),row.names = F)




# is there a difference in means between high probability nest vs high probability survive?
final_mod_pres<-glm(presence~ndvi+cor_txt+Highmarsh+HIMARSH+ent_txt+pca+uvvr_mean, 
                    data=dat,
                    family = binomial(link="logit"))

final_mod_surv<-glm(fate~ndvi+cor_txt+Highmarsh+HIMARSH+ent_txt+pca+uvvr_mean, 
    data=dat_s,
    family = binomial(link="logit"))

dat$pred<-final_mod_pres$fitted.values






## PCA of remote sensing layers
#select variables and scale
scale_s<-apply(surv_train%>%dplyr::select(cor_txt,ent_txt,precip,pca,ndvi,uvvr_mean),2,scale,center=T,scale=T)
scale_p<-apply(pres_train%>%dplyr::select(cor_txt,ent_txt,precip,pca,ndvi,uvvr_mean),2,scale,center=T,scale=T)



#1. get pca object
pca_s<-princomp(scale_s)
pca_p<-princomp(scale_p)

#2. look at eigen values via scree plot to choose a k rank for the data 
#(number of important variables, PCs that explain 80% variation)
#gives cumulative importance of each PC
summary(pca_p)
summary(pca_s) 
# just plots variances explained of PCs (skree plot, eyeball important PCs)
plot(pca_p) 
plot(pca_s)   

fviz_eig(pca_s, geom = "bar", bar_width = 0.3) + ggtitle("")


#3. interpret PCs by looking at loadings of old variables projected on to new
# from step 2, looks like we want the first 4 dimensions
cor_p<-fviz_pca_var(pca_p, col.circle = "black") + ggtitle("")
cor_s<-fviz_pca_var(pca_s, col.circle = "black") + ggtitle("")
cor_p+cor_s

#4. look at data plotted onto the PCs
p<-fviz_pca_biplot(pca_p, label = "var", 
                habillage = pres_train[,24] #choose what variable to label the observations by
) +
  ggtitle("")

s<-fviz_pca_biplot(pca_s, label = "var", 
                habillage = surv_train[,11] #choose what variable to label the observations by
) +
  ggtitle("")
p+s
ggsave(paste0(path_out,"Final_outputs/habitat_PCA_surv_pres_b.png"),
       width=10,height=6,dpi=300,units = "in")



#5. can relate the PCs to an outcome variable like probability of nesting

g<-glm(fate~uvvr_mean+pca+cor_txt+ent_txt+ndvi,data = dat, family = binomial(link = "logit"))
dat2$preds<-g$fitted.values

