
library(tidyverse)
library(car)#vif
library(PerformanceAnalytics)#correlation plots
library(mgcv)
library(AICcmodavg)

### Set up
# -------------------------------------------
source("05a_Set_Model_Parameters.R")

#check which parameters are being used
#using subset of data with UVVR? or All data using just NAIP veg predictors?
predictors
#using background points or random veg points as absences?
ab_type




## Look for outliers in X and Y and check distributions
#-------------------------------------------------------------------

#outliers in...
#Predictors: Transformation is an option when you dont want to lose any data
#Response: choose a statistical method that uses a probability distribution that allows greater variation for large mean values
# Gamma for continuous or poisson/negative binomial for counts

#Diagnostic tools:
#Cooks statistic for linear regression (change in regression parameters as each observation is individually omitted)
#if there are multiple outliers with similar values, they wont be detected



#list all predictor variables in the full model
full_mod<-mod_list_pres[[1]]
full_mod_surv<-mod_list_surv[[1]]
terms<-unlist(strsplit(as.character(full_mod$terms[[3]]),split=" * ", fixed=TRUE))%>%
  strsplit(split=" + ", fixed=TRUE)%>%
  unlist()
vars<-terms[-c(1,4)]


pres_preds<-pres_train[,vars]
surv_preds<-surv_train[,vars]

# can cause over dispersion in poisson glms
# when outcome is binary, outliers are better handled

#Boxplot, histogram, scatterplot are graphical tools for outlier detection

for(i in 1:(length(vars))){
  par(mfrow = c(1, 3),mar=c(5,4,7,2)) #bottom, left, top, right
  hist(pres_preds[,i], main = "Histogram")
  boxplot(pres_preds[,i], main = "Boxplot")
  qqnorm(pres_preds[,i], main = "Normal Q-Q plot")
  mtext(paste0(vars[i]), side = 3, line = - 2, outer = TRUE)
}


#leverage
#i_n <- influence(mod)$hat # calculate the influence of data points in the model
#train[which.max(i_n),]

#cooks distance
# R code
#c_d <- cooks.distance(mod)
#train[which.max(c_d),]


#NDVI - right skewed, 0 inflated, outliers on tail
#cor_ext- zero inflated, outliers on either side of zero (could be normal if not for zeros)
#Highmarsh - if the nest is in highmarsh habitat, pretty even divide, ~300 difference
#HIMARSH - proportion of high marsh in 15 meter buffer * maybe we could extend this buffer and look for upland land cover types. Is the nearest edge forest? or open?
#ent_txt, zero inflated, otherwise looks normally distibuted
#pca - seems normally distributed but with 3 modes
#latitude-  seems normal
#UVVR mean - right skewed, 0 inflated, outliers on tail



#might need to rescale PCA and latitude to be on a similar degree of magnitude as other variables
par(mfrow = c(1, 3))
hist(surv_train$pca)
hist(scale(surv_train$pca,scale=T,center = F))
hist(scale(surv_train$pca,scale=T,center = T))
#try scaling without centering




## 2. Is there collinearity among the covariates?
#-------------------------------------------------------------
#ignored collinearity leads to nothing being significant

# a) plot pairwise scatterplots, (alternatively look at correlation coefficients, or a PCA biplot of all covariates)
PerformanceAnalytics::chart.Correlation(pres_preds, histogram = TRUE, method = "spearman")
#entropy, pca, and ndvi correlated
#also binary high marsh and high marsh proportion


# b) Look at variance inflation factor (larger if the covariate R2 is large, most variation in the covariate is explained by variation in all the other covariates)
# remove variables with the highest VIF until all remaining variable have below a select threshold
# even moderate collinearity can mask weak ecological signals so using a low threshold like 3 or 2 is important.

# Full model
# add to VIFs to model selection table (and AIC)
vif_out<-vif(mod_list_pres[[1]])
model_selection[1,"response"]<-"placement"
model_selection[1,"mod_name"]<-"Full"
model_selection[1,"mod_function"]<-deparse1(mod_list_pres[[1]]$formula)
model_selection[1,"max_VIF"]<-max(vif_out)
model_selection[1,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[1,"AIC"]<-mod_list_pres[[1]]$aic

vif_out<-vif(mod_list_surv[[1]])
model_selection[2,"response"]<-"success"
model_selection[2,"mod_name"]<-"Full"
model_selection[2,"mod_function"]<-deparse1(mod_list_surv[[1]]$formula)
model_selection[2,"max_VIF"]<-max(vif_out)
model_selection[2,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[2,"AIC"]<-mod_list_surv[[1]]$aic


summary(mod_list_pres[[1]])

#in response to VIF, success model has high VIF for pca, placement model high VIF for HIMARSH
  # model with only binary high marsh
mod_list_pres[[2]]<-update(mod_list_pres[[1]], .~.-HIMARSH)
  # model without pca
mod_list_surv[[2]]<-update(mod_list_surv[[1]], .~.-pca)
    # add to VIFs to model selection table
vif_out<-vif(mod_list_pres[[2]])
model_selection[3,"response"]<-"placement"
model_selection[3,"mod_name"]<-"High Marsh Presence"
model_selection[3,"mod_function"]<-deparse1(mod_list_pres[[2]]$formula)
model_selection[3,"max_VIF"]<-max(vif_out)
model_selection[3,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[3,"AIC"]<-mod_list_pres[[2]]$aic

vif_out<-vif(mod_list_surv[[2]])
model_selection[4,"response"]<-"success"
model_selection[4,"mod_name"]<-"No PC"
model_selection[4,"mod_function"]<-deparse1(mod_list_surv[[2]]$formula)
model_selection[4,"max_VIF"]<-max(vif_out)
model_selection[4,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[4,"AIC"]<-mod_list_surv[[2]]$aic

  # model with only high marsh proportions
mod_list_pres[[3]]<-update(mod_list_pres[[1]], .~.-Highmarsh)
mod_list_surv[[3]]<-update(mod_list_surv[[2]], .~.-Highmarsh)
    # add to VIFs to model selection table
vif_out<-vif(mod_list_pres[[3]])
model_selection[5,"response"]<-"placement"
model_selection[5,"mod_name"]<-"High Marsh Proportion"
model_selection[5,"mod_function"]<-deparse1(mod_list_pres[[3]]$formula)
model_selection[5,"max_VIF"]<-max(vif_out)
model_selection[5,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[5,"AIC"]<-mod_list_pres[[3]]$aic

vif_out<-vif(mod_list_surv[[3]])
model_selection[6,"response"]<-"success"
model_selection[6,"mod_name"]<-"No PC, Only High Marsh Proportion"
model_selection[6,"mod_function"]<-deparse1(mod_list_surv[[3]]$formula)
model_selection[6,"max_VIF"]<-max(vif_out)
model_selection[6,"VIF_var"]<-names(vif_out[which(vif_out==max(vif_out))])
model_selection[6,"AIC"]<-mod_list_surv[[3]]$aic


#All VIFs below 2 now. 



#label the list of models
names(mod_list_pres)<-c("Full","High Marsh Presence","High Marsh Proportion")
names(mod_list_surv)<-c("Full","No PC","No PC, Only High Marsh Proportion")

