
library(AICcmodavg)
library(tidyverse)
library(terra)
library(sf)
library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(usdm)
###########################################################################################
# look for associations between:
# 1. vegetation characteristics and nesting probabilities (nest placement and nest success)
#    Do females select for certain veg characteristics associated with nest success?
# 2. nest characteristics and nest success probability
#    Do the choices females make after selecting a nesting site aid in nest success?
###########################################################################################

# if we don't find evidence for these, are extreme flooding events driving success and masking any historical habitat drivers?
# do we see a trend in the quality of nesting habitat over time? Changes in nesting niches over time? 
  # is a change in niche due to adaptation, still within the range limits of the species, or a result of available habitat?
  # How much salt marsh/high marsh have we lost over this period? Are sparrows losing suitable habitat?



### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"C:/Users/10788/Desktop/SaltMarsh/Data/"
path_out<-"C:/Users/10788/Desktop/SaltMarsh/Outputs/"



### Surrounding Vegetation and Nest Predictions
#--------------------------------------------------------
# What habitat variables influence where females choose nest locations?

# Run an analysis with all data points, looking at surrounding veg characteristics

# look at fine resolution characteristics (Demo veg data-1m) and larger resolution characteristics (survey data-50m)
# expect a stronger relationship with larger resolution data since predictions seem to distinguish larger, site level differences

## 1. Data cleaning to transfer data into the appropriate classes for following steps

#Simple function to transfer LongLat to UTM with setted zone 
LongLatToUTM<-function(x,y,zone){
  require(sp)
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

#transfer rapid data LongLat to UTM: 
rapid_dat10=rapid_dat9
utm<-LongLatToUTM(rapid_dat10$Long,rapid_dat10$Lat,"18")
rapid_dat10$Easting<-utm$X
rapid_dat10$Northing<-utm$Y

#extract nesting probability based on rapid vegetation data location
out_list1<-list()
for(i in 1:length(pres)) {
    out_list1[[i]]<-terra::extract(pres[[i]], cbind(utm$X,utm$Y))
}
#combined the extracted presence data from eight separate layer into one
pres2<-out_list1[[1]]
for(i in 1:length(utm$X)) {
  pres2[i,]=0
  for (j in 1:length(out_list1)) {
    if(is.na(out_list1[[j]][i,])==FALSE)
      pres2[i,]<-out_list1[[j]][i,]
  }
}
rapid_dat10$presence<-pres2$lyr1

#extract nesting survival based on rapid vegetation data location
out_list2<-list()
for(i in 1:length(surv)) {
  out_list2[[i]]<-terra::extract(surv[[i]], cbind(utm$X,utm$Y))
}

#combined the extracted presence data from eight separate layer into one
surv2<-out_list2[[1]]
for(i in 1:length(utm$X)) {
  surv2[i,]=0
  for (j in 1:length(out_list2)) {
    if(is.na(out_list2[[j]][i,])==FALSE)
      surv2[i,]<-out_list2[[j]][i,]
  }
}
rapid_dat10$survival<-surv2$lyr1

## 2. Correlation plots, visualization of the importance of each veg species

#test plot
ggplot(rapid_dat10,aes(x=presence,y=alt_pres))+
  geom_point()+
  geom_smooth(method="gam")
#test plot
ggplot(rapid_dat10,aes(y=presence))+
  geom_boxplot(aes(fill=dom_species))

#store presence/absence data for all the groups
rapid_dat10.pres<-rapid_dat10[,c(1,43:55,60:61)]%>%
  melt(id=c("id","presence","survival"))%>%
  as.data.frame()
#store precentage data for all the groups
rapid_dat10.pct<-rapid_dat10[,c(1,19:31,60:61)]%>%
  melt(id=c("id","presence","survival"),na.rm = TRUE)%>%
  as.data.frame()
rapid_dat10.pct.clean<-rapid_dat10.pct[-which(rapid_dat10.pct$value==0),]

#Wilcoxon plots
ggplot(rapid_dat10.pres)+
  geom_boxplot(aes(x=as.factor(value),y=presence,fill=as.factor(value)))+
  facet_wrap(~variable)+
  stat_compare_means(aes(x=as.factor(value),y=presence,fill=as.factor(value)),
                      label.x = 1,label.y=0.5,paired=FALSE)
ggplot(rapid_dat10.pres)+
  geom_boxplot(aes(x=as.factor(value),y=survival,fill=as.factor(value)))+
  facet_wrap(~variable)+
  stat_compare_means(aes(x=as.factor(value),y=survival,fill=as.factor(value)),
                     label.x = 1,label.y=0.5,paired=FALSE)
#gam plots
ggplot(rapid_dat10.pct.clean,aes(x=value,y=presence))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_wrap(~variable,scales="free")
ggplot(rapid_dat10.pct.clean,aes(x=value,y=survival))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_wrap(~variable,scales="free")
                          
# 3. Run a principal components analysis on the percent cover variables                             
rapid_dat.pca<-rapid_dat10[,c(19:31,60:61)]
rapid_dat.nor<-scale(rapid_dat.pca)                     
rapid_dat.corr<-cor(rapid_dat.nor) 
ggcorrplot(rapid_dat.corr)                         
rapid_dat.princ<-princomp(na.omit(rapid_dat.corr))                             
summary(rapid_dat.princ)

### Nest characteristics and Nest Predictions
#--------------------------------------------------------
# What variables might be important for nest survival after the female chooses a nest location?
# Run an analysis with just nest points, looking at nest construction variables

# 1. Correlation of continuous variables with predicted nest probabilites




# 2. scatter plots (continuous variables)/boxplots (categorical variables)
#ggplot(dat_all%>%,aes(x=presence,y=Geom_histogram() + facet_wrap() 

#ggplot(dat_all%>%pivot_longer(cols = c(tll_mode,dom_pc,alt,patens,juncus,distichlis), names_to = "var", values_to = ""),
#       aes(x=presence, fill =)+Geom_density(alpha=.3) 
       
#       Aes(x=var1,y=var2)+geom_boxplot() or/+ geom_jitter()-jitter just adds the points 

       
# 3. Run a principal components analysis on the percent cover variables? See if we can condense some of the species into fewer variables.


# 4. regression for predicted probabilities (also use binary pres/abs predictions with ANOVA?)
      # choose a global model with all variables
       #choose subsets of the global model
      #compare sets for the 2 resolutions of veg data (survey and demo) and the nest data (but focus on just doing survey veg data for now)
