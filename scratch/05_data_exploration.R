library(tidyverse)
library(car)#vif
library(patchwork)
library(sf) #spatial data
library(PerformanceAnalytics)#correlation plots
library(gstat) #variograms
library(lme4) #mixed models
library(tseries) #autocorrelation
library(mgcv)
library(AICcmodavg)

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. parameters
# what prediction resolution?
reso<-30

# use random background points or veg plots as absences?
ab_type<-"b"
#ab_type<-"v" 

# pick which predictors and response to use for models
predictors <- "uvvr"
#predictors <- "no uvvr" #removes observations mising uvvr data


# veg class codes
veg_codes<-data.frame(veg_code=c(1:2,4:9),
                      veg_class=c("HIMARSH","LOMARSH","MUD","PHRG","POOL","STRM","TERRBRD","UPLND"))

## 3. Load observations with predictors
# list of predictor files for each zone as a raster stack
load(paste0(path_out,"predictor_files_all_zones_",reso,"m.rds"))

# point predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_local.csv"))%>%
  filter(bp%in%c("p",ab_type))%>%
  mutate(presence=ifelse(bp=="p",1,0))%>%
  left_join(veg_codes,by="veg_class")

# buffered predictors
dat<-read.csv(paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"))%>%
  dplyr::select(id,HIMARSH,LOMARSH,POOL,PHRG,STRM,MUD,TERRBRD)%>% #remove UPLND, keep phrag since directly management related
  right_join(dat,by="id")%>%
  mutate(veg_code=as.factor(veg_code),
         #create binary variables for each veg class
         HIMARSH_class=ifelse(veg_class=="HIMARSH",1,0),
         LOMARSH_class=ifelse(veg_class=="LOMARSH",1,0),
         TERRBRD_class=ifelse(veg_class=="TERRBRD",1,0),
         STRM_class=ifelse(veg_class=="STRM",1,0),
         MUD_class==ifelse(veg_class=="MUD",1,0),
         PHRG_class=ifelse(veg_class=="PHRG",1,0))


# filter nests to thinned data
filt<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_10_20_fail_thin_1m_dist_errors_removed.shp"))
dat<-filter(dat,id%in%filt$id | bp%in%c("b","v"))

## 4. Remove missing values from predictors
#data with all vegetation class and NAIP predictors
dat_comp<-dat[complete.cases(dat[ , c("veg_code","ent_txt")]),]
dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]<-dat_comp[,c("HIMARSH","LOMARSH","PHRG","STRM","MUD","TERRBRD")]%>%
  replace(is.na(.),0)
#associated predictors
form<-y~ndvi+cor_txt+veg_code+HIMARSH+LOMARSH+PHRG+MUD+STRM+TERRBRD+pca+latitude 

#data with all predictors including UVVR
dat_comp_all<-dat_comp[complete.cases(dat_comp[ , c("uvvr_mean")]),]
#associated predictors
form_all<-y~ndvi+cor_txt+veg_class+HIMARSH+LOMARSH+PHRG+MUD+STRM+TERRBRD+pca+uvvr_mean+latitude



## 5. pick which predictors to use for models
# don't use UVVR data (increases data availability)
if(predictors=="no uvvr" & analysis=="presence"){
  mod_dat<-dat_comp
  mod_form<-form
}

if(predictors=="no uvvr" & analysis=="survival"){
  mod_dat<-dat_comp[!is.na(dat_comp$fate),] #filter just the nests with known fates
  mod_form<-form
}
# use UVVR data (limits data availability by more than half)
if(predictors=="uvvr" & analysis=="presence"){
  mod_dat<-dat_comp_all
  mod_form<-form_all
}

if(predictors=="uvvr" & analysis=="survival"){
  mod_dat<-dat_comp_all[!is.na(dat_comp_all$fate),] #filter just the nests with known fates
  mod_form<-form_all
}




## 6. Check for balance among sample groups 
nests<-mod_dat[is.na(mod_dat$site) == F,]


#only 1 obs in 2010, remove
ggplot(nests,aes(x=Year))+
  geom_histogram()
mod_dat<-filter(mod_dat,is.na(Year)|Year!=2010)

#remove sites with less than 10 observations
ggplot(nests,aes(x=site))+
  geom_histogram(stat="count")

site_n<-summarise(group_by(nests,site),count=n())

mod_dat<-filter(mod_dat,is.na(site)|mod_dat$site %in% site_n[site_n$count>=10,]$site)




## 6. divide into testing and training data 80 training/20 testing
table(mod_dat$presence)

dat_split<-mod_dat %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
#check ratio of data
map_int(dat_split, nrow)
#check balance in each split dataset
table(dat_split[["train"]]$presence) #retains same balance
table(dat_split[["test"]]$presence)
#split into separate objects
train<-dat_split[["train"]]
test<-dat_split[["test"]]




### Build models

##  pick which response to use for models
if(analysis=="survival"){
  y<-train$fate
}
if(analysis=="presence"){
  y<-train$presence
}


## 8. Run full model
mod<-glm(mod_form, 
         data=train,
         family = binomial(link="logit"))




### Data Exploration (Zuur et al)
#------------------------------------------------------------

## 1. Look for outliers in X and Y
#----------------------------------------

#outliers in...
#Predictors: Transformation is an option when you dont want to lose any data
#Response: choose a statistical method that uses a probability distribution that allows greater variation for large mean values
# Gamma for continuous or poisson/negative binomial for counts

#Diagnostic tools:
#Cooks statistic for linear regression (change in regression parameters as each observation is individually omitted)
#if there are multiple outliers with similar values, they wont be detected



#look at just the predictors
vars<-c("HIMARSH","LOMARSH","PHRG","EDGE","cor_txt","ndvi","pca","ent_txt","uvvr_mean","latitude")
preds<-mod_dat[,vars]

# can cause over dispersion in poisson glms
# when outcome is binary, outliers are better handled

#Boxplot, histogram, scatterplot are graphical tools for outlier detection

for(i in 1:(length(vars))){
par(mfrow = c(1, 3),mar=c(5,4,7,2)) #bottom, left, top, right
hist(preds[,i], main = "Histogram")
boxplot(preds[,i], main = "Boxplot")
qqnorm(preds[,i], main = "Normal Q-Q plot")
mtext(paste0(vars[i]), side = 3, line = - 2, outer = TRUE)
}
#leverage
i_n <- influence(mod)$hat # calculate the influence of data points in the model
train[which.max(i_n),]

#cooks distance
# R code
c_d <- cooks.distance(mod)
train[which.max(c_d),]


#All predictors are non-normal- Veg cover classes have high zeros and maxes, uvvr and NDVI and low marsh right skewed, ent zero inflated, 





## 2. Is there collinearity among the covariates?
#-------------------------------------------------------------
#ignored collinearity leads to nothing being significant

# a) plot pairwise scatterplots, (alternatively look at correlation coefficients, or a PCA biplot of all covariates)
chart.Correlation(preds, histogram = TRUE, method = "spearman")
#entropy still correlated with PCA and NDVI, use correlation instead


# b) Look at variance inflation factor (larger if the covariate R2 is large, most variation in the covariate is explained by variation in all the other covariates)
  # remove variables with the highest VIF until all remaining variable have below a select threshold
  # even moderate collinearity can mask weak ecological signals so using a low threshold like 3 or 2 is important.

# Full model
vif(mod)      
summary(mod) 
Anova(mod,type=3)

#model with veg class
mod_form
mod_class<-glm(y~ndvi + cor_txt + veg_class + pca + uvvr_mean + latitude, 
               data=train,
               family = binomial(link="logit"))
vif(mod_class)
summary(mod_class)    
Anova(mod_class,type=3)

#model with veg proportions
mod_prop<-glm(y~ndvi + cor_txt +  HIMARSH + LOMARSH + PHRG + EDGE + pca + uvvr_mean + latitude, 
              data=train,
              family = binomial(link="logit"))
vif(mod_prop)
summary(mod_prop)  
Anova(mod_prop,type=3)


#which model looks better?
#anova compares models to see if adding variables significantly reduces noise (resid deviance), order smallest to most complex models
anova(mod_class,mod_prop,mod, test="Chisq")
anova(mod_class,mod_prop, test="Chisq") #looks like neither is significantly different

aictab(list(mod,mod_class,mod_prop),modnames = c("Full","Classes","Proportions")) # full model fits the best...



## 3. What are the relationships between Y and X?
#----------------------------------------------------------
ggplot(train%>%
         pivot_longer(c("ndvi","cor_txt","HIMARSH","LOMARSH","PHRG","EDGE","pca","uvvr_mean"),
                      names_to = "variables",values_to = "values"),
       aes(x=values,y=presence))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_wrap(~variables,scales = "free")

#PHRG and EDGE are U shaped, ndvi quadratic?




## 4. Are response observations independent?
#----------------------------------------------------------------------
# inflates type 1 error or false positives and raises p values in regression

#Good to check dependence before and after modeling
# plot response against spatial and temporal variables
# check that model residuals do not have dependence structure

#formal way to check is with autocorrelation
  #plot auto correlation funcitons (ACF) for regular time series (pearsons correlation between time series and lagged time series)
  #plot variograms for irregularly spaced time series and spatial data
    # straight diagonal line indicates perfect correlation. Want it to reach a sill where the linear trend devolves


# a) variograms https://csegrecorder.com/articles/view/the-variogram-basics-a-visual-intro-to-useful-geostatistical-concepts 

  # spatial autocorrelation
train$res<-mod$residuals  
dat2<-st_as_sf(train, coords=c("longitude","latitude"))
dat2<-arrange(dat2,Year)

gram<-variogram(res~1,data=dat2)  #filter(dat2,site=="ER")
plot(gram)

fit.gram<-fit.variogram(gram, model = vgm(psill=4,model= "Exp",range=1,nugget=1.5))
plot(gram,fit.gram)#correlation in residuals seems to diminish after 1 meter distance??



  # temporal autocorrelation
dat3<-dat2%>%group_by(site,Year)%>%summarise(res=mean(res))%>%
  arrange(site,Year)%>%
  mutate(lag1=lag(res,n=1),
         lag2=lag(res,n=2))

plot(dat3$res,dat3$lag2,xlab="t",ylab="t-1")
cor(dat3$res,dat3$lag2)


#plot residuals over time and lat
plot(train$res~train$latitude)
plot(train$res~train$Year)
stripchart(train$res~train$site,pch=1,method="jitter")

#solutions:
#model spatial or temporal relationships
#use lagged response variable as covariate
#mixed effects models
#include the variable in the model
#nest data in hierarchical structure 

# b) compare a structured model to one without structure using a hypothesis test like ANOVA (JUST FOR NEST SURVIVAL)

if(analysis=="survival"){
  
}





## 5. Should we consider interactions?
#---------------------------------------------------
# a) Does the effect of veg class depend on habitat continuity (veg proportion)?
mod.class.cont<-glm(y~ndvi + cor_txt +  HIMARSH + LOMARSH*LOMARSH_class + PHRG*PHRG_class + EDGE*EDGE_class + pca + uvvr_mean +latitude, 
                              data=train,
                              family = binomial(link="logit"))
# nest presence
  # High marsh, low marsh, mud, stream locations that are buffered by more high marsh more associated with nesting
  # Phrag and upland areas less buffered by high marsh are more associated with nesting 
hi_prop<-ggplot(mod_dat,aes(x=as.numeric(presence),y=HIMARSH,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Proportion High Marsh (15 m radius)",) + 
  theme_classic()
  # Terrestrial borders that are buffered by more low marsh more associated with nesting
brd_prop1<-ggplot(mod_dat,aes(x=as.numeric(presence),y=LOMARSH,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Proportion Low Marsh (15 m radius)",) + 
  theme_classic()
  # nesting on the terrestrial border tends to extend deeper in the terrestrial border- is this transitioning to high marsh?
brd_prop2<-ggplot(mod_dat,aes(x=as.numeric(presence),y=TERRBRD,colour=as.factor(EDGE_class),shape=as.factor(EDGE_class),linetype=as.factor(EDGE_class)))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Proportion Terrestrial Border (15 m radius)",) + 
  theme_classic()
  # nests in phragmites are in larger patches of phrag
brd_prop3<-ggplot(mod_dat,aes(x=as.numeric(presence),y=PHRG,colour=as.factor(PHRG_class),shape=as.factor(PHRG_class),linetype=as.factor(PHRG_class)))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Proportion Phragmites (15 m radius)",) + 
  theme_classic()
  # hypotheses combining the above plots:
  ## High marsh, low marsh, mud, stream locations that are buffered by more high marsh more associated with nesting (More nests at edges of high marsh)
  ## Terrestrial borders that are buffered by more low marsh more associated with nesting (More nests where low marsh meets upland)
ggsave(hi_prop,filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_proportion_interaction.png"), width = 7, height = 6, dpi = "retina")
ggsave((brd_prop1/brd_prop2/brd_prop3),filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/upland_proportion_interaction.png"), width = 8, height = 12, dpi = "retina")



#nest survival
# successful nests in phragmites are in sparser patches of phrag
phrag_prop_suc<-ggplot(mod_dat%>%filter(!(is.na(fate))),aes(x=as.numeric(fate),y=PHRG,colour=as.factor(PHRG_class),shape=as.factor(PHRG_class),linetype=as.factor(PHRG_class)))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "Proportion Phragmites (15 m radius)",) + 
  theme_classic()

# successful nests on the terrestrial border are closer to the marsh
brd_prop_suc<-ggplot(mod_dat%>%filter(!(is.na(fate))),aes(x=as.numeric(fate),y=TERRBRD,colour=as.factor(LOMARSH_class),shape=as.factor(LOMARSH_class),linetype=as.factor(LOMARSH_class)))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "Proportion Terrestrial Border (15 m radius)",) + 
  theme_classic()

# Disrupted areas in high marsh (mudflats and streams) have more successful nesting
hi_prop_suc1<-ggplot(mod_dat%>%filter(!(is.na(fate))),aes(x=as.numeric(fate),y=HIMARSH,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Failed","Fledged"))+
  labs(x = "Nest Success", y = "Proportion High Marsh (15 m radius)",) + 
  theme_classic()
# successful nest not in High marsh tend to be close to high marsh
hi_prop_suc2<-ggplot(mod_dat%>%filter(!(is.na(fate))),aes(x=as.numeric(fate),y=HIMARSH,colour=as.factor(HIMARSH_class),shape=as.factor(HIMARSH_class),linetype=as.factor(HIMARSH_class)))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE) +
  #scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Site", y = "Proportion Low Marsh (15 m radius)",) + 
  theme_classic()

ggsave((hi_prop_suc1/hi_prop_suc2),filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/himarsh_proportion_success_interaction.png"), width = 7, height = 11, dpi = "retina")
ggsave(brd_prop_suc,filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/upland_proportion_success_interaction.png"), width = 7, height = 6, dpi = "retina")
ggsave(phrag_prop_suc,filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/phrag_proportion_success_interaction.png"), width = 7, height = 6, dpi = "retina")


# b) UVVR/NDVI depend on veg class? - seems like ndvi, pca, and uvvr all do, but not correlation
mod.class.int<-glm(y~ndvi*veg_code + cor_txt + pca*veg_code + uvvr_mean*veg_code+latitude, 
                    data=train,
                    family = binomial(link="logit"))

#just look at high marsh, low marsh, and phragmites
  #nest presence
int1<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG")),aes(x=as.numeric(presence),y=ndvi,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "NDVI",) + 
  theme_classic()
  #nest survival
int2<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG"),!(is.na(fate))),aes(x=as.numeric(fate),y=ndvi,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "NDVI",) + 
  theme_classic()
  #presence
int3<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG")),aes(x=as.numeric(presence),y=uvvr_mean,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Unveg:Veg Ratio",) + 
  theme_classic()
#nest survival
int4<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG"),!(is.na(fate))),aes(x=as.numeric(fate),y=uvvr_mean,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "Unveg:Veg Ratio",) + 
  theme_classic()
  #presence
int5<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG")),aes(x=as.numeric(presence),y=pca,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "Reflectance PC",) + 
  theme_classic()
#nest survival
int6<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG"),!(is.na(fate))),aes(x=as.numeric(fate),y=pca,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "Reflectance PC",) + 
  theme_classic()
  #presence
int7<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG")),aes(x=as.numeric(presence),y=cor_txt,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Absent","Present"))+
  labs(x = "Nest Site", y = "NDVI correlation",) + 
  theme_classic()
#nest survival
int8<-ggplot(mod_dat%>%filter(veg_class%in%c("HIMARSH","LOMARSH","PHRG"),!(is.na(fate))),aes(x=as.numeric(fate),y=cor_txt,colour=veg_class,shape=veg_class,linetype=veg_class))+ 
  geom_point(position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", size = 1, se = FALSE, col = "black") +
  scale_colour_manual(values = c("green", "blue","firebrick")) +
  scale_x_continuous(breaks=c(0,1),labels = c("Fail","Fledge"))+
  labs(x = "Nest Success", y = "NDVI correlation",) + 
  theme_classic()

int_ndvi<-(int1+int2)+ plot_layout(guides = "collect")
int_uvvr<-(int3+int4)+ plot_layout(guides = "collect")
int_pca<-(int5+int6)+ plot_layout(guides = "collect")
int_cor<-(int7+int8)+ plot_layout(guides = "collect")

#ggsave((int_ndvi/int_uvvr),filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/ndvi_uvvr_veg_interaction.png"), width = 7, height = 11, dpi = "retina")
#ggsave((int_pca/int_cor),filename=paste0(path_out,"Intermediate_outputs/Data_Model_Exploration/pca_cor_veg_interaction.png"), width = 7, height = 11, dpi = "retina")

# does nest association with UVVR depend on surrounding site characteristics (human modifications, upland habitat, slope)?
if(analysis=="survival"){
mod.yr<-lmer(y~ndvi*veg_class + cor_txt*veg_class + pca*veg_class + uvvr_mean*veg_class+latitude + (1|Year), 
                           data=train,
                           family = binomial(link="logit"))
mod.site<-lmer(y~ndvi*veg_class + cor_txt*veg_class + pca*veg_class + uvvr_mean*veg_class+latitude + (1|Year), 
               data=train,
               family = binomial(link="logit"))

mod.site.yr<-lmer(y~ndvi*veg_class + cor_txt*veg_class + pca*veg_class + uvvr_mean*veg_class + latitude + Year + (1|site), 
                  data=train,
                  family = binomial(link="logit"))
}
  
  
  
## 6. Model Selection and Evaluation
#---------------------------------------------------
# also try a model without significant predictors
mod.class.reduced<-glm(y~ndvi + uvvr_mean + veg_class + latitude, 
                                data=train,
                                family = binomial(link="logit"))
mod.prop.reduced<-glm(y~ndvi + uvvr_mean + veg_class + HIMARSH + LOMARSH + PHRG + EDGE +latitude, 
                 data=train,
                 family = binomial(link="logit"))

# gather all models
mods<-list(mod,mod_prop,mod_class,mod.prop.reduced,mod.class.reduced,mod.class.cont,mod.class.int)
mod.names<-c("Full","Veg Proportions","Veg Classes","Veg Proportions Reduced","Veg Classes Reduced","Veg Class*Proportion","Veg Class Interactions")


ranks<-aictab(mods,mod.names,rank=T)
anova(mod.class.reduced,mod.prop.reduced,mod_class,mod_prop,mod,mod.class.cont,mod.class.int,test="Chisq")

# select the top ranked model based on AIC
final_mod<-mods[[as.numeric(rownames(ranks[1,]))]]


final_mod_eval <- evaluate(filter(test,bp=="p"), filter(test,bp=="b"), final_mod)
