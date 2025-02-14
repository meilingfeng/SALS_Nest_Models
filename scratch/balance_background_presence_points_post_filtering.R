#We only filtered incomplete records from the GLM (complete) dataset, so only check imbalance in the presence dataset for these models
# for GLMs data availability is comparatively more important, especially for rare species (prevalence <0.1)

#  if(sum(pres_dat$presence)/nrow(pres_dat)<0.3&sum(pres_dat$presence)/nrow(pres_dat)>0.1){

# randomly select the same number of background points as the complete nest records
# sample the same sample size as each region
#  reg_size<-as.vector(table(pres_dat[pres_dat$presence==1,]$Region))

#  temp1<-list()

# set.seed(123)
#  for(i in 1:length(reg_size)){
#    temp1[[i]]<-pres_dat[pres_dat$Region==regions[[i]]&pres_dat$presence==0,]
#    temp1[[i]]<-temp1[[i]][sample(1:nrow(temp1[[i]]), reg_size[i]), ]
#  }
# temp1[[length(reg_size)+1]]<-pres_dat[pres_dat$presence==1,]
#  pres_dat <- do.call("rbind", temp1)
#  pres_dat<-dplyr::filter(pres_dat,!is.na(id))

#  surv_dat<-pres_dat[!is.na(pres_dat$fate),] 
#  surv_dat<-surv_dat[!(is.na(surv_dat$time_since_tide)),]
#  }
