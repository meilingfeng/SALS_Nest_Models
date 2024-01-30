library(tidyverse)
library(sf)
library(terra)
library(exactextractr)

##--------------------------------------------------------------------------------
#Summarize prediction surfaces of nesting probability within SALS extent polygons and priority sites
##--------------------------------------------------------------------------------

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. Load nest prediction rasters
pres_list<-unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Placement"),~list.files(.,pattern = "pres_BRTpreds_30mb.tif$",full.names=T)))
surv_list<-unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Success"),~list.files(.,pattern = "surv_BRTpreds_30mb.tif$",full.names=T)))

#read as raster layers
pres<-map(pres_list,rast)
surv<-map(surv_list,rast)

## 3. Thresholds that maximize sensitivity and specificity
load(paste0(path_out,"Final_outputs/Nest_Predictions/final_BRT_mods.RDS"))

## 4. SALS priority sites layer - based on high quality high marsh- where populations are stable or increasing and state level site recommendations
priority<-st_read(paste0(dat_path,"SALS_priority_marshes/SALS_Priority_Marshes_Attributed_2022.shp"))%>%
  dplyr::select(id=MarshName,area_m=area)

## 5. SALS extent layer
sf_use_s2(FALSE)
extent<-st_read(paste0(dat_path,"SALS_Breeding_Range/SALS_breeding_19Apr23.shp"))%>%
  st_set_crs(st_crs("EPSG:4269"))%>%#ask mo what coord system this is in
  dplyr::select(id=TARGET_FID,area_m)

## 6. make sure coordinate systems align
priority<-st_transform(priority,crs(pres[[1]]))#%>%
  #mutate(area_m=st_area(priority))  # calculated area in arcpro since it was being weird in R
extent<-st_transform(extent,crs(pres[[1]]))#%>%
  #mutate(area_m=st_area(extent))

## 7. define empty lists for raster summary outputs
priority_list<-list()
extent_list<-list()
priority_final<-list()
extent_final<-list()




##Calculate Average nesting and success probability (average of cell predictions within each polygon)
#----------------------------------------------------------------------------
preds<-list(pres,surv)
mark<-c("pres","surv")

#For both nest presence and nest success...
for(k in 1:length(preds)){
pred_type<-preds[[k]]#selects whether to use presence or success predictions

# and for each regional zone (1-8) along the east coast...
for(i in 1:length(pres)) {
  #summarize the number of cells weighted by proportion of coverage within each site polygon (coverage fraction)
  priority_list[[i]]<-exact_extract(pred_type[[i]], priority, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(prop=n/sum(n,na.rm=T),
           wavg=round(sum(prop*value,na.rm=T),2),
           value=round(value,2))%>%
    filter(value!="NaN")%>%
    summarise(wavg=mean(wavg,na.rm=T),
              sd=round(sqrt(sum(((value-wavg)^2)*prop,na.rm=T)),2),
              max=max(value,na.rm = T),
              min=min(value,na.rm = T))%>%
    ungroup()
  
  extent_list[[i]]<-exact_extract(pred_type[[i]], extent, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                    summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(prop=n/sum(n,na.rm=T),
           wavg=round(sum(prop*value,na.rm=T),2),
           value=round(value,2))%>%
    filter(value!="NaN")%>%
    summarise(wavg=mean(wavg,na.rm=T),
              sd=round(sqrt(sum(((value-wavg)^2)*prop,na.rm=T)),2),
              max=max(value,na.rm = T),
              min=min(value,na.rm = T))%>%
    ungroup()
}


#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
priority_final[[k]]<-do.call("rbind",priority_list)%>%
  group_by(id)%>%
  summarize(across(c(wavg,sd),\(x) mean(x,na.rm=T)),
            min=min(min,na.rm=T),
            max=max(max,na.rm = T))%>%
  dplyr::select(id,everything())
extent_final[[k]]<-do.call("rbind",extent_list)%>%
  group_by(id)%>%
  summarize(across(c(wavg,sd),\(x) mean(x,na.rm=T)),
            min=min(min,na.rm=T),
            max=max(max,na.rm = T))%>%
  dplyr::select(id,everything())

#mark the column names as nest presence or success
colnames(priority_final[[k]])[-1] <- paste(colnames(priority_final[[k]])[-1], mark[k], sep = '_')
colnames(extent_final[[k]])[-1] <- paste(colnames(extent_final[[k]])[-1], mark[k], sep = '_')

}

#combine nest presence and success into one dataframe
priority_sum<-do.call("cbind",priority_final)%>%
  distinct(id,.keep_all=T)%>%
  dplyr::select(-id.1)
extent_sum<-do.call("cbind",extent_final)%>%
  distinct(id,.keep_all=T)%>%
  dplyr::select(-id.1)




## Calculate percent area of nest presence and success (proportion of polygon area with cells above maxSS threshold)***need to remove sites on the prediction edge getting only a few values and sum sites across zones that are split between them (Parker’s Marsh Natural Area Preserve to Scarsborough Neck,)
#--------------------------------------------------------------------------------------------
thr.brt<-c(thr.p.brt,thr.s.brt)
#For both nest presence and nest success...
for(k in 1:length(preds)){
  pred_type<-preds[[k]]#selects whether to use presence or success predictions
  thr<-thr.brt[k]
#for each regional zone (1-8) along the east coast...
for(i in 1:length(pred_type)) {
  #transform probabilities into presence/absence
  temp<-pred_type[[i]]
  temp[temp>=thr] <-1
  temp[temp<thr] <-0
  #summarize the number of cells weighted by proportion of coverage for each probability value within each marsh polygon (coverage fraction)
  priority_list[[i]]<-exact_extract(temp, priority, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each nest buffer(id)
    filter(value==1)%>%
    left_join(priority,by="id")%>%
    mutate(n=(n*900),#area of a 30x30 m cell
           pct=as.numeric(round((n/area_m)*100,2)),
           n=round(n,2))%>%
    ungroup()%>%
    dplyr::select(id,area_pred=n,area_m,pct)
  
  extent_list[[i]]<-exact_extract(temp, extent, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                    summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count fraction by value and summing the weighted values in each nest buffer(id)
    filter(value==1)%>%
    left_join(mutate(extent,id=as.double(id)),by="id")%>%
    mutate(n=(n*900),#area of a 30x30 m cell
           pct=as.numeric(round((n/area_m)*100,2)),
           n=round(n,2))%>%
    ungroup()%>%
    dplyr::select(id,area_pred=n,area_m,pct)
  
  #write binary distributions to a raster file
  writeRaster(temp,paste0(path_out,"Final_outputs/Nest_Predictions/Placement/z",i,"_",mark[k],"_BRTbinary_30m_b.tif"),overwrite=T)
}
  
  #empty region dataframes indicate no SALS nests in that region
  #combine values for nests across all regions into 1 dataframe
  priority_final[[k]]<-do.call("rbind",priority_list)%>%
    group_by(id)%>%
    summarize(across(everything(),\(x) sum(x,na.rm=T)))
  extent_final[[k]]<-do.call("rbind",extent_list)%>%
    group_by(id)%>%
    summarize(across(everything(),\(x) sum(x,na.rm=T)))
  
  #mark the column names as nest presence or success
  colnames(priority_final[[k]])[-1] <- paste(colnames(priority_final[[k]])[-1], mark[k], sep = '_')
  colnames(extent_final[[k]])[-1] <- paste(colnames(extent_final[[k]])[-1], mark[k], sep = '_')
}


#combine nest presence and success into one dataframe (different lengths this time because removed sites that did not pass threshold for presence vs survival)
priority_sum2<-full_join(priority_final[[1]],priority_final[[2]],by="id")%>%
  distinct(id,.keep_all=T)
extent_sum2<-full_join(extent_final[[1]],extent_final[[2]],by="id")%>%
  distinct(id,.keep_all=T)


#combine continuous summaries and presence summaries into one table
priority_sum3<-full_join(priority_sum,priority_sum2,by="id")
extent_sum3<-full_join(extent_sum,extent_sum2,by="id")

#write to file
write.csv(priority_sum3,paste0(path_out,"Final_outputs/Nest_Predictions/SALS_Priority_Sites_Prediction_Summary_1_1_24.csv"),row.names = F)
write.csv(extent_sum3,paste0(path_out,"Final_outputs/Nest_Predictions/SALS_Extent_Prediction_Summary_1_1_24.csv"),row.names = F)
