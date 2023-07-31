#####################################################################################################################################
## Summarize environmental data in each buffer zone - 5,10,15m - for finer resolution data
######################################################################################################################################


## Set file path to data and outputs
# -------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
#dat_path<-"/home/FCAM/mlfeng/Data/"
path_out<-"D:/Nest_Models/Outputs/"


# Load all nest locations and environmental data
if(!exists("nests_buff")){
  source("04a_Format_Load_Data_Nests_and_Predictors.R")
}

#list for the 8 region outputs
out_list<-list()


for(j in 1:length(nests_buff)){
  
  #list for the 3 buffer areas
  veg_prop<-list()
  
  ## Summarize predictors in buffers
  #---------------------------------------------------
  ## a. Marsh Vegetation Classes
  #-------------------------------
  area<-round(res(vg_cls[[1]])[1]^2)
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(vg_cls)) {
    #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
    out_list[[i]]<-exact_extract(vg_cls[[i]], st_transform(nests_buff[[j]],crs(vg_cls[[1]])), #set coord system of nest points to raster layer
                                 function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='id')%>%
      #convert cell count to area by multiplying count by cell area. 
      #then calculate the proportional area of each veg class
      group_by(id)%>%
      mutate(area=n*area,
             prop_area=round(area/sum(area),digits=2),
             region=zones[i])%>%
      ungroup()%>%
      #join additional attributes
      left_join(nests%>%st_drop_geometry(),by='id')%>%
      left_join(veg_class,by="value")%>%
      dplyr::select(-c("value","area","n"))%>%
      distinct(id,prop_area,veg_class,.keep_all=T)%>%
      pivot_wider(names_from= "veg_class",values_from = "prop_area")%>%
      unnest() #removes list cols
  }
  
  #empty region dataframes indicate no SALS nests in that region
  # align columns across dataframes 
  nms <- veg_class$veg_class   # Vector of columns you want in the data.frames (the names of each vegetation class)
  nms<-c(nms,"NA") #NA indicates outside of marsh area (no vegetation class within nest buffer)
  
  for(i in 1:length(out_list)){
    
    Missing <- setdiff(nms, names(out_list[[i]]))  # Find names of missing columns by comparing dataframe columns to vector of desired column names
    out_list[[i]][Missing] <- 0                    # Add missing columns, fill observations with '0's (0% cover)
    
  }
  
  #combine vegetation values for nests across all regions into 1 dataframe
  veg_prop[[j]]<-do.call("rbind",out_list)%>%
    distinct(id,.keep_all=T)
  #rename the "NA" vegetation variable column (last column) to "MISSING"
  colnames(veg_prop[[j]])[(ncol(veg_prop[[j]]))]<-"MISSING"
  
  
  
  
  ## b. NDVI
  #-------------------------------
  ndvi_buff<-list()
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(ndvi)) {
    #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
    out_list[[i]]<-exact_extract(ndvi[[i]], st_transform(nests_buff[[j]],crs(ndvi[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='id')%>%
      #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
      group_by(id)%>%
      mutate(weighted=n*value)%>%
      summarise(ndvi=round(sum(weighted,na.rm=T),digits=5))%>%
      ungroup()
  }
  
  #empty region dataframes indicate no SALS nests in that region
  #combine values for nests across all regions into 1 dataframe
  ndvi_buff[[j]]<-do.call("rbind",out_list)%>%
    distinct(id,.keep_all=T)
  
  
  
  
  ## c. PCA
  #-------------------------------
  pca_buff<-list()
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(pca)) {
    #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
    out_list[[i]]<-exact_extract(pca[[i]], st_transform(nests_buff[[j]],crs(pca[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='id')%>%
      #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
      group_by(id)%>%
      mutate(weighted=n*value)%>%
      summarise(pca=round(sum(weighted,na.rm=T),digits=5))%>%
      ungroup()
  }
  
  #empty region dataframes indicate no SALS nests in that region
  #combine values for nests across all regions into 1 dataframe
  pca_buff[[j]]<-do.call("rbind",out_list)%>%
    distinct(id,.keep_all=T)
  
  
  
  
  ## d. Entropy
  #-------------------------------
  entro_buff<-list()
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(ndvi)) {
    #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
    out_list[[i]]<-exact_extract(txt_entro[[i]], st_transform(nests_buff[[j]],crs(ndvi[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='id')%>%
      #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
      group_by(id)%>%
      mutate(weighted=n*value)%>%
      summarise(ent_txt=round(sum(weighted,na.rm=T),digits=5))%>%
      ungroup()
  }
  
  #empty region dataframes indicate no SALS nests in that region
  #combine values for nests across all regions into 1 dataframe
  entro_buff[[j]]<-do.call("rbind",out_list)%>%
    distinct(id,.keep_all=T)
  
  
  
  ## e. Correlation
  #-------------------------------
  corr_buff<-list()
  #for each regional zone (1-8) along the east coast...
  for(i in 1:length(ndvi)) {
    #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
    out_list[[i]]<-exact_extract(txt_corr[[i]], st_transform(nests_buff[[j]],crs(ndvi[[1]])), function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                                 summarize_df=T,include_cols='id')%>%
      #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
      group_by(id)%>%
      mutate(weighted=n*value)%>%
      summarise(cor_txt=round(sum(weighted,na.rm=T),digits=5))%>%
      ungroup()
  }
  
  #empty region dataframes indicate no SALS nests in that region
  #combine values for nests across all regions into 1 dataframe
  corr_buff[[j]]<-do.call("rbind",out_list)%>%
    distinct(id,.keep_all=T)
  
  
  
  
  ## i. join variables at each nest buffer into 1 table
  #-------------------------------
  final_dat[[j]]<-left_join(veg_prop[[j]],ndvi_buff[[j]], by='id')%>%
    left_join(pca_buff[[j]],by='id')%>%
    left_join(entro_buff[[j]],by='id')%>%
    left_join(corr_buff[[j]],by='id')%>%
    dplyr::select(id,bp,region,site,Year,fate, 
                  HIMARSH,LOMARSH,POOL,PHRG,MISSING,STRM,MUD,UPLND,TERRBRD,
                  ndvi,pca,ent_txt,cor_txt)%>%
    distinct(id,.keep_all=T)
  names(final_dat[[j]])[-c(1:6)] <- paste0(names(final_dat[[j]])[-c(1:6)], "_",buffs[j])
  
}

final_dat2<-left_join(final_dat[[1]],final_dat[[2]][-c(1:6)], by = c("id","bp","region","site","Year","fate"))%>%
  left_join(final_dat[[3]][-c(1:6)], by = c("id","bp","region","site","Year","fate"))

#save(uvvr_mean,uvvr_diff2,ndvi_buff,pca_buff,entro_buff,corr_buff, file=paste0(path_out,"Final_outputs/4c_final_files.rds"))
write.csv(final_dat,paste0(path_out,"Final_outputs/SALS_nest_vars_buff15.csv"),row.names = F)

