
library(sf)
library(tidyverse)
library(tidyterra)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package
library(exactextractr)

dat_path<-"D:/Data/"
path_out<-"D:/Outputs/"

if(!exists("nests")){
  source("03a_Format_Data_Nests_and_Predictors.R")
}


## Get Appropriate Buffer Area for Nests
# ----------------------------------------

# find minimum distance between nest (nearest neighbor distance)
dist.mat <- st_distance(nests) # Great Circle distance since in lat/lon

#calc mean distance between nests within each site

#give each point its record ID
dimnames(dist.mat)<-list(nests$id,nests$id)
t<-as.data.frame(dist.mat)

#list unique site names
sites<-unique(nests$site)

#split distance matrix into sites
t3<-t2<-list()
for(i in 1:length(sites)){
t2[[i]]<-t[startsWith(rownames(t), sites[i]), startsWith(colnames(t), sites[i])]
}

#split again into sites x year
t3.2<-t2.2<-list()

grouped_ids<-group_by(nests,site,Year)%>%
  st_drop_geometry()%>%
  dplyr::select(id,site,Year)%>%
  group_split(.keep = F)

for(i in 1:length(grouped_ids)){
  t2.2[[i]]<-t[rownames(t) %in% grouped_ids[[i]]$id, colnames(t) %in% grouped_ids[[i]]$id]
}

#remove sites with only 1 nest (no distances to other points in the marsh to calculate, so df is empty)
t3<-t2[sapply(t2, class) == "data.frame"]
t3.2<-t2.2[sapply(t2.2, class) == "data.frame"]

t4<-sapply(t3, function(x){apply(x,1,mean)})
dist_mean<-unlist(t4)
dist_mean<-data.frame(id=names(dist_mean),mean_dist=dist_mean)

t4.2<-sapply(t3.2, function(x){apply(x,1,mean)})
dist_mean2<-unlist(t4.2)
dist_mean2<-data.frame(id=names(dist_mean2),mean_dist=dist_mean2)


# Calculate nearest distance
nn.dist <- sapply(t3,function(z){
  apply(z, 1, function(x) {
  return(sort(x, partial = 2)[2])
})
})%>%
  unlist()

nn.dist2 <- sapply(t3.2,function(z){
  apply(z, 1, function(x) {
    return(sort(x, partial = 2)[2])
  })
})%>%
  unlist()

# Get index for nearest distance
nn.index <- sapply(t3, function(z){
  apply(z, 1, function(x) { 
    names(x[order(x, decreasing=F)[2]]) 
  })
})%>%
  unlist()

nn.index2 <- sapply(t3.2, function(z){
  apply(z, 1, function(x) { 
    names(x[order(x, decreasing=F)[2]])
  })
})%>%
  unlist()

dist_min<-data.frame(id=names(nn.dist),min_dist=nn.dist, neighbor=nn.index)
dist_min2<-data.frame(id=names(nn.dist2),min_dist=nn.dist2, neighbor=nn.index2)

nest_dists<-left_join(st_drop_geometry(nests),dist_mean,by="id")%>%
  left_join(dist_min,by="id")%>%
  filter(mean_dist<2000)

nest_dists2<-left_join(st_drop_geometry(nests),dist_mean2,by="id")%>%
  left_join(dist_min2,by="id")%>%
  filter(mean_dist<2000)

#Nests with mean distances from other nests greater than 2000 m look like typos in arcpro

par(mfcol=c(2,2))
hist(nest_dists2$mean_dist,xlab="Mean Nest Distance within Sites*Years")
abline(v=mean(nest_dists2$mean_dist),col="red") #mean of the mean distance btween nests at a site within a season is 234 m

hist(nest_dists2$min_dist,xlab="Closest Nest Distance within Sites*Years")
abline(v=mean(nest_dists2$min_dist),col="red") # mean closest distance between nests at a site within a season is 34 m

hist(nest_dists$mean_dist,xlab="Mean Nest Distance within Sites")
abline(v=mean(nest_dists$mean_dist),col="red") # 258 m per site across years

hist(nest_dists$min_dist,xlab="Closest Nest Distance within Sites")
abline(v=mean(nest_dists$min_dist),col="red") # 16 m per site across years

# would go with min distance over mean since mean is influenced by the size of the marsh patch/site... 30 m resolution?

#look at mean across sites/years
nest_dists_sum<-nest_dists%>%group_by(site)%>%
  summarize(mean_dist=mean(mean_dist),
            min_dist=mean(min_dist))
nest_dists_sum2<-nest_dists2%>%group_by(site,Year)%>%
  summarize(mean_dist=mean(mean_dist),
            min_dist=mean(min_dist),
            lat=mean(latitude))

ggplot(nest_dists_sum2,aes(x=Year,y=min_dist,fill=lat))+
  geom_point()

ggplot(nest_dists_sum,aes(x=mean_dist))+
  geom_histogram()

## Summarize environmental data in each buffer zone (15m around each nest = 30m to match coarsest resolution of environmental data -UVVR)
# -----------------------------------------------------------------------------------------------------------------

## a. Marsh Vegetation Classes
#-------------------------------
#set coord system of nest points to raster layer
nests_buff<-st_transform(nests_buff,crs(vg_cls[[1]]))
#for each regional zone (1-8) along the east coast...
for(i in 1:length(vg_cls)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(vg_cls[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell count to area by multiplying count by cell area. 
    #then calculate the proportional area of each veg class
    group_by(id)%>%
    distinct(.keep_all = T)%>%
    mutate(area=n*area,
           prop_area=round(area/sum(area),digits=2),
           region=zones[i])%>%
    ungroup()%>%
    #join additional attributes
    left_join(nests%>%st_drop_geometry(),by='id')%>%
    left_join(veg_class,by="value")%>%
    dplyr::select(-c("value","area","n"))%>%
    distinct(.keep_all=T)%>%
    pivot_wider(names_from= "veg_class",values_from = "prop_area")
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
veg_prop<-do.call("rbind",out_list)
#rename the "NA" vegetation variable column (last column) to "MISSING"
colnames(veg_prop)[ncol(veg_prop)]<-"MISSING"





## b. Mean UVVR 
#-------------------------------
#set coord system of nest points to raster layer
nests_buff<-st_transform(nests_buff,crs(uvvr))
uvvr_mean<-exact_extract(uvvr, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                         summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_mean<-group_by(uvvr_mean,id)%>%
  mutate(weighted=n*value)%>%
  summarise(uvvr_mean=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()


## c. UVVR change 
#-------------------------------
uvvr_diff2<-exact_extract(uvvr_diff, nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                          summarize_df=T,include_cols='id')
#convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
uvvr_diff2<-group_by(uvvr_diff2,id)%>%
  mutate(weighted=n*value)%>%
  summarise(uvvr_diff=round(sum(weighted,na.rm=T),digits=5))%>%
  ungroup()



## d. NDVI
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(ndvi[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(ndvi=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
ndvi_buff<-do.call("rbind",out_list)


## e. PCA
#-------------------------------

#for each regional zone (1-8) along the east coast...
for(i in 1:length(pca)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(pca[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(pca=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
pca_buff<-do.call("rbind",out_list)




## f. Homogeneity
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_homo[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(hom_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
homo_buff<-do.call("rbind",out_list)



## g. Entropy
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_entro[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(ent_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
entro_buff<-do.call("rbind",out_list)



## h. Correlation
#-------------------------------
#for each regional zone (1-8) along the east coast...
for(i in 1:length(ndvi)) {
  #summarize the number of cells weighted by proportion of coverage within each nest buffer (coverage fraction)
  out_list[[i]]<-exact_extract(txt_corr[[i]], nests_buff, function(df) summarize(group_by(df,value,id),n=sum(coverage_fraction),.groups='drop'),
                               summarize_df=T,include_cols='id')%>%
    #convert cell coverage to weighted average by multiplying count by value and summing the weighted values in each nest buffer(id)
    group_by(id)%>%
    mutate(weighted=n*value)%>%
    summarise(cor_txt=round(sum(weighted,na.rm=T),digits=5))%>%
    ungroup()
}

#empty region dataframes indicate no SALS nests in that region
#combine values for nests across all regions into 1 dataframe
corr_buff<-do.call("rbind",out_list)




## i. join UVVR and Proportion of vegetation classes in each nest buffer into 1 table
#-------------------------------
final_dat<-left_join(veg_prop,uvvr_mean, by='id')%>%
  left_join(uvvr_diff2,by='id')%>%
  left_join(ndvi_buff,by='id')%>%
  left_join(pca_buff,by='id')%>%
  left_join(homo_buff,by='id')%>%
  left_join(entro_buff,by='id')%>%
  left_join(corr_buff,by='id')%>%
  dplyr::select(id,region,site,Year,fate, 
                HIMARSH,LOMARSH,POOL,PHRG,MISSING,STRM,MUD,UPLND,TERRBRD,
                uvvr_mean,uvvr_diff,
                ndvi,pca,hom_txt,ent_txt,cor_txt)%>%
  distinct(.keep_all = T)


write.csv(final_dat,paste0(dat_path,"SALS_nest_vars_buff15.csv"),row.names = F)

