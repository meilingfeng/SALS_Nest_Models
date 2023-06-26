library(sf)
library(tidyverse)
library(rgdal) #package for geospatial analyses
library(terra)#updated version of raster package

###################################################################################################
# Check the minimum distance between nests to detect any remaining coordinate issues 
    # (nests that are super far away from any other nests are likely typos)
# use the minimum distance between nests as the resolution of the final prediction surface
###################################################################################################



## 1. Load data focal species nest shapefile
#----------------------------------------------------------------
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"


nests<-st_read(paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020.shp"))%>%
  rename(id=name)



## 2. Calculate minimum distance between nests within each site
# ----------------------------------------------------

# find minimum distances between nests (nearest neighbor distance)
dist.mat <- st_distance(nests) # defaults to using Great Circle distance calculation since data is in lat/lon


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

#split again into sites by year
t3.2<-t2.2<-list()

grouped_ids<-group_by(nests,site,Year)%>%
  st_drop_geometry()%>%
  dplyr::select(id,site,Year)%>%
  group_split(.keep = F)

for(i in 1:length(grouped_ids)){
  t2.2[[i]]<-t[rownames(t) %in% grouped_ids[[i]]$id, colnames(t) %in% grouped_ids[[i]]$id]
}

#remove sites with only 1 nest (no distances to other points in the site to calculate, so df is empty)
t3<-t2[sapply(t2, class) == "data.frame"]
t3.2<-t2.2[sapply(t2.2, class) == "data.frame"]


#calculate mean distance to other nests
t4<-sapply(t3, function(x){apply(x,1,mean)})
dist_mean<-unlist(t4)
dist_mean<-data.frame(id=names(dist_mean),mean_dist=dist_mean)

t4.2<-sapply(t3.2, function(x){apply(x,1,mean)})
dist_mean2<-unlist(t4.2)
dist_mean2<-data.frame(id=names(dist_mean2),mean_dist=dist_mean2)


# Calculate nearest distance (minimum distance to another nest)
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


#join mean and min distance into one table for each nest point
#Remove nests with a mean distance from other nests greater than 1000 m, likely a typo
nest_dists<-left_join(st_drop_geometry(nests),dist_mean,by="id")%>%
  left_join(dist_min,by="id")%>%
  filter(min_dist<1000)

nest_dists2<-left_join(st_drop_geometry(nests),dist_mean2,by="id")%>%
  left_join(dist_min2,by="id")%>%
  filter(min_dist<1000)


# some plots of nest distances summarized within sites and years
par(mfcol=c(2,2))
hist(nest_dists2$mean_dist,xlab="Mean Nest Distance within Sites*Years")
abline(v=mean(nest_dists2$mean_dist),col="red") #mean of the mean distance btween nests at a site within a season is 234 m

hist(nest_dists2$min_dist,xlab="Closest Nest Distance within Sites*Years")
abline(v=mean(nest_dists2$min_dist),col="red") # mean closest distance between nests at a site within a season is 34 m

hist(nest_dists$mean_dist,xlab="Mean Nest Distance within Sites")
abline(v=mean(nest_dists$mean_dist),col="red") # 258 m per site across years

hist(nest_dists$min_dist,xlab="Closest Nest Distance within Sites")
abline(v=mean(nest_dists$min_dist),col="red") # 16 m per site across years


# Use the minimum distance per site instead of mean distance since mean is influenced by the size of the marsh patch/site

#look at mean across sites/years
nest_dists_sum<-nest_dists%>%group_by(site)%>%
  summarize(mean_dist=mean(mean_dist),
            min_dist=mean(min_dist))
nest_dists_sum2<-nest_dists2%>%group_by(site,Year)%>%
  summarize(mean_dist=mean(mean_dist),
            min_dist=mean(min_dist),
            lat=mean(latitude))

ggplot(nest_dists_sum2,aes(x=Year,y=min_dist,color=lat))+
  geom_point()

ggplot(nest_dists_sum,aes(x=mean_dist))+
  geom_histogram()


#remove nests that are far away from nests within the same site (likely coordinate errors)
nests2<-nests%>%
  filter(!(id%in%dist_min2[dist_min2$min_dist>1000,]$id))


# write this new filtered nest data to file. This is the final nest dataset used for modeling.
st_write(nests2,paste0(path_out,"Final_outputs/Nest_locations/SALS_nests_2010_2020_dist_err_removed.shp"),delete_layer = T)
