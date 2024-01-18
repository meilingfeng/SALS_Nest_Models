

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)

## Plotting spatial data 
library(spData)
library(tmap)
library(viridis)

### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"

## 2. Load nest prediction rasters
pres<-map(unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Placement"),~list.files(.,pattern = "pres_BRTpreds_30mb.tif$",full.names=T))),rast)
surv<-map(unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Success"),~list.files(.,pattern = "surv_BRTpreds_30mb.tif$",full.names=T))),rast)

pres_bi<-map(unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Placement"),~list.files(.,pattern = "pres_BRTbinary_30m_b.tif$",full.names=T))),rast)
surv_bi<-map(unlist(map(paste0(path_out,"Final_outputs/Nest_Predictions/Success"),~list.files(.,pattern = "surv_BRTbinary_30m_b.tif$",full.names=T))),rast)

## 3. SALS priority sites layer - based on high quality high marsh- where populations are stable or increasing and state level site recommendations
priority<-st_read(paste0(dat_path,"SALS_priority_marshes/SALS_Priority_Marshes_Attributed_2022.shp"))%>%
  dplyr::select(id=MarshName,area_m=area)

## 4. SALS extent layer
sf_use_s2(FALSE)
extent<-st_read(paste0(dat_path,"SALS_Breeding_Range/SALS_breeding_19Apr23.shp"))%>%
  st_set_crs(st_crs("EPSG:4269"))%>%
  dplyr::select(id=TARGET_FID,area_m)

## 6. make sure coordinate systems align
priority<-st_transform(priority,crs(pres[[1]]))
extent<-st_transform(extent,crs(pres[[1]]))

# East coast 
data(us_states)
states<-us_states%>%filter(REGION=="Norteast")%>%st_transform(crs(pres[[1]]))


# site summaries
priority_sums<-read.csv(paste0(path_out,"Final_outputs/Nest_Predictions/SALS_Priority_Sites_Prediction_Summary_1_1_24.csv"))
extent_sums<-read.csv(paste0(path_out,"Final_outputs/Nest_Predictions/SALS_Extent_Prediction_Summary_1_1_24.csv"))

  # join to attribute tables
priority<-left_join(priority,priority_sums,by="id")
extent<-left_join(extent,extent_sums,by="id")

tmaptools::palette_explorer()

## 3. Thresholds that maximize sensitivity and specificity
load(paste0(path_out,"Final_outputs/Nest_Predictions/final_BRT_mods.RDS"))




#Summary stats
#---------------------------------------------------------------------
priority<-priority%>%
  mutate(pres_value=ifelse(wavg_pres>=thr.p.brt,1,0),
         pres_coverage=ifelse(pct_pres>=50,1,0),
         pres_best=ifelse(max_pres>=thr.p.brt,1,0))

extent<-extent%>%
  mutate(pres_value=ifelse(wavg_pres>=thr.s.brt,1,0),
         pres_coverage=ifelse(pct_pres>=50,1,0),
         pres_best=ifelse(max_pres>=thr.s.brt,1,0))

100*(sum(priority$pres_value,na.rm=T)/nrow(priority))#percent of priority marshes with average nesting probability above threshold
100*(sum(priority$pres_coverage,na.rm=T)/nrow(priority))#percent of priority marshes with more than 50% nesting area
100*(sum(priority$pres_best,na.rm=T)/nrow(priority))# percent with nesting probability at threshold
mean(priority$wavg_pres,na.rm=T)#average probability across priority marshes

100*(sum(extent$pres_value,na.rm=T)/nrow(extent))#percent of all marshes with average nesting probability above threshold
100*(sum(extent$pres_coverage,na.rm=T)/nrow(extent))#percent of all marshes with more than 50% nesting area
100*(sum(extent$pres_best,na.rm=T)/nrow(extent))# percent with nesting probability at threshold
mean(extent$wavg_pres,na.rm=T)#average probability across all marshes

#total suitable habitat area within priority marshes
p_area<-sum(priority_sums$area_pred_pres,na.rm=T)/10000
#total area of priority marshes
priority_area<-sum(priority_sums$area_m_pres,na.rm=T)/10000
#total suitable habitat area
total_area<-(sum(unlist(map(pres_bi,function(x){as.numeric(global(x==1,sum,na.rm=T))})))*900)/10000
#total habitat area outside priority marshes
op_area<-total_area-p_area

# percent nesting habitat covered by priority areas (hectares)
100*(p_area/total_area);p_area
total_area

# percent of priority areas covered by nesting habitat (hectares)
100*(p_area/priority_area);priority_area

#total predicted nesting habitat
freq(pres_bi[[1]])
freq(pres_bi[[2]])
freq(pres_bi[[3]])
freq(pres_bi[[4]])
freq(pres_bi[[5]])
freq(pres_bi[[6]])
freq(pres_bi[[7]])
freq(pres_bi[[8]])
((49717+6700+29004+23286+65031+6223+6989+6698)*900)/10000

# land ownership within nesting habitat
pad<-vect("D:/Misc_data/PADUS_Land_Ownership/PADUS3_0Combined_Region1.shp")%>%
  terra::project(crs(pres[[1]]))
pad2<-tidyterra::mutate(pad,owner=ifelse(((Own_Type=="PVT"|Own_Type=="UNK") & (Mang_Type=="PVT"|Mang_Type=="UNK"))|GAP_Sts==4,1,0))

pvt_count<-c()
unprot_count<-c()
outs<-list()
for(i in 1:length(pres_bi)){
t<-pres_bi[[i]]
t[t==0]<-NA
pad_r<-rasterize(pad2, t,field="owner",fun=max,background=2)
pad_r2<-mask(pad_r,t)
pvt_count[i]<-ncell(pad_r2[pad_r2==1])
unprot_count[i]<-ncell(pad_r2[pad_r2==2])
outs[[i]]<-pad_r2
}

pct_pvt<-100*((sum(pvt_count)*900)/total_area) #percent of nest habitat with private ownership and management
pct_unprot<-100*((sum(unprot_count)*900)/total_area) #percent of nest habitat with unknown protection or land management
100-(pct_pvt+pct_unprot) # pct of nest habitat on public protected land

#writeRaster(pad_r2,paste0(path_out,"Final_outputs/Nest_Predictions/owner.tif"))

wthn_count<-c()
out_count<-c()
for(i in 1:length(outs)){
within<-mask(outs[[i]],vect(priority))
outside<-mask(outs[[i]],vect(priority),inverse=T)
wthn_count[i]<-ncell(within[within==2|within==1]) #count of unprotected/private habitat cells in priority marshes
out_count[i]<-ncell(outside[outside==2|outside==1]) #count of unprotected/private habitat cells outside priority marshes
}

pct_prior<-100*((sum(wthn_count)*900)/p_area) #percent of priority nest habitat with unknown or private ownership and management
pct_op<-100*((sum(out_count)*900)/op_area) #percent of opportunity nest habitat with unknown or private ownership and management

100-pct_prior;(p_area-(sum(wthn_count)*900))/4046.8627 #percent of priority nest habitat on protected public land;acres
100-pct_op;(sum(out_count)*900)/4046.8627 #percent of opportunity nest habitat on protected public land;acres

(p_area-(sum(wthn_count)*900))/(total_area-((sum(wthn_count)*900)+(sum(out_count)*900)))# percent of public nest habitat in priority marshes
((total_area-((sum(wthn_count)*900)+(sum(out_count)*900)))-(p_area-(sum(wthn_count)*900)))/4046.8627
## Map of observations
#--------------------------------------------------

tmap_mode("view")

tm_basemap("Stamen.TerrainBackground")+
  tm_shape(extent)+
  tm_fill(col="wavg_pres",palette = "viridis")+
  tm_shape(priority%>%filter(pres==1))+
  tm_borders(col="red")

tmap_mode("plot")
tmap_options(bg.color = "skyblue", legend.text.color = "black", inner.margins = c(0, .02, .02, .02))
tmap_style("classic")
land<-st_transform(ne_download(scale = 10, type = "land", category = "physical",returnclass = "sf"),crs(extent))

t<-priority%>%filter(pres==1&!is.na(wavg_pres))
tm_shape(land,bbox = (ext(t[1,])+c(20000,20000,10000,10000)))+
  tm_fill(col="white")+
  tm_shape(extent)+
  tm_fill(col="wavg_pres",palette = "viridis",title = "Average Nesting\nProbability")+
  tm_shape(t)+
  tm_borders(col="red",lwd=2)



## combine mulitple panes

tmap_arrange(map_world + tm_shape(north_america) + tm_fill(), 
             tm_shape(north_america) + tm_fill())

## say we want different colours for different countries
tm_shape(north_america) + tm_fill(col = "name_long")

# to explore the inbuilt colours, look at.
##tmaptools::palette_explorer()

## It is easy to add a compass and scale bar

na_map <- tm_shape(north_america) + 
  tm_fill(col = "name_long", title = "Country name") ## note legend title
na_map + tm_compass(position = c("right", "bottom"), type = "4star") +
  tm_scale_bar()
# We might also want to change the layout (equivalent to theme in ggplot)

na_map + tm_compass(position = c("right", "bottom"), type = "4star") +
  tm_scale_bar() + tm_layout(bg.color = "aliceblue")
## there are also built in themes
na_map + tm_style("classic")

#Map of predictors
#-------------------------------------------------




#variable contributions
#--------------------------------------------