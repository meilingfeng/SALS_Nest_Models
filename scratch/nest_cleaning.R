#Converting DMS to DD
#dat3<-dat2%>%
#  filter(Easting>6700000)%>%
#  mutate(Long=paste0(substr(Easting,1,2),"d",
#                            substr(Easting,3,4),"m",
#                            substr(Easting,5,6),"s",
#                            "W"),
#         Long=(as.numeric(char2dms(Long, chd='d', chm='m', chs='s'))),
#         Lat=paste0(substr(Northing,1,2),"d",
#                     substr(Northing,3,4),"m",
#                     substr(Northing,5,6),"s",
#                     "N"),
#         Lat=as.numeric(char2dms(Lat, chd='d', chm='m', chs='s')))%>%
#  rbind(dat2[!(dat2$Easting>6700000)|is.na(dat2$Easting),])

