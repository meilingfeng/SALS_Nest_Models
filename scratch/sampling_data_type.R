# First set the resampling method and make the output data type match the input data type
#####
#if the data is numeric (decimals), use average resampling
#if the data is nominal (classes), use mode resampling (only the veg class layer)
# on spat raster data types https://rdrr.io/cran/raster/man/dataType.html 
if(grepl("INT",rast(input)@cpp[["dataType"]],fixed=T)){
  rs_method<-"mode"
  dat_type<-"UInt16"
}
if(grepl("FLT",rast(input)@cpp[["dataType"]])){
  rs_method<-"average"
  dat_type<-"Float32"
}
