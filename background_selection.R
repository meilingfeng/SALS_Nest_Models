library(tidyverse)
library(dismo)
library(terra)
library(sf)


## 2. Select random background points in each zone
#----------------------------------------------------------------------------
# select n random points 
n<-10000
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)
bg <- randomPoints(mask, n, tryf = 3)
#And inspect the results by plotting

# set up the plotting area for two maps
par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)
