#Wallace vignette https://wallaceecomod.github.io/vignettes/wallace_vignette.html
#Wallace website  https://wallaceecomod.github.io/

# To use Maxent in Wallace
# download maxent program https://biodiversityinformatics.amnh.org/open_source/maxent/
# locate maxent.jar, which is the Maxent program itself, in the downloaded folder
# install.packages('dismo')
# find the directory path to dismo/java by running system.file('java', package="dismo") 



# load the package
library(wallace)
# run the app
run_wallace()

library(dismo)
vignette('sdm', 'dismo')

#install.packages(c('raster','rgdal','dismo','rJava'))

