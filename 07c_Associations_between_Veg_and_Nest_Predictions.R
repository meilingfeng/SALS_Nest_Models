

library(tidyverse)
library(terra)
library(sf)



###########################################################################################
# look for associations between:
# 1. vegetation characteristics and nesting probabilities (nest placement and nest success)
#    Do females select for certain veg characteristics associated with nest success?
# 2. nest characteristics and nest success probability
#    Do the choices females make after selecting a nesting site aid in nest success?
###########################################################################################

# if we don't find evidence for these, are extreme flooding events driving success and masking any historical habitat drivers?
# do we see a trend in the quality of nesting habitat over time? Changes in nesting niches over time? 
  # is a change in niche due to adaptation, still within the range limits of the species, or a result of available habitat?
  # How much salt marsh/high marsh have we lost over this period? Are sparrows losing suitable habitat?



### Set up
# -------------------------------------------
## 1. file paths
dat_path<-"D:/Nest_Models/Data/"
path_out<-"D:/Nest_Models/Outputs/"



### Surrounding Vegetation and Nest Predictions
#--------------------------------------------------------
# What habitat variables influence where females choose nest locations?

# Run an analysis with all data points, looking at surrounding veg characteristics

# look at fine resolution characteristics (Demo veg data-1m) and larger resolution characteristics (survey data-50m)
# expect a stronger relationship with larger resolution data since predictions seem to distinguish larger, site level differences

# 1. Correlation of continuous variables with predicted nest probabilites




# 2. scatter plots (continuous variables)/boxplots (categorical variables)
#ggplot(dat_all%>%,aes(x=presence,y=Geom_histogram() + facet_wrap() 
                      
#ggplot(dat_all%>%pivot_longer(cols = c(tll_mode,dom_pc,alt,patens,juncus,distichlis), names_to = "var", values_to = ""),
#                             aes(x=presence, fill =)+Geom_density(alpha=.3) 
                             
#aes(x=var1,y=var2)+geom_boxplot() or + geom_jitter()-jitter just adds the points 
                             
# 3. Run a principal components analysis on the percent cover variables? See if we can condense some of the species into fewer variables.
                             

                             
                             
                             
                             

### Nest characteristics and Nest Predictions
#--------------------------------------------------------
# What variables might be important for nest survival after the female chooses a nest location?
# Run an analysis with just nest points, looking at nest construction variables

# 1. Correlation of continuous variables with predicted nest probabilites




# 2. scatter plots (continuous variables)/boxplots (categorical variables)
#ggplot(dat_all%>%,aes(x=presence,y=Geom_histogram() + facet_wrap() 

#ggplot(dat_all%>%pivot_longer(cols = c(tll_mode,dom_pc,alt,patens,juncus,distichlis), names_to = "var", values_to = ""),
#       aes(x=presence, fill =)+Geom_density(alpha=.3) 
       
#       Aes(x=var1,y=var2)+geom_boxplot() or/+ geom_jitter()-jitter just adds the points 

       
# 3. Run a principal components analysis on the percent cover variables? See if we can condense some of the species into fewer variables.


# 4. regression for predicted probabilities (also use binary pres/abs predictions with ANOVA?)
      # choose a global model with all variables
       #choose subsets of the global model
      #compare sets for the 2 resolutions of veg data (survey and demo) and the nest data (but focus on just doing survey veg data for now)
