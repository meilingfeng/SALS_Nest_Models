# Visual of Discrimination evaluation (important if making presence absence maps) - histogram of presences to absences 
library(dismo)
library(gbm)
par(mfrow=c(3,2))

for (mod in 1:length(model.names)){
  
  ifelse(mod==3, leg<-T,leg<-F)
  
  presence.absence.hist(d.pres[[5]], which.model=mod,add.legend=leg,legend.cex=0.6,N.bars=15, opt.methods=c("MaxSens+Spec","MaxKappa"),add.opt.legend = leg,main="Presence")
  
  presence.absence.hist(d.surv[[5]], which.model=mod,add.legend=leg,legend.cex=0.6,N.bars=15, opt.methods=c("MaxSens+Spec","MaxKappa"),add.opt.legend = leg,main="Survival")
  mtext(model.names[mod], side = 2, line = 26.3, cex = 1.6)
}
