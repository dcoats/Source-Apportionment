setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment")
#save(simevmeans,file="evmeanscv01")
#save(sim1means,file="jagsmeanscv01")
#save(pmfsim,file="pmfsimcv01")
load("evmeanscv01")
load("jagsmeanscv01")
load("pmfsimcv01")
ev01<-simevmeans
jags01<-sim1means
pmf01<-pmfsim
ev01
jags01
pmf01

#save(simevmeans,file="evmeanscv2")
#save(sim1means,file="jagsmeanscv2")
#save(sim2means,file="jagsnopriorscv2")
#save(pmfsim,file="pmfsimcv2")
load("evmeanscv2")
load("jagsmeanscv2")
load("jagsnopriorscv2")
load("pmfsimcv2")
ev2<-simevmeans
jags2<-sim1means
jagsnopriors2<-sim2means
pmf2<-pmfsim
ev2
jags2
jagsnopriors2
pmf2

#save(simevmeans,file="evmeanscv5")
#save(sim1means,file="jagsmeanscv5")
#save(pmfsim,file="pmfsimcv5")
load("evmeanscv5")
load("jagsmeanscv5")
load("pmfsimcv5")
ev5<-simevmeans
jags5<-sim1means
pmf5<-pmfsim
ev5
jags5
pmf5





#########################################################
#########################################################
######              Creating Fleets 2             #######
#########################################################
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment")

#save(sim1means,file="mod2jagsmeanscv2")
load("mod2jagsmeanscv2")
jags2<-sim1means
#save(simevmeans,file="mod2evmeanscv2")
load("mod2evmeanscv2")
ev2<-simevmeans
#save(pmfsim,file="mod2pmfmeanscv2")
load("mod2pmfmeanscv2")
pmf2<-pmfsim
####
jags2
ev2
pmf2

#save(sim1means,file="mod2jagsmeanscv35")
load("mod2jagsmeanscv35")
jags35<-sim1means
#save(pmfsim,file="mod2pmfmeanscv35")
load("mod2pmfmeanscv35")
pmf35<-pmfsim

####
jags35
pmf35


#save(sim1means,file="mod2jagsmeanscv05")
load("mod2jagsmeanscv05")
jags05<-sim1means

####
jags05



##############################################################################
##############################################################################
###########                 Model 3 Summary                      #############
##############################################################################
##############################################################################

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment")

#save(sim1means,file="mod3jagsmeanscv05")
load("mod3jagsmeanscv05")
jags05<-sim1means
#save(simevmeans,file="mod3evmeanscv05")
load("mod3evmeanscv05")
ev05<-simevmeans
#save(pmfsim,file="mod2pmfmeanscv2")
load("mod2pmfmeanscv2")
pmf05<-pmfsim
####
jags05
ev05
pmf05

