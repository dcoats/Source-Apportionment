setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment")
save(simevmeans,file="evmeans5low")
save(sim2means,file="jagsmeans5low")
save(pmfsim,file="pmfsim5low")
load("evmeans5low")
load("jagsmeans5low")
load("pmfsim5low")
evlow<-simevmeans
jagslow<-sim2means
pmflow<-pmfsim
evlow
jagslow
pmflow