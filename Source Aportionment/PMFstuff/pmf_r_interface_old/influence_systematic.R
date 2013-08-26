




##
##  MORG 9 SOURCE -- LOOK AT VARYING AMOUNTS OF ERRORS ON UNCERTAINTIES
##


  concs <- morgc
  uncs <- morgu
  mynamesmorg <- dimnames(morgc)[[2]]
  #names <- mynames
  newdatesmorg <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
  newdatesmorg


numfacs <- i <- 9
#datesnew <- getdates()
#newdates <- datesnew[apply(sldat==-999.9,1,sum)==0]
#mynames<-c("Na","Mg","Al","Si","Ph",
#		"S","Cl","K","Ca","Ti",
#		"V","Cr","Mn","Fe","Co",
#		"Ni","Cu","Zn","Ga","As",
#		"Se","Br","Rb","Sr","Y",
#		"Zr","Mo","Pd","Ag","Cd",
#		"In","Sn","Sb","Ba",
#		"La","Au","Hg","Tl","Pb",
#		"U","OC","EC","SO","NO")

#  write(t(morgc),file="morgc.txt",ncol=ncol(morgc))
#  write(t(morgu),file="morgu.txt",ncol=ncol(morgu))


runpmfplain("default_nogkey.ini","morgc.txt", "morgu.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    
lamORIGMORG<-read.table(mypaste("testlam",numfacs,".txt"))
facORIGMORG<-read.table(mypaste("testfac",numfacs,".txt"))


cvvals <- c(.01,.05,.15)
reperrorvals <- 50

lamarrayerrsMORG <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsMORG <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

for (cvi in 2:length(cvvals))
{
  startrep <- 1
  #
  startrep <- reperri
  startrep
  for (reperri in startrep:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    if (cvi> 1) uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
                                          nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","morgc.txt", "uncerrs.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsMORG[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsMORG[,,cvi,reperri] <- as.matrix(fac)
     reperri
  }
}

date()

cvi <- 1
faccompsall1MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsall1MORG <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse1MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp1MORG <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall1MORG[i,] <- bestcormatcors(facarrayerrsMORG[,,cvi,i],facORIGMORG)
  lamcompsall1MORG[i,] <- bestrevcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG)
  faccompsallrmse1MORG[i,] <- bestrmsematcors(facarrayerrsMORG[,,cvi,i],facORIGMORG) ## now does aae
  lamcompsallexp1MORG[i,] <- bestexpcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG,
                                               facarrayerrsMORG[,,cvi,i],facORIGMORG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1MORG))
boxplot(data.frame(lamcompsall1MORG))

date()

cvi <- 2
faccompsall2MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsall2MORG <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse2MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp2MORG <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall2MORG[i,] <- bestcormatcors(facarrayerrsMORG[,,cvi,i],facORIGMORG)
  lamcompsall2MORG[i,] <- bestrevcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG)
  faccompsallrmse2MORG[i,] <- bestrmsematcors(facarrayerrsMORG[,,cvi,i],facORIGMORG) ## now does aae
  lamcompsallexp2MORG[i,] <- bestexpcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG,
                                               facarrayerrsMORG[,,cvi,i],facORIGMORG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2MORG))
boxplot(data.frame(lamcompsall2MORG))

date()

cvi <- 3
faccompsall3MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsall3MORG <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse3MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp3MORG <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall3MORG[i,] <- bestcormatcors(facarrayerrsMORG[,,cvi,i],facORIGMORG)
  lamcompsall3MORG[i,] <- bestrevcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG)
  faccompsallrmse3MORG[i,] <- bestrmsematcors(facarrayerrsMORG[,,cvi,i],facORIGMORG) ## now does aae
  lamcompsallexp3MORG[i,] <- bestexpcormatcors(lamarrayerrsMORG[,,cvi,i],lamORIGMORG,
                                               facarrayerrsMORG[,,cvi,i],facORIGMORG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3MORG))
boxplot(data.frame(lamcompsall3MORG))

date()






par(mfcol=c(2,3))
boxplot(data.frame(faccompsall2MORG),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsall2MORG),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall3MORG),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsall3MORG),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall1MORG),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsall1MORG),ylim=c(0,1),xlab="Profiles")


facnames9MORG  <- c("Mobile2\n(high Cu)","Zinc\nSmelt.","Summer\nSecond.","Wood1","Wood2",
                    "Mobile1","SO/S/\nNa/EC","Steel\nMill","Winter\nSecond.")

pdf("MORGUncertainty.pdf",width=12,height=9)
set.cex <- .65
par(mfcol=c(2,3))
maxy <- max(c( max(data.frame(t(t(faccompsallrmse2MORG)/apply(facORIGMORG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3MORG)/apply(facORIGMORG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1MORG)/apply(facORIGMORG,2,mean)) ))  )) 
boxplot(data.frame(t(t(faccompsallrmse2MORG)/apply(facORIGMORG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse2MORG)/apply(facORIGMORG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp2MORG),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp2MORG)) )) ,3)    ))

boxplot(data.frame(t(t(faccompsallrmse3MORG)/apply(facORIGMORG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
title(paste("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse3MORG)/apply(facORIGMORG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp3MORG),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp3MORG)) )) ,3)    ))

boxplot(data.frame(t(t(faccompsallrmse1MORG)/apply(facORIGMORG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
#title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
title(paste("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse1MORG)/apply(facORIGMORG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp1MORG),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp1MORG)) )) ,3)    ))
dev.off()



newmeans1MORG <- matrix(NA,reperrorvals,numfacs)
newmeans2MORG <- matrix(NA,reperrorvals,numfacs)
newmeans3MORG <- matrix(NA,reperrorvals,numfacs)
for (reperri in 1:reperrorvals)
{
  matchbool <- getaae(facarrayerrsMORG[,,1,reperri],facORIGMORG)==faccompsallrmse1MORG[reperri,]
  newmeans1MORG[reperri,] <- apply(facarrayerrsMORG[,,1,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsMORG[,,2,reperri],facORIGMORG)==faccompsallrmse2MORG[reperri,]
  newmeans2MORG[reperri,] <- apply(facarrayerrsMORG[,,2,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsMORG[,,3,reperri],facORIGMORG)==faccompsallrmse3MORG[reperri,]
  newmeans3MORG[reperri,] <- apply(facarrayerrsMORG[,,3,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
}


pdf("MORGContributions.pdf",width=12,height=5)
par(mfrow=c(1,3))
junk <- apply(facarrayerrsMORG[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)
boxplot(data.frame(newmeans2MORG), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .05*true)"))
points(1:9,apply(facORIGMORG,2,mean),col="red",pch="-",cex=3)

boxplot(data.frame(newmeans3MORG), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .15*true)"))
points(1:9,apply(facORIGMORG,2,mean),col="red",pch="-",cex=3)

boxplot(data.frame(newmeans1MORG), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9MORG,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
points(1:9,apply(facORIGMORG,2,mean),col="red",pch="-",cex=3)
dev.off()


#WORK AREA#####
#
faccompsallrmse1MORG[1,]
getaae(facarrayerrsMORG[,,1,1],facORIGMORG)
matchbool <- getaae(facarrayerrsMORG[,,1,1],facORIGMORG)==faccompsallrmse1MORG[1,]
apply(facarrayerrsMORG[,,1,1],2,mean)[apply(matchbool*(1:9),2,sum)]
###############




























##
##  METALS ONLY -- LOOK AT VARYING AMOUNTS OF ERRORS ON UNCERTAINTIES
##


numfacs <- i <- 9
datesnew <- getdates()
newdates <- datesnew[apply(sldat==-999.9,1,sum)==0]
mynames<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")

#  write(t(morgc),file="morgc.txt",ncol=ncol(morgc))
#  write(t(morgu),file="morgu.txt",ncol=ncol(morgu))
  concs <- sldat2
  uncs <- slunc2

#mynames <- dimnames(morgc)[[2]]
#mynames
names <- mynames
#newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
#newdates

runpmfplain("default_nogkey.ini","sldat2.txt", "slunc2.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    
lamORIG<-read.table(mypaste("testlam",numfacs,".txt"))
facORIG<-read.table(mypaste("testfac",numfacs,".txt"))

facnames9 <- c("Lead\nSmelter","Copper\nSmelter","Fireworks","Summer\nSecondary","Zinc\nSmelter",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecondary")
facnames9 <- c("Lead\nSmelt.","Copper\nSmelt.","Fire-\nworks","Summer\nSecond.","Zinc\nSmelt.",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecond.")
facnames9MORG  <- c("Mobile\n(high Cu)","Zinc\nSmelt.","Summer\nSecond.","Wood1","Wood2",
                    "Mobile1","SO/S/\nNa/EC","Steel\nMill","Winter\nSecond.")



cvvals <- c(.01,.05,.15)
reperrorvals <- 50

lamarrayerrsMS <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsMS <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

for (cvi in 1:length(cvvals))
{
  for (reperri in 1:reperrorvals)
  #for (reperri in 15:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    if (cvi> 1) uncerrs <- uncs * matrix(rnorm(ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs),byrow=T)
    #if (cvi> 1) uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    if (cvi==1) uncerrs <- uncs * matrix( runif(ncol(uncs),.5,1.5) ,
                                          nrow(uncs),ncol(uncs),byrow=T) #cvi=1
    #if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
    #                                      nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","sldat2.txt", "uncerrs.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsMS[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsMS[,,cvi,reperri] <- as.matrix(fac)
  }
  reperri
}

## Note cvi of 1 is unif(.5,1.5), cvi of 2 is norm(1,.05), cvi of 2 is norm(1, .15)

#cor(facarrayerrs[,,3,1],facarrayerrs[,,3,2])
## Compare facs and lams across time series
#faccomps2 <- matrix(NA,reperri,reperri)
#for (i in 2:(reperri)) for (j in 1:(i-1)) faccomps2[i,j] <- bestcormat(facarrayerrs[,,cvi,i],facarrayerrs[,,cvi,j])
#lamcomps2 <- matrix(NA,reperri,reperri)
#for (i in 2:(reperri)) for (j in 1:(i-1)) lamcomps2[i,j] <- bestrevcormat(lamarrayerrs[,,cvi,i],lamarrayerrs[,,cvi,j])
#
#round(faccomps,3)
#round(lamcomps,3)
date()

cvi <- 1
faccompsall1MS <- matrix(NA,reperrorvals,numfacs)
lamcompsall1MS <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse1MS <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp1MS <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall1MS[i,] <- bestcormatcors(facarrayerrsMS[,,cvi,i],facORIG)
  lamcompsall1MS[i,] <- bestrevcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG)
  faccompsallrmse1MS[i,] <- bestrmsematcors(facarrayerrsMS[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp1MS[i,] <- bestexpcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG,facarrayerrsMS[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1MS))
boxplot(data.frame(lamcompsall1MS))

date()

cvi <- 2
faccompsall2MS <- matrix(NA,reperrorvals,numfacs)
lamcompsall2MS <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse2MS <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp2MS <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall2MS[i,] <- bestcormatcors(facarrayerrsMS[,,cvi,i],facORIG)
  lamcompsall2MS[i,] <- bestrevcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG)
  faccompsallrmse2MS[i,] <- bestrmsematcors(facarrayerrsMS[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp2MS[i,] <- bestexpcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG,facarrayerrsMS[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2MS))
boxplot(data.frame(lamcompsall2MS))


date()

cvi <- 3
faccompsall3MS <- matrix(NA,reperrorvals,numfacs)
lamcompsall3MS <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse3MS <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp3MS <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall3MS[i,] <- bestcormatcors(facarrayerrsMS[,,cvi,i],facORIG)
  lamcompsall3MS[i,] <- bestrevcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG)
  faccompsallrmse3MS[i,] <- bestrmsematcors(facarrayerrsMS[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp3MS[i,] <- bestexpcormatcors(lamarrayerrsMS[,,cvi,i],lamORIG,facarrayerrsMS[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3MS))
boxplot(data.frame(lamcompsall3MS))






par(mfcol=c(2,3))
boxplot(data.frame(faccompsall2MS),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsall2MS),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall3MS),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsall3MS),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall1MS),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsall1MS),ylim=c(0,1),xlab="Profiles")


pdf("MetalsUncertaintyS.pdf",width=12,height=9)
set.cex <- .65
par(mfcol=c(2,3))
maxy <- max(c( max(data.frame(t(t(faccompsallrmse2MS)/apply(facORIG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3MS)/apply(facORIG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1MS)/apply(facORIG,2,mean)) ))  )) 
#boxplot(data.frame(faccompsallrmse2MS),ylim=c(0,2.5),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse2MS)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
title(paste("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse2MS)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp2MS),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp2MS)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse3),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse3MS)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
title(paste("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse3MS)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp3MS),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp3MS)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse1M),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse1M)/apply(facUNWT,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse1MS)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
title(paste("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse1MS)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp1MS),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp1MS)) )) ,3)    ))
dev.off()



newmeans1MS <- matrix(NA,reperrorvals,numfacs)
newmeans2MS <- matrix(NA,reperrorvals,numfacs)
newmeans3MS <- matrix(NA,reperrorvals,numfacs)
for (reperri in 1:reperrorvals)
{
  matchbool <- getaae(facarrayerrsMS[,,1,reperri],facORIG)==faccompsallrmse1MS[reperri,]
  newmeans1M[reperri,] <- apply(facarrayerrsMS[,,1,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsMS[,,2,reperri],facORIG)==faccompsallrmse2MS[reperri,]
  newmeans2M[reperri,] <- apply(facarrayerrsMS[,,2,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsMS[,,3,reperri],facORIG)==faccompsallrmse3MS[reperri,]
  newmeans3M[reperri,] <- apply(facarrayerrsMS[,,3,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
}


pdf("MetalsContributions.pdf",width=12,height=5)
par(mfrow=c(1,3))
junk <- apply(facarrayerrsMS[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)
boxplot(data.frame(newmeans2MS), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .05*true)"))
title("(a)")
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)

boxplot(data.frame(newmeans3MS), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .15*true)"))
title("(b)")
#points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)

boxplot(data.frame(newmeans1MS), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
title("(c)")
#points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)
dev.off()


pdf("MetalsContributions.pdf",width=12,height=5)
par(mfrow=c(1,3))
junk <- apply(facarrayerrsMS[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)
boxplot(data.frame(t( apply(facarrayerrsMS[,,2,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .05*true)"))
points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)

boxplot(data.frame(t( apply(facarrayerrsMS[,,3,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .15*true)"))
points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)

boxplot(data.frame(t( apply(facarrayerrsMS[,,1,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)
dev.off()






