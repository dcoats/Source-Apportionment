##
##  METALS ONLY -- LOOK AT VARYING AMOUNTS OF ERRORS ON UNCERTAINTIES
##

## RE-DONE ANALYSIS WITH STLOUIS DATA (LHT version)

lhtdat <- read.table("LHTconc.txt",header=T)
lhtunc <- read.table("LHTunc.txt",header=T)

jaydat <- read.table("../receptor/FromJay20071019/ToWilliam/PMcoarse_input_conc.txt",header=T)
jayunc <- read.table("../receptor/FromJay20071019/ToWilliam/PMcoarse_input_unc.txt",header=T)

#signoiselht <- apply(lhtdat/lhtunc,2,mean)
#signoisejay <- apply(jaydat[,-1]/jayunc[,-1],2,mean)
#signoiselht
#signoisejay

dropvals <- c(9:14,18,21,23,25,26,31,33,35,36,38,39,43,45,46,49,52,53,55,58,59)

#signoisejay[-dropvals] / signoiselht[-1]
## ADJUST UNC in LHT using the ratio above, then re-run .... DIDN'T AFFECT MUCH ##




#newdates <- jaydat[,1]
#datesnew <- as.Date(newdates, "%b/%d/%Y")
#datesnew
lhtdat2 <- jaydat[,-1]
lhtunc2 <- jayunc[,-1]

lhtdat2 <- lhtdat2[,-c(1,dropvals)]
lhtunc2 <- lhtunc2[,-c(1,dropvals)]

#lhtunc2 <- lhtunc2 / matrix( rep( (signoisejay[-dropvals] / signoiselht[-1]) , nrow(lhtunc2) ),
#                             nrow(lhtunc2), ncol(lhtunc2), byrow=T)
#signoiselht <- apply(lhtdat2/lhtunc2,2,mean)
#signoisejay <- apply(jaydat[,-1]/jayunc[,-1],2,mean)
#signoiselht
#signoisejay

#signoisejay[-dropvals] / signoiselht



write(t(lhtdat2), file="lhtdat2.txt", ncol=ncol(lhtdat2))
write(t(lhtunc2), file="lhtunc2.txt", ncol=ncol(lhtunc2))

mynames <- c( "OC1" ,"OC2", "OC3", "OC4", "OP", "EC1", "EC2", "EC3", "SO4",  "NO3", 
 "NH4",  "AL", "AS", "BA", "CA", "CO", "CR", "CU", "FE", "HG",
 "K" ,"MN", "NI" ,"P", "PB", "RB" ,"SE", "SI", "SR", "TI",
 "V", "ZN", "ZR")
 

numfacs <- i <- 10
  concs <- lhtdat2
  uncs <- lhtunc2

#mynames <- dimnames(morgc)[[2]]
#mynames
names <- mynames
#newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
#newdates

runpmfplain("default_nogkey.ini","lhtdat2.txt", "lhtunc2.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

lamORIG<-read.table(mypaste("testlam",numfacs,".txt"))
facORIG<-read.table(mypaste("testfac",numfacs,".txt"))

lam <- lamORIG
fac <- facORIG

apply(facORIG,2,mean)
sum(apply(facORIG,2,mean))
 apply(facORIG,2,mean) / sum(apply(facORIG,2,mean))

apply(facTEMP,2,mean)
sum(apply(facTEMP,2,mean))
 apply(facTEMP,2,mean) / sum(apply(facTEMP,2,mean))

#           mass        %       LHT%   LHTFkey%   jaydat
#1: lead    .13        .011     .013   .008       .011 lead      .011
#2: copper  .16        .013     .005   .011       .3457 sulfate  .3464
#3: sulfate  4.65      .380     .33    .375       .015 copper    same below
#4: gas        .65     .053     .16    .057       .1106 
#5: zinc    .16        .013     .013   .009       .0159 zinc
#6: C-rich SO4 1.79    .146     .20               .1111
#7: soil    .39        .032     .04               .0337 soil
#8: diesel  .97        .079     .02               .06203 steel
#9: steel   .57        .046     .07               .0387 diesel
#10: nitrate 2.79      .228     .15               .2549 nitrate
## NOTE after tweaking ion uncertainties, #2 and #3 switch

# For Fkey, set row 9 columns 1, 2, 5 to the value of 5
# For Fkey, set row 9 columns 1, 3, 5 to the value of 5  --- TWEAKED

myFkey <- matrix(0,33,10)
#myFkey[9,c(1,2,5)] <- 5
myFkey[9,c(1,3,5)] <- 5  ## TWEAKED
myFkey

write(t(myFkey),file="myFkey.txt",ncol=ncol(myFkey))
#write( (myFkey),file="myFkey.txt",ncol=nrow(myFkey))
write(t(facORIG),file="startfac.txt",ncol=ncol(facORIG))
write(t(lamORIG),file="startlam.txt",ncol=ncol(lamORIG))


runpmfplain("default_Fkey.ini","lhtdat2.txt", "lhtunc2.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    
lamTEMP<-read.table(mypaste("testlam",numfacs,".txt"))
facTEMP<-read.table(mypaste("testfac",numfacs,".txt"))

lam <- lamTEMP
fac <- facTEMP

apply(facTEMP,2,mean)
sum(apply(facTEMP,2,mean))
 apply(facTEMP,2,mean) / sum(apply(facTEMP,2,mean))





par(mfrow=c(2,1))
for (j in 1:2) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 3:4) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 5:6) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 7:8) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 9:10) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")

par(mfrow=c(2,1))
for (j in 1:2) plot(1:nrow(fac),fac[,j],type="l")
for (j in 3:4) plot(1:nrow(fac),fac[,j],type="l")
for (j in 5:6) plot(1:nrow(fac),fac[,j],type="l")
for (j in 7:8) plot(1:nrow(fac),fac[,j],type="l")
for (j in 9:10) plot(1:nrow(fac),fac[,j],type="l")


lamplots2old(lam,"gkey",names=mynames) 
facplots3old(fac,"gkeyf",newdates)
facplots(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 



lamplotslogps(lamORIG,"gkey",names=mynames,facnames[reorderfacs]) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)


facnames9 <- c("Lead\nSmelter","Copper\nSmelter","Fireworks","Summer\nSecondary","Zinc\nSmelter",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecondary")
facnames9 <- c("Lead\nSmelt.","Copper\nSmelt.","Fire-\nworks","Summer\nSecond.","Zinc\nSmelt.",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecond.")
facnames9MORG  <- c("Mobile\n(high Cu)","Zinc\nSmelt.","Summer\nSecond.","Wood1","Wood2",
                    "Mobile1","SO/S/\nNa/EC","Steel\nMill","Winter\nSecond.")



cvvals <- c(.01,.05,.15)
reperrorvals <- 50

lamarrayerrsM <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsM <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

for (cvi in 1:length(cvvals))
{
  for (reperri in 1:reperrorvals)
  #for (reperri in 15:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    if (cvi> 1) uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
                                          nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","sldat2.txt", "uncerrs.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsM[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsM[,,cvi,reperri] <- as.matrix(fac)
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
faccompsall1M <- matrix(NA,reperrorvals,numfacs)
lamcompsall1M <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse1M <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp1M <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall1M[i,] <- bestcormatcors(facarrayerrsM[,,cvi,i],facORIG)
  lamcompsall1M[i,] <- bestrevcormatcors(lamarrayerrsM[,,cvi,i],lamORIG)
  faccompsallrmse1M[i,] <- bestrmsematcors(facarrayerrsM[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp1M[i,] <- bestexpcormatcors(lamarrayerrsM[,,cvi,i],lamORIG,facarrayerrsM[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1M))
boxplot(data.frame(lamcompsall1M))

date()

cvi <- 2
faccompsall2M <- matrix(NA,reperrorvals,numfacs)
lamcompsall2M <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse2M <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp2M <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall2M[i,] <- bestcormatcors(facarrayerrsM[,,cvi,i],facORIG)
  lamcompsall2M[i,] <- bestrevcormatcors(lamarrayerrsM[,,cvi,i],lamORIG)
  faccompsallrmse2M[i,] <- bestrmsematcors(facarrayerrsM[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp2M[i,] <- bestexpcormatcors(lamarrayerrsM[,,cvi,i],lamORIG,facarrayerrsM[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2M))
boxplot(data.frame(lamcompsall2M))


date()

cvi <- 3
faccompsall3M <- matrix(NA,reperrorvals,numfacs)
lamcompsall3M <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse3M <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp3M <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall3M[i,] <- bestcormatcors(facarrayerrsM[,,cvi,i],facORIG)
  lamcompsall3M[i,] <- bestrevcormatcors(lamarrayerrsM[,,cvi,i],lamORIG)
  faccompsallrmse3M[i,] <- bestrmsematcors(facarrayerrsM[,,cvi,i],facORIG) ## now does aae
  lamcompsallexp3M[i,] <- bestexpcormatcors(lamarrayerrsM[,,cvi,i],lamORIG,facarrayerrsM[,,cvi,i],facORIG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3M))
boxplot(data.frame(lamcompsall3M))






par(mfcol=c(2,3))
boxplot(data.frame(faccompsall2M),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsall2M),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall3M),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsall3M),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall1M),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsall1M),ylim=c(0,1),xlab="Profiles")


pdf("MetalsUncertainty.pdf",width=12,height=9)
set.cex <- .65
par(mfcol=c(2,3))
maxy <- max(c( max(data.frame(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) ))  )) 
#boxplot(data.frame(faccompsallrmse2M),ylim=c(0,2.5),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
title(paste("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp2M),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp2M)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse3),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
title(paste("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp3M),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp3M)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse1M),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse1M)/apply(facUNWT,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
title(paste("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp1M),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp1M)) )) ,3)    ))
dev.off()



newmeans1M <- matrix(NA,reperrorvals,numfacs)
newmeans2M <- matrix(NA,reperrorvals,numfacs)
newmeans3M <- matrix(NA,reperrorvals,numfacs)
for (reperri in 1:reperrorvals)
{
  matchbool <- getaae(facarrayerrsM[,,1,reperri],facORIG)==faccompsallrmse1M[reperri,]
  newmeans1M[reperri,] <- apply(facarrayerrsM[,,1,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsM[,,2,reperri],facORIG)==faccompsallrmse2M[reperri,]
  newmeans2M[reperri,] <- apply(facarrayerrsM[,,2,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
  matchbool <- getaae(facarrayerrsM[,,3,reperri],facORIG)==faccompsallrmse3M[reperri,]
  newmeans3M[reperri,] <- apply(facarrayerrsM[,,3,reperri],2,mean)[apply(matchbool*(1:9),2,sum)]
}

set.cex <- 1.2

#pdf("MetalsContributions.pdf",width=12,height=5)
#par(mfrow=c(1,3))
pdf("MetalsContributions2.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
junk <- apply(facarrayerrsM[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)
boxplot(data.frame(newmeans2M), xaxt="n", ylim=c(0,ymaxcont),cex=set.cex)  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .05*true)"))
title("(a)")
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)
dev.off()

pdf("MetalsContributions3.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans3M), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .15*true)"))
title("(b)")
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)
dev.off()

pdf("MetalsContributions1.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans1M), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
title("(c)")
points(1:9,apply(facORIG,2,mean),pch=4,cex=2)
dev.off()



maxy <- max(c( max(data.frame(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) ))  )) 

pdf("MetalsRAAE2.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse2M)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(a)")
dev.off()


pdf("MetalsRAAE3.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(b)")
dev.off()


pdf("MetalsRAAE1.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse1M)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(c)")
dev.off()



allmin <- min(c(c(as.matrix(t(t(lamcompsallexp2M)) )), 
                c(as.matrix(t(t(lamcompsallexp3M)) )), 
                c(as.matrix(t(t(lamcompsallexp1M)) ))))

pdf("MetalsEMC2.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(lamcompsallexp2M),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp2M)) )), col="gray", lwd=2)
title("(a)")
dev.off()


pdf("MetalsEMC3.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(lamcompsallexp3M),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp3M)) )), col="gray", lwd=2)
title("(b)")
dev.off()


pdf("MetalsEMC1.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
set.cex <- 1.2
boxplot(data.frame(lamcompsallexp1M),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp1M)) )), col="gray", lwd=2)
title("(c)")
dev.off()


