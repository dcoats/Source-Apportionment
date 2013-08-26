##
##  METALS ONLY -- LOOK AT VARYING AMOUNTS OF ERRORS ON UNCERTAINTIES
##

## RE-DONE ANALYSIS WITH STLOUIS DATA (LHT version)

lhtdat <- read.table("LHTconc.txt",header=T)
lhtunc <- read.table("LHTunc.txt",header=T)

#jaydat <- read.table("../receptor/FromJay20071019/ToWilliam/PMcoarse_input_conc.txt",header=T)
#jayunc <- read.table("../receptor/FromJay20071019/ToWilliam/PMcoarse_input_unc.txt",header=T)

#signoiselht <- apply(lhtdat/lhtunc,2,mean)
#signoisejay <- apply(jaydat[,-1]/jayunc[,-1],2,mean)
#signoiselht
#signoisejay

#dropvals <- c(9:14,18,21,23,25,26,31,33,35,36,38,39,43,45,46,49,52,53,55,58,59)

#signoisejay[-dropvals] / signoiselht[-1]
## ADJUST UNC in LHT using the ratio above, then re-run .... DIDN'T AFFECT MUCH ##




newdates <- lhtdat[,1]
datesnew <- as.Date(newdates, "%d%b%Y")

lhtdat2 <- lhtdat[,-1]
lhtunc2 <- lhtunc[,-1]

#lhtunc2 <- lhtunc2 / matrix( rep( (signoisejay[-dropvals] / signoiselht[-1]) , nrow(lhtunc2) ),
#                             nrow(lhtunc2), ncol(lhtunc2), byrow=T)
#signoiselht <- apply(lhtdat2/lhtunc2,2,mean)
#signoisejay <- apply(jaydat[,-1]/jayunc[,-1],2,mean)
#signoiselht
#signoisejay

#signoisejay[-dropvals] / signoiselht

uncs <- lhtunc2
concs <- lhtdat2



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

%Latest result: 11/28/07
             mimic   LHT     NC
#1   lead    .011    .013    .011
#2  sulfate  .369    .33     .360
#3  copper   .012    .005    .016
#4  gas      .075    .16     .133
#5  zinc     .014    .013    .011
#6  C-r. sul .142    .20     .143
#7  soil     .032    .04     .032
#8  steel    .051    .07     .042
#9  diesel   .066    .02     .040
#10 nitrate  .228    .15     .213

compdata <- matrix(
  c( .33 , .20 , .16 , .15 , .07 , .04 , .02 , .013, .013, .005, 
     .369, .142, .075, .228, .051, .032, .066, .014, .011, .012),
      2,10, byrow=T)
facnames10reorder <- 
 c("Secondary\nSulfate","Carbon-rich\nSulfate","Gasoline\nExhaust",
   "Secondary\nNitrate","Steel\nProcessing","Soil","Diesel/\nrailroad",
   "Zinc\nSmelting","Lead\nSmelting","Copper\nSmelting")

#names(compdata) <- list(facnames10reorder, c("LHT2006","Current"))
compdata


### Make Figure comparing with LHT ###
set.cex <- 1
pdf("CompareWithLHT.pdf",width=10,height=6)
par(mar=c(6,3,2,2))
 barplot(compdata, beside = TRUE,
         col = c("black", "lightgray"),
         legend = c("Lee, Hopke, and Turner (2006)","current"), 
         ylim = c(0, .4), yaxt="n")
par(mgp=c(3,.4,0))
axis(1,labels=facnames10reorder,at=seq(2,29,3),cex.axis=set.cex,las=2,tick=F)
par(mgp=c(3,1,0))
axis(2,labels=c("0%","10%","20%","30%","40%"),at=seq(0,.4,.1),las=1)
abline(h=0)
dev.off()


############################
### FKEYING BELOW IGNORE ###
############################

apply(facTEMP,2,mean)
sum(apply(facTEMP,2,mean))
 apply(facTEMP,2,mean) / sum(apply(facTEMP,2,mean))

#           mass        %       LHT%   Fkey%   %(rev)
#1: lead    .13        .011     .013   .008    .011 1 lead
#2: copper  .16        .013     .005   .011    .013 3
#3: sulfate  4.65      .380     .33    .375    .369 2 sulfate
#4: gas        .65     .053     .16    .057    .054 4
#5: zinc    .16        .013     .013   .009    .013 5
#6: C-rich SO4 1.79    .146     .20            .150 6
#7: soil    .39        .032     .04            .032 7
#8: diesel  .97        .079     .02            .082 9
#9: steel   .57        .046     .07            .046 8
#10: nitrate 2.79      .228     .15            .228 10
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

############################
############################
############################




par(mfrow=c(2,1))
for (j in 8:9) {plot(1:nrow(lam),lam[,j],xaxt="n") 
           axis(1,at=1:nrow(lam),labels=mynames,cex.axis=.4)}

par(mfrow=c(2,1))

for (j in 1:2) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 3:4) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 5:6) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 7:8) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
for (j in 9:10) plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")

mynames <- c( "OC1" ,"OC2", "OC3", "OC4", "OP", "EC1", "EC2", "EC3", "SO4",  "NO3", 
 "NH4",  "AL", "AS", "BA", "CA", "CO", "CR", "CU", "FE", "HG",
 "K" ,"MN", "NI" ,"P", "PB", "RB" ,"SE", "SI", "SR", "TI",
 "V", "ZN", "ZR")

specnames <- mynames
lams <- lam[,1]

lamplotsnew <- function(lams, specnames)
{
  b <- barplot(lams,xaxt="n",xlab="",ylim=c(0,1))
  axis(1,at= b, labels=specnames, cex.axis=.7)
}

lamplotsnewlog <- function(lams, specnames)
{
  lams[lams < .0001] <- .0001
  b <- barplot(lams,xaxt="n",yaxt="n",xlab="",ylim=c(0.0001,1),log="y")
  axis(1,at= b, labels=specnames, cex.axis=.6)
  axis(2,at=c(.0001,.001,.01,.1,1),labels=c(".0001",".001",".01",".1","1"),las=2,cex=.3)
}

lamplotsnewlog(lam[,1],specnames)

#par(mfrow=c(5,1))
layout(matrix(1:10,5,2))
#layout.show(5)
for (i in 1:10) 
{
  par(mar=c(2,3,.5,1))
  lamplotsnew(lam[,i],specnames)
}

for (j in 1:10) 
{
  par(mar=c(2,3,.5,1))
  plot(as.Date(newdates,"%d%B%Y"),fac[,j],type="l")
}




lamplots2old(lam,"gkey",names=mynames) 
facplots3old(fac,"gkeyf",newdates)
facplots(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


lamplotslogps(lamORIG,"gkey",names=mynames,facnames[reorderfacs]) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)

%Latest result: 11/28/07
             mimic   LHT     NC
#1   lead    .011    .013    .011
#2  sulfate  .369    .33     .360
#3  copper   .012    .005    .016
#4  gas      .075    .16     .133
#5  zinc     .014    .013    .011
#6  C-r. sul .142    .20     .143
#7  soil     .032    .04     .032
#8  steel    .051    .07     .042
#9  diesel   .066    .02     .040
#10 nitrate  .228    .15     .213

facnames10 <- c("Lead\nSmelting","Secondary\nSulfate","Copper\nSmelting","Gasoline\nExhaust",
   "Zinc\nSmelting","Carbon-rich\nSulfate","Soil","Steel\nProcessing",
   "Diesel/\nrailroad","Secondary\nNitrate")

facnames9 <- c("Lead\nSmelter","Copper\nSmelter","Fireworks","Summer\nSecondary","Zinc\nSmelter",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecondary")
facnames9 <- c("Lead\nSmelt.","Copper\nSmelt.","Fire-\nworks","Summer\nSecond.","Zinc\nSmelt.",
               "Steel\nMill","Mobile","OC/SO4","Winter\nSecond.")
facnames9MORG  <- c("Mobile\n(high Cu)","Zinc\nSmelt.","Summer\nSecond.","Wood1","Wood2",
                    "Mobile1","SO/S/\nNa/EC","Steel\nMill","Winter\nSecond.")

quants <- seq(0,1,.025)
quants <- c(.01,.25,.5,.75,.99)
round(rbind(quants,
 exp( qnorm(quants, -.5*log(.25^2+1) , sqrt(log(.25^2+1)))),
 exp( qnorm(quants, -.5*log(.5^2+1) , sqrt(log(.5^2+1)))),
 exp( qnorm(quants, -.5*log(.75^2+1) , sqrt(log(.75^2+1))))  ), 2)



### New study evaluating perturbations of measurements (as opposed to uncertainties)

cvvals <- c(.01,.25,.50,.75)
reperrorvals <- 50

lamarrayerrsMM <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsMM <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

for (cvi in 1:length(cvvals))
{
  for (reperri in 1:reperrorvals)
  #for (reperri in 1:5)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    
    #if (cvi> 1) uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    #if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
    #                                      nrow(uncs),ncol(uncs)) #cvi=1
    if (cvi==1) concerrs <- concs * matrix( exp( rnorm(nrow(concs)*ncol(concs), -.5*log(.25^2+1) , sqrt(log(.25^2+1))) ) , nrow(concs),ncol(concs))
    if (cvi==2) concerrs <- concs * matrix( exp( rnorm(nrow(concs)*ncol(concs), -.5*log(.50^2+1) , sqrt(log(.50^2+1))) ) , nrow(concs),ncol(concs))
    if (cvi==3) concerrs <- concs * matrix( exp( rnorm(nrow(concs)*ncol(concs), -.5*log(.75^2+1) , sqrt(log(.75^2+1))) ) , nrow(concs),ncol(concs))
    if (cvi==4) concerrs <- concs * matrix( exp( rnorm(nrow(concs)*ncol(concs), -.5*log(.01^2+1) , sqrt(log(.01^2+1))) ) , nrow(concs),ncol(concs))
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(concerrs),file="concerrs.txt",ncol=ncol(concerrs))
    runpmfplain("default_nogkey.ini","concerrs.txt", "lhtunc2.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    


     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsMM[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsMM[,,cvi,reperri] <- as.matrix(fac)
     write(paste(cvi,"--",reperri), "repdat.txt", 1)
  }
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



### ALTERNATIVE...FASTER ###
cvi <- 1
faccompsall1MM <- matrix(NA,reperrorvals,numfacs)
lamcompsall1MM <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse1MM <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp1MM <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall1MM[i,] <- apply( cor(facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
  lamcompsall1MM[i,] <- apply( revcormat(lamarrayerrsMM[,,cvi,i],lamORIG) , 2, max)
  faccompsallrmse1MM[i,] <- apply( getaae(facarrayerrsMM[,,cvi,i],facORIG) , 2, min) ## now does aae
  lamcompsallexp1MM[i,] <- apply( expcormat(lamarrayerrsMM[,,cvi,i],lamORIG,
                                           facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1MM))
boxplot(data.frame(lamcompsall1MM))

date()

cvi <- 2
faccompsall2MM <- matrix(NA,reperrorvals,numfacs)
lamcompsall2MM <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse2MM <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp2MM <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall2MM[i,] <- apply( cor(facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
  lamcompsall2MM[i,] <- apply( revcormat(lamarrayerrsMM[,,cvi,i],lamORIG) , 2, max)
  faccompsallrmse2MM[i,] <- apply( getaae(facarrayerrsMM[,,cvi,i],facORIG) , 2, min) ## now does aae
  lamcompsallexp2MM[i,] <- apply( expcormat(lamarrayerrsMM[,,cvi,i],lamORIG,
                                           facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2MM))
boxplot(data.frame(lamcompsall2MM))


date()

cvi <- 3
faccompsall3MM <- matrix(NA,reperrorvals,numfacs)
lamcompsall3MM <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse3MM <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp3MM <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall3MM[i,] <- apply( cor(facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
  lamcompsall3MM[i,] <- apply( revcormat(lamarrayerrsMM[,,cvi,i],lamORIG) , 2, max)
  faccompsallrmse3MM[i,] <- apply( getaae(facarrayerrsMM[,,cvi,i],facORIG) , 2, min) ## now does aae
  lamcompsallexp3MM[i,] <- apply( expcormat(lamarrayerrsMM[,,cvi,i],lamORIG,
                                           facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3MM))
boxplot(data.frame(lamcompsall3MM))


cvi <- 4
faccompsall4MM <- matrix(NA,reperrorvals,numfacs)
lamcompsall4MM <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse4MM <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp4MM <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall4MM[i,] <- apply( cor(facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
  lamcompsall4MM[i,] <- apply( revcormat(lamarrayerrsMM[,,cvi,i],lamORIG) , 2, max)
  faccompsallrmse4MM[i,] <- apply( getaae(facarrayerrsMM[,,cvi,i],facORIG) , 2, min) ## now does aae
  lamcompsallexp4MM[i,] <- apply( expcormat(lamarrayerrsMM[,,cvi,i],lamORIG,
                                           facarrayerrsMM[,,cvi,i],facORIG) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall4MM))
boxplot(data.frame(lamcompsall4MM))







#par(mfcol=c(2,3))
#boxplot(data.frame(faccompsall2M),ylim=c(0,1),xlab="Contributions")
#title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
#boxplot(data.frame(lamcompsall2M),ylim=c(0,1),xlab="Profiles")

#boxplot(data.frame(faccompsall3M),ylim=c(0,1),xlab="Contributions")
#title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
#boxplot(data.frame(lamcompsall3M),ylim=c(0,1),xlab="Profiles")

#boxplot(data.frame(faccompsall1M),ylim=c(0,1),xlab="Contributions")
#title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
#boxplot(data.frame(lamcompsall1M),ylim=c(0,1),xlab="Profiles")


par(mfcol=c(2,3))
boxplot(data.frame(faccompsall1MM),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.75*true, 1.33*true)")
boxplot(data.frame(lamcompsall1MM),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall2MM),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 2*true)")
boxplot(data.frame(lamcompsall2MM),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall3MM),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.75*true, 1.33*true)")
boxplot(data.frame(lamcompsall3MM),ylim=c(0,1),xlab="Profiles")


pdf("MetalsUncertainty.pdf",width=12,height=9)
set.cex <- .65
par(mfcol=c(2,3))
maxy <- max(c( max(data.frame(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) ))  )) 
#boxplot(data.frame(faccompsallrmse2MM),ylim=c(0,2.5),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
title(paste("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp2MM),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp2MM)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse3),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse3M)/apply(facORIG,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
title(paste("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp3MM),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp3MM)) )) ,3)    ))

#boxplot(data.frame(faccompsallrmse1M),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse1M)/apply(facUNWT,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
#title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
title(paste("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)\n Mean Relative AAE =",
       round(mean(   as.matrix(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) )) ,3)    ))
boxplot(data.frame(lamcompsallexp1MM),ylim=c(0,1),xlab="Profiles",xaxt="n")
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("\n\nMean Revised Profile Correlation =",
       round(mean(   as.matrix(t(t(lamcompsallexp1MM)) )) ,3)    ))
dev.off()



newmeans1MM <- matrix(NA,reperrorvals,numfacs)
newmeans2MM <- matrix(NA,reperrorvals,numfacs)
newmeans3MM <- matrix(NA,reperrorvals,numfacs)
newmeans4MM <- matrix(NA,reperrorvals,numfacs)
for (reperri in 1:reperrorvals)
{
  matchbool <- t(t(getaae(facarrayerrsMM[,,1,reperri],facORIG))==faccompsallrmse1MM[reperri,])
  newmeans1MM[reperri,] <- apply(facarrayerrsMM[,,1,reperri],2,mean)[apply(matchbool*(1:10),2,sum)]
  matchbool <- t(t(getaae(facarrayerrsMM[,,2,reperri],facORIG))==faccompsallrmse2MM[reperri,])
  newmeans2MM[reperri,] <- apply(facarrayerrsMM[,,2,reperri],2,mean)[apply(matchbool*(1:10),2,sum)]
  matchbool <- t(t(getaae(facarrayerrsMM[,,3,reperri],facORIG))==faccompsallrmse3MM[reperri,])
  newmeans3MM[reperri,] <- apply(facarrayerrsMM[,,3,reperri],2,mean)[apply(matchbool*(1:10),2,sum)]
  matchbool <- t(t(getaae(facarrayerrsMM[,,4,reperri],facORIG))==faccompsallrmse4MM[reperri,])
  newmeans4MM[reperri,] <- apply(facarrayerrsMM[,,4,reperri],2,mean)[apply(matchbool*(1:10),2,sum)]
}

set.cex <- .85
junk <- apply(facarrayerrsMM[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)

#pdf("MetalsContributions4newM.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans4MM), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
title("(a1)")
points(1:10,apply(facORIG,2,mean),pch=4,cex=2)
#dev.off()

#pdf("MetalsContributions1newM.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans1MM), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
title("(a)")
points(1:10,apply(facORIG,2,mean),pch=4,cex=2)
#dev.off()

#pdf("MetalsContributions2newM.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans2MM), xaxt="n", ylim=c(0,ymaxcont),cex=set.cex)  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .05*true)"))
title("(b)")
points(1:10,apply(facORIG,2,mean),pch=4,cex=2)
#dev.off()

#pdf("MetalsContributions3newM.pdf",width=10,height=6)
par(mar=c(4,3,2,2))
boxplot(data.frame(newmeans3MM), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
#title(paste("Contributions\nUncertainties from N(true, .15*true)"))
title("(c)")
points(1:10,apply(facORIG,2,mean),pch=4,cex=2)
#dev.off()



maxy <- max(c( max(data.frame(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) )), 
               max(data.frame(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) ),na.rm=T), 
               max(data.frame(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) ))  )) 

pdf("MetalsRAAE2newM.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse2MM)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(b)")
dev.off()


pdf("MetalsRAAE3newM.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse3MM)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(c)")
dev.off()


pdf("MetalsRAAE1new.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) ),
      log="y",ylim=c(0.0025,maxy),xlab="",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%"),at=c(.005,.01,.05,.1,.5,1,5),las=2)
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(faccompsallrmse1MM)/apply(facORIG,2,mean)) )), col="gray", lwd=2)
title("(a)")
dev.off()

## Mean level for all is 10%, 18%, and 25% for 1,2,&3

## Checking the computational approach
scaledat <- matrix(NA,50,10)
for (ii in 1:50) scaledat[ii,] <- faccompsallrmse1MM[ii,] / apply(facORIG,2,mean)


allmin <- min(c(c(as.matrix(t(t(lamcompsallexp2MM)) )), 
                c(as.matrix(t(t(lamcompsallexp3MM)) )), 
                c(as.matrix(t(t(lamcompsallexp1MM)) ))))

pdf("MetalsEMC2M.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(lamcompsallexp2MM),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp2MM)) )), col="gray", lwd=2)
title("(b)")
dev.off()


pdf("MetalsEMC3M.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(lamcompsallexp3MM),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp3MM)) )), col="gray", lwd=2)
title("(c)")
dev.off()


pdf("MetalsEMC1M.pdf",width=10,height=6)
par(mar=c(4,4,2,2))
#set.cex <- 1.2
boxplot(data.frame(lamcompsallexp1MM),ylim=c(allmin,1),xlab="",xaxt="n")
par(mgp=c(3,2,0))
axis(1,labels=facnames10,at=1:10,cex.axis=set.cex)
abline(h= mean(   as.matrix(t(t(lamcompsallexp1MM)) )), col="gray", lwd=2)
title("(a)")
dev.off()

## Overall means: 99%, 97%, 90%
