lhtdatANC <- read.table("LHTconcALLNOCLEAN.txt",header=T,na.strings=".")
lhtuncANC <- read.table("LHTuncALLNOCLEAN.txt",header=T)

geomean <- function(x){ (prod(x,na.rm=T))^(1/length(na.omit(x))) }
geomean2 <- function(x){ (prod(x[x>0],na.rm=T))^(1/length(na.omit(x[x>0]))) }
geomean2 <- function(x){ 
  junk <- na.omit(x)
  junk <- junk[junk>0]
  exp(mean(log(junk)))
}
geomean2(c(1,.5,NA,.25))

junk <- as.matrix(lhtdatANC[,-1])
apply(junk,2,geomean2)
geomean(junk[,1])



newdates <- lhtdatANC[,1]
datesnew <- as.Date(newdates, "%d%b%Y")

lhtdat2ANC <- lhtdatANC[,-c(1,49)]  # remove Sulfur 
lhtunc2ANC <- lhtuncANC[,-c(1,49)]

uncs <- lhtunc2ANC
concs <- lhtdat2ANC



write(t(lhtdat2ANC), file="lhtdat2ANC.txt", ncol=ncol(lhtdat2ANC))
write(t(lhtunc2ANC), file="lhtunc2ANC.txt", ncol=ncol(lhtunc2ANC))

mynames <- c( "OC1" ,"OC2", "OC3", "OC4", "OP", "EC1", "EC2", "EC3", "SO4",  "NO3", 
 "NH4",  "AL", "AS", "BA", "CA", "CO", "CR", "CU", "FE", "HG",
 "K" ,"MN", "NI" ,"P", "PB", "RB" ,"SE", "SI", "SR", "TI",
 "V", "ZN", "ZR")
 mynames <- c( "OC1" ,"OC2", "OC3", "OC4", "OP", "EC1", "EC2", "EC3", "SO4",  "NO3", 
 "NH4",  "AL", "AS", "BA", "CA", "CO", "CR", "CU", "FE", "HG",
 "K" ,"MN", "NI" ,"P", "PB", "RB" ,"SE", "SI", "SR", "TI",
 "V", "ZN", "ZR",  
    "AG", "AU", "BR", "CD", "CL", "GA", "IN", "LA", "MG", "MO", "NA", "PD", "SB", "SN", #"SU", 
    "TL", "UR", "YT")
 

numfacs <- i <- 10
  concs <- lhtdat2ANC
  uncs <- lhtunc2ANC

#mynames <- dimnames(morgc)[[2]]
#mynames
names <- mynames
#newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
#newdates

runpmfplain("default_nogkey.ini","lhtdat2ANC.txt", "lhtunc2ANC.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

lamORIGANC<-read.table(mypaste("testlam",numfacs,".txt"))
facORIGANC<-read.table(mypaste("testfac",numfacs,".txt"))

lam <- lamORIGANC
fac <- facORIGANC

apply(facORIGANC,2,mean)
sum(apply(facORIGANC,2,mean))
 apply(facORIGANC,2,mean) / sum(apply(facORIGANC,2,mean))


%Latest result: 11/28/07
             mimic   LHT     NC     ANC
#1   lead    .011    .013    .011   .018
#2  sulfate  .369    .33     .360   .388  #3
#3  copper   .012    .005    .016   .024  #2
#4  gas      .075    .16     .133   .142
#5  zinc     .014    .013    .011   .012
#6  C-r. sul .142    .20     .143   .018  #8
#7  soil     .032    .04     .032   .036  #6
#8  steel    .051    .07     .042   .075  #7
#9  diesel   .066    .02     .040   .076  
#10 nitrate  .228    .15     .213   .211

specnames <- mynames

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





cvvals <- c(.01,.25,.5,.75)
reperrorvals <- 50

lamarrayerrsMM <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsMANC <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

for (cvi in 1:length(cvvals))
{
  for (reperri in 1:reperrorvals)
  #for (reperri in 15:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    #if (cvi> 1) uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    #if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
    #                                      nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    if (cvi==1) concerrs <- concs * matrix( exp( rnorm(nrow(uncs)*ncol(uncs), -.5*log(.25^2+1) , sqrt(log(.25^2+1))) ) , nrow(uncs),ncol(uncs))
    if (cvi==2) concerrs <- concs * matrix( exp( rnorm(nrow(uncs)*ncol(uncs), -.5*log(.50^2+1) , sqrt(log(.50^2+1))) ) , nrow(uncs),ncol(uncs))
    if (cvi==3) concerrs <- concs * matrix( exp( rnorm(nrow(uncs)*ncol(uncs), -.5*log(.75^2+1) , sqrt(log(.75^2+1))) ) , nrow(uncs),ncol(uncs))
    if (cvi==4) concerrs <- concs * matrix( exp( rnorm(nrow(uncs)*ncol(uncs), -.5*log(.01^2+1) , sqrt(log(.01^2+1))) ) , nrow(uncs),ncol(uncs))

    write(t(concerrs),file="concerrs.txt",ncol=ncol(uncerrs))
 
    runpmfplain("default_nogkey.ini","concerrs.txt", "lhtunc2ANC.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    


     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsMANC[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsMANC[,,cvi,reperri] <- as.matrix(fac)
     write(paste(cvi,"--",reperri), "repdat.txt", 1)
  }
}

## Note cvi of 1 is unif(.5,1.5), cvi of 2 is norm(1,.05), cvi of 2 is norm(1, .15)



### ALTERNATIVE...FASTER ###
cvi <- 1
faccompsall1MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsall1MANC <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse1MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp1MANC <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall1MANC[i,] <- apply( cor(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
  lamcompsall1MANC[i,] <- apply( revcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC) , 2, max)
  faccompsallrmse1MANC[i,] <- apply( getaae(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, min) ## now does aae
  lamcompsallexp1MANC[i,] <- apply( expcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC,
                                           facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1MANC))
boxplot(data.frame(lamcompsall1MANC))

date()

cvi <- 2
faccompsall2MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsall2MANC <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse2MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp2MANC <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall2MANC[i,] <- apply( cor(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
  lamcompsall2MANC[i,] <- apply( revcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC) , 2, max)
  faccompsallrmse2MANC[i,] <- apply( getaae(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, min) ## now does aae
  lamcompsallexp2MANC[i,] <- apply( expcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC,
                                           facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2MANC))
boxplot(data.frame(lamcompsall2MANC))


date()



cvi <- 3
faccompsall3MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsall3MANC <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse3MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp3MANC <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall3MANC[i,] <- apply( cor(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
  lamcompsall3MANC[i,] <- apply( revcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC) , 2, max)
  faccompsallrmse3MANC[i,] <- apply( getaae(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, min) ## now does aae
  lamcompsallexp3MANC[i,] <- apply( expcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC,
                                           facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3MANC))
boxplot(data.frame(lamcompsall3MANC))


date()



cvi <- 4
faccompsall4MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsall4MANC <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse4MANC <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp4MANC <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall4MANC[i,] <- apply( cor(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
  lamcompsall4MANC[i,] <- apply( revcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC) , 2, max)
  faccompsallrmse4MANC[i,] <- apply( getaae(facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, min) ## now does aae
  lamcompsallexp4MANC[i,] <- apply( expcormat(lamarrayerrsMANC[,,cvi,i],lamORIGANC,
                                           facarrayerrsMANC[,,cvi,i],facORIGANC) , 2, max)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall4MANC))
boxplot(data.frame(lamcompsall4MANC))


date()




par(mfcol=c(2,3))
boxplot(data.frame(faccompsall2MNC),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsall2MNC),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall3MNC),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsall3MNC),ylim=c(0,1),xlab="Profiles")

boxplot(data.frame(faccompsall1MNC),ylim=c(0,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsall1MNC),ylim=c(0,1),xlab="Profiles")

