junk <- matrix(rnorm(200),20,10)
junk2 <- rbind(junk,junk[20,]+rnorm(10,0,.01))
diag(junk %*% solve(t(junk) %*% junk) %*% t(junk))
diag(junk2 %*% solve(t(junk2) %*% junk2) %*% t(junk2))


sldat <- read.table("sldat.txt")
slunc <- read.table("slunc.txt")

sldat2 <- sldat
slunc2 <- slunc
sldat2[sldat== -999.9] <- NA
slunc2[sldat== -999.9] <- NA
sldat2 <- na.omit(sldat2)
slunc2 <- na.omit(slunc2)
datesnew <- getdates()

#alldat<-read.table("data_for_dates.txt",header=T)
#datesnew<-alldat[,1]
#datesnew<-as.Date(datesnew)

#namesALL <- names
#namesMET <- names[1:44]
names <- namesMET

newdates <- datesnew[apply(sldat==-999.9,1,sum)==0]
dimnames(sldat2)[[2]] <- namesMET
newdatesall <- getdates()[apply(sldat==-999.9,1,sum)==0]
#par(mfrow=c(1,1))
#plot(newdatesall,apply(sldat2,1,sum),type="l")
#junk <- apply(sldat2,1,sum)
#names(junk) <- newdatesall
#t(t(junk))  -- 4Jul2001 and 4&5July2002 are very large
#summary(sldat2)
#apply(sldat2,2,mean)
plot(newdatesall,sldat2$K,type="l")
junk <- sldat2$K
names(junk) <- newdatesall
t(t(junk))


Lambda <- lam
unc <- slunc2

getinf <- function(Lambda,unc,weight=1)
# Gets influence statistics (h_i) using Lambda, uncertainties, and weights
{
  h <- matrix(NA,nrow(unc),ncol(unc))
  for (i in 1:nrow(unc))
  {
    xw <- sqrt(diag(1/c(unc[i,]^2/weight))) %*% as.matrix(Lambda)
    h[i,] <- diag(xw %*% solve(t(xw) %*% xw) %*% t(xw))
  }
  h
}

wt <- rep(1,44)
wt[42] <- 100  # let the weight for EC be 100 times that of the other species
hmat <- getinf(lam,slunc2,weight=wt)
meanhmat <- round(apply(hmat,2,mean),5) 
names(meanhmat) <- names
meanhmat * 1000


Lambda <- lam
#Lambda <- lam[1:44,]

if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

wtadj <- rep(1,ncol(sldat2))
#wtadj <- c(rep(100,10),rep(1,34))
wt <- c(unc[i,]^2/wtadj)

obs <- 200
#lm.sl <- lm(as.double(sldat2[obs,])[1:43] ~ l1 + l2 + l3 + l4 + l5 + l6 - 1)#, weights=wt)
lm.sl <- lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, weights=wt)
inf.sl <- influence.measures(lm.sl)
summary(inf.sl)
#dffits(lm.sl)
#rstudent(lm.sl)
sumdfbeta <- apply(abs(inf.sl[[1]][,1:ncol(Lambda)]),1,sum)
cor(cbind(inf.sl[[1]][,7:10],rstudent(lm.sl),sumdfbeta))
pairs(cbind(inf.sl[[1]][,7:10],rstudent(lm.sl),sumdfbeta))

 
dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, weights=wt))

apply(abs(dffitvals),2,median)

dfplot<-function(dffitvals,namespec=names){
	par(cex=1)
      temp <- (apply(abs(dffitvals),2,median))
	plot(temp,ylim=c(0,max(temp)+.05),ylab="Median of Abs(dffits)",
	xlab="Species",xaxt="n",type="n")
	par(cex=.6)
	axis(1,at=1:ncol(dffitvals),label=namespec,las=2)
	###  makes a box around significant values
	par(cex=1)	
	for(ii in 1:length(temp)){
	  if(temp[ii]>=quantile(temp,.9,na.rm=TRUE)){
	  rect(ii-.5,par("usr")[3],ii+.5,par("usr")[4]-.01,col="grey90",border="transparent")
	}}
	points(temp,bty="o") #redraws the points on top of the shaded boxes
	#end of box maxing
}
par(mfrow=c(1,1))
dfplot( dffitvals )

dfplot2<- function(dffitvals,namespec=names,maxval=log10(max(abs(dffitvals))+1)){
   par(mfrow=c(3,4))
   for (i in 1:ncol(dffitvals)) 
   {
      plot(newdates,log10(abs(dffitvals[,i])+1),type="l",ylim=c(0,maxval))
      #plot(1:nrow(dffitvals),dffitvals[,i],ylim=c(-maxval,maxval))
      title(namespec[i])
   }
}

maxval <- log10(max(abs(dffitvals))+1)
dfplot2( dffitvals[,1:12], namespec=names[1:12],maxval )
dfplot2( dffitvals[,13:24], namespec=names[13:24],maxval )
dfplot2( dffitvals[,25:36], namespec=names[25:36],maxval )
dfplot2( dffitvals[,37:44], namespec=names[37:44],maxval )


   
## TRY TO ADJUST WEIGHTS SO THAT ALL SPECIES HAVE EQUAL INFLUENCE (d_i) DEFINED BY temp
###

wtadj <- rep(1,ncol(sldat2))
dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, 
          weights= 1/ (c(unc[obs,]^2/wtadj)) ))

dfplot(dffitvals)
temp <-  apply(abs(dffitvals),2,median)



## RUN THE CODE BELOW MULTIPLE TIMES TO FINE TUNE WEIGHT ADJUSTMENTS
for (jj in 1:10)
{
  temp <-  apply(abs(dffitvals),2,median)
  wtadj <- wtadj * (temp[44]/temp[1:44])
  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, 
     weights= 1/c(unc[obs,]^2/wtadj) ))
  par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration",jj))
}
## END OF WEIGHT ADJUSTMENT CODE



wtadjfin <- wtadj
rbind(names,round(wtadjfin))

maxval <- log10(max(abs(dffitvals))+1)
dfplot2( dffitvals[,1:12], namespec=names[1:12],maxval )
dfplot2( dffitvals[,13:24], namespec=names[13:24],maxval )
dfplot2( dffitvals[,25:36], namespec=names[25:36],maxval )
dfplot2( dffitvals[,37:44], namespec=names[37:44],maxval )





## ANOTHER EXAMPLE: ADJUST WEIGHTS SO THAT SPECIES 41:44 OF SUM(d_i) THAT IS K=5 TIMES LARGER
## THAN SUM(d_i) FOR SPECIES 1:40.
K <- 5

wtadj <- rep(1,ncol(sldat2))
dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, 
          weights= 1/ (c(unc[obs,]^2/wtadj)) ))

dfplot(dffitvals)
temp <-  apply(abs(dffitvals),2,median)
dold <- temp

sum(temp[1:40])
sum(temp[41:44])

## RUN THE CODE BELOW MULTIPLE TIMES TO FINE TUNE WEIGHT ADJUSTMENTS
for (jj in 1:5)
{
  jj
  temp <-  apply(abs(dffitvals),2,median)
  wtadj[41:44] <- wtadj[41:44] * (sum(temp[1:40])/sum(temp[41:44])*K)^2
      # tried the squared adjustment to speed up
  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 - 1, 
  weights= 1/c(unc[obs,]^2/wtadj)))
  par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration",jj,"-- ratio is",sum(temp[41:44])/sum(temp[1:40]),
              "(should be K=",K,")"))
}
## END OF WEIGHT ADJUSTMENT CODE




wtadjfin <- wtadj
rbind(names,round(wtadjfin))

maxval <- log10(max(abs(dffitvals))+1)
dfplot2( dffitvals[,1:12], namespec=names[1:12],maxval )
dfplot2( dffitvals[,13:24], namespec=names[13:24],maxval )
dfplot2( dffitvals[,25:36], namespec=names[25:36],maxval )
dfplot2( dffitvals[,37:44], namespec=names[37:44],maxval )





###
###
###
### ADAPT APPROACH TO USE WITH PMF
###
###
###

## Equal influence for all species

i <- 3
  write(t(slunc2),file="slunc2.txt",ncol=ncol(slunc2))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
runpmfplain("default_nogkey.ini","sldat2.txt","slunc2.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))
lamUNWT <- lam
facUNWT <- fac
lamplots(lam,"gkey") 
facplots2(fac,"gkeyf",newdates)
weekendplots2(i,"gkeyw",newdates)
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
## not working: system(mypaste("copy"," gkeyn6ADJ.jpg"," gkeyn6UNWT.jpg"))

Lambda <- lam
if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(slunc2[obs,]^2) ))

par(mfrow=c(4,5))
dfplot(dffitvals)

tempold <- apply(abs(dffitvals),2,median)

wtadj <- rep(1,ncol(sldat2))
for (jjj in 1:5)
{
  temp <-  apply(abs(dffitvals),2,median)
  wtadj <- wtadj * (temp[44]/temp[1:44])
  #wtadj <- wtadj * (1/temp[1:44])
  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
     weights= 1/c(unc[obs,]^2/wtadj) ))
  par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration 0: fine-tune \#:",jjj))
}
temp <-  apply(abs(dffitvals),2,median)
oldtemp <- temp
#newuncall <- newunc
#newunc <- unc
tempunc <- matrix(NA,nrow(newunc),ncol(newunc))
for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( unc[obs,]^2 / wtadj )
newunc <- tempunc
newunc.ready <- TRUE

## RUN THE CODE BELOW MULTIPLE TIMES TO FINE TUNE WEIGHT ADJUSTMENTS
for (jj in 1:3)
{
  #if (newunc.ready==FALSE)
  #{
  #  tempunc <- matrix(NA,nrow(newunc),ncol(newunc))
  #  #for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( newunc[obs,]^2 / sqrt(1/temp[1:44]) )
  #  for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( newunc[obs,]^2 / sqrt(temp[44]/temp[1:44]) )
  #  newunc <- tempunc
  #}
  if (min(newunc) < 0.00001) newunc <- newunc*0.00001/min(newunc)
  write(t(newunc),file="newunc.txt",ncol=ncol(newunc))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
  runpmfplain("default_nogkey.ini","sldat2.txt", "newunc.txt" ,mypaste("testfac",i,".txt"),
     mypaste("testlam",i,".txt"),
     mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
  lam<-read.table(mypaste("testlam",i,".txt"))
  Lambda <- lam
  if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
  if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
  if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
  if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
  if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
  if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
  if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
  if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
  if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(newunc[obs,]^2) ))
  #par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration",jj))

  wtadj <- rep(1,ncol(sldat2))
  for (jjj in 1:5)
  {
    temp <-  apply(abs(dffitvals),2,median)
    #wtadj <- wtadj * (1/temp[1:44])
    wtadj <- wtadj * (temp[44]/temp[1:44])
    dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
    for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
       dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
       weights= 1/c(newunc[obs,]^2/wtadj) ))
    #par(mfrow=c(1,1))
    dfplot(dffitvals)
    title(paste("iteration",jj,": fine-tune \#:",jjj)) 
  }

  temp <-  apply(abs(dffitvals),2,median)
  tempunc <- matrix(NA,nrow(newunc),ncol(newunc))
  for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( newunc[obs,]^2 / wtadj )
  newunc <- tempunc
  #newunc.ready <- TRUE
  lam<-read.table(mypaste("testlam",i,".txt"))
  fac<-read.table(mypaste("testfac",i,".txt"))
  lamjj <- lam
  facjj <- fac
  lamplots(lamjj,"gkey") 
  facplots2(facjj,"gkeyf",newdates)
  weekendplots2(i,"gkeyw",newdates)
  system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
}
## END OF WEIGHT ADJUSTMENT CODE

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 3/8/2006 @@@@@@@@@@@@@@@@@@@@@@@@@@@@


###
###
### Weight the last four (41:44) is K times the weight of the first 40 species 
###
###

K <- 1
i <- 7
  write(t(slunc2),file="slunc2.txt",ncol=ncol(slunc2))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
runpmfplain("default_nogkey.ini","sldat2.txt","slunc2.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))
lamUNWT <- lam
facUNWT <- fac
lamplots(lam,"gkey") 
facplots2(fac,"gkeyf",newdates)
weekendplots2(i,"gkeyw",newdates)
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
system(mypaste("multipanel_gkey2COL.bat"),invisible=TRUE) 
    #see License in this directory for the imagemagick license
## not working: system(mypaste("copy"," gkeyn6ADJ.jpg"," gkeyn6UNWT.jpg"))

Lambda <- lam
if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(slunc2[obs,]^2) ))

tempold <- apply(abs(dffitvals),2,median)

par(mfrow=c(3,4))
dfplot(dffitvals)
  title(paste("iteration 0-- ratio is",round(sum(tempold[41:44])/sum(tempold[1:40]),3),
              "(K=",K,")"))

wtadj <- rep(1,ncol(sldat2))
for (jjj in 1:5)
{
  temp <-  apply(abs(dffitvals),2,median)
  wtadj[41:44] <- wtadj[41:44] * (sum(temp[1:40])/sum(temp[41:44]) *K)
  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
     weights= 1/c(slunc2[obs,]^2/wtadj) ))
  #par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration 0: fine-tune \#:",jjj,"-- \n ratio is",round(sum(temp[41:44])/sum(temp[1:40]),3),
              "(K=",K,")"))
}
temp <-  apply(abs(dffitvals),2,median)
oldtemp <- temp
tempunc <- matrix(NA,nrow(newunc),ncol(newunc))
for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( slunc2[obs,]^2 / wtadj )
newunc <- tempunc

## RUN THE CODE BELOW MULTIPLE TIMES TO FINE TUNE WEIGHT ADJUSTMENTS
for (jj in 5:6)
{
  if (min(newunc) < 0.00001) newunc <- newunc*0.00001/min(newunc)
  write(t(newunc),file="newunc.txt",ncol=ncol(newunc))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
  runpmfplain("default_nogkey.ini","sldat2.txt", "newunc.txt" ,mypaste("testfac",i,".txt"),
     mypaste("testlam",i,".txt"),
     mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
  lam<-read.table(mypaste("testlam",i,".txt"))
  Lambda <- lam
  if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
  if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
  if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
  if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
  if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
  if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
  if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
  if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
  if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

  dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
  for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
     dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(newunc[obs,]^2) ))
  #par(mfrow=c(1,1))
  dfplot(dffitvals)
  title(paste("iteration",jj,"-- \n ratio is",round(sum(temp[41:44])/sum(temp[1:40]),3),
              "(K=",K,")"))

  wtadj <- rep(1,ncol(sldat2))
  for (jjj in 1:5)
  {
    temp <-  apply(abs(dffitvals),2,median)
    wtadj[41:44] <- wtadj[41:44] * (sum(temp[1:40])/sum(temp[41:44]) *K)
    dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
    for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
       dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
       weights= 1/c(newunc[obs,]^2/wtadj) ))
    #par(mfrow=c(1,1))
    dfplot(dffitvals)
    title(paste("iteration",jj,": fine-tune \#:",jjj,"-- \n ratio is",round(sum(temp[41:44])/sum(temp[1:40]),3),
              "(K=",K,")"))
  }

  temp <-  apply(abs(dffitvals),2,median)
  tempunc <- matrix(NA,nrow(newunc),ncol(newunc))
  for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( newunc[obs,]^2 / wtadj )
  newunc <- tempunc
  lam<-read.table(mypaste("testlam",i,".txt"))
  fac<-read.table(mypaste("testfac",i,".txt"))
  lamjj <- lam
  facjj <- fac
  lamplots(lamjj,"gkey") 
  facplots2(facjj,"gkeyf",newdates)
  weekendplots2(i,"gkeyw",newdates)
  system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
}
## END OF WEIGHT ADJUSTMENT CODE

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

## In the example above with K=1, sources 1, 5, 6, 7 remained almost the same 
## with sources 2, 3, 4, changing dramatically



## JUNK HERE
plot(newdates,sldat2[,8],type="l")
axis(1,labels=c(as.Date(newdates[sldat2[,8] > quantile(sldat2[,8],.998)])),
  at=newdates[sldat2[,8] > quantile(sldat2[,8],.998)])
#  at=(1:length(newdates))[sldat2[,8] > quantile(sldat2[,8],.998)])










##
## NEXT UP...PUT TOGETHER DATA SET AND DO THE WEIGHT CHANGING ANALYSIS
##

metals <- cbind(newdates,sldat2,slunc2)
dimnames(metals)

write(t(metals),file="../receptor/StLouis/metals.txt",ncol=ncol(metals))
write(dimnames(metals)[[2]],file="../receptor/StLouis/metalsnames.txt")


##
## RUN makesmalldata.sas HERE
##

morgc <- read.table("../receptor/StLouis/metorgconc.txt",header=T)
morgc[1,]

morgu <- read.table("../receptor/StLouis/metorgunc.txt",header=T)
morgu[1,]

apply(morgc/morgu,2,mean)[order(apply(morgc/morgu,2,mean),decreasing=T)]
           S           Fe           Zn           Ca           Cu            K           SO levoglucosan           NO 
 170.6435927   92.3724838   59.4710936   49.4597856   37.3135322   31.1763945   26.9999999   21.1673179   21.0000000 
          Si           Br           Pb           OC      hopste3         pah3           Mn         pah4      hopste2 
  20.3247754   19.7120566   17.0797703    9.7699059    9.6726184    8.9171688    8.2905716    8.0910245    7.6229539 
        pah6     hopste15           Cl         pah8         pah1     hopste12      hopste1     hopste10     hopste11 
   7.3956107    6.6477327    6.3610001    6.1564463    6.0624211    5.9074266    5.6551534    5.5896064    5.5207107 
     diacid2      hopste5         pah5      hopste8     hopste14      hopste4     hopste13     hopste16      hopste6 
   5.5063097    5.4730373    5.4602777    5.3458110    5.3445311    5.2831736    5.1904542    5.0760199    5.0622803 
     diacid5      diacid7      diacid6      hopste7      hopste9         pah7      diacid1     diacid12      diacid9 
   4.9692095    4.9681717    4.9607002    4.9590303    4.9299417    4.9117023    4.9063253    4.8835649    4.7308005 
        pah9     diacid10      diacid4         pah2           Se     diacid11           Al      diacid3           EC 
   4.7241379    4.7102456    4.7099480    4.6625772    4.4194309    4.3972082    4.3938818    3.4072562    3.3522497 
     diacid8           Mg           Ph           Na           Sr           Hg           Ba           As           Zr 
   3.3064518    2.9044524    2.3099402    2.2243234    2.2034354    1.8036320    1.6203091    1.0022440    0.6990304 
          Cd           Sn           Sb           Rb           Ni           Co           Ag           Mo           La 
   0.6681382    0.6319733    0.5932903    0.4986152    0.4598682    0.3415297    0.3050448    0.2658790    0.2595371 
          Tl           Pd           Cr            Y            V            U           Ga           In           Ti 
   0.2447029    0.2420411    0.2078984    0.1971828    0.1930859    0.1488386    0.1400491    0.1369765    0.1179009 
          Au 
   0.0642095 

apply(morgc/morgu,2,median)[order(apply(morgc/morgu,2,median),decreasing=T)]
           S           Fe           Zn           Ca            K           SO           Cu           NO           Si 
152.16721512  92.14067194  49.33381903  46.39104172  27.77953836  27.00000000  23.18395752  20.99999998  19.16032949 
levoglucosan           Br           Pb           OC         pah1           Mn      hopste3      hopste1      hopste4 
 16.48493987  15.37183402  13.62616498   9.55974395   6.60776723   6.53412841   5.01030998   5.00000000   5.00000000 
     hopste5      hopste6      hopste7      hopste8      hopste9     hopste10     hopste11     hopste12     hopste13 
  5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000 
    hopste14     hopste15     hopste16         pah4         pah5         pah6         pah7         pah8         pah9 
  5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000   5.00000000 
     hopste2      diacid6      diacid5      diacid1      diacid7     diacid12         pah2         pah3      diacid2 
  5.00000000   4.99186570   4.97486422   4.95382554   4.93391557   4.93073525   4.91245804   4.87275192   4.86060855 
    diacid10      diacid9      diacid4     diacid11           Al           Se      diacid3      diacid8           Mg 
  4.82300538   4.80835641   4.79201101   4.65543809   4.44523242   4.19617365   3.59325771   3.49455218   3.47616444 
          EC           Sr           Na           As           Zr           Ba           Rb           Ph           Co 
  2.96296296   0.83422594   0.82707829   0.48125697   0.43251542   0.43057616   0.36251669   0.32150478   0.30151533 
          Tl           Sn           La           Cd           Ag            Y           Mo           Pd            U 
  0.20009853   0.18558469   0.18171079   0.15579874   0.13885607   0.11297337   0.09783403   0.09406202   0.06858147 
          Sb           Ni           Cl           Ti            V           Cr           Ga           In           Au 
  0.04112810   0.01777814   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000 
          Hg 
  0.00000000 

eigen(cor(morgc))$values
 [1] 20.655175693  7.007834257  4.796415092  4.288046553  3.412492896  2.862458680  2.790727851  2.532384129  2.162560953
[10]  1.965205574  1.708056063  1.630111632  1.528992299  1.356194672  1.325995172  1.292210200  1.192688470  1.177312946
[19]  1.170465674  1.089461211  1.053951462  0.952079662  0.872096604  0.850648430  0.818864591  0.752327437  0.715170471
[28]  0.643897159  0.625879365  0.615555839  0.561357873  0.550999704  0.515264800  0.479624116  0.462106543  0.415649682
[37]  0.394326699  0.380585063  0.346597145  0.323879517  0.312697679  0.259378216  0.238284226  0.230008236  0.225910624
[46]  0.201418544  0.176052825  0.174728615  0.162179333  0.152833515  0.147509210  0.133348732  0.119613186  0.116024782
[55]  0.109790517  0.105077798  0.099220370  0.083056307  0.075314720  0.069513376  0.065465725  0.060501264  0.052080353
[64]  0.050397340  0.045302442  0.034892055  0.031221151  0.030516130  0.025812937  0.022680880  0.018160514  0.016672251
[73]  0.015267506  0.013738576  0.010729652  0.009848249  0.009579311  0.006473771  0.006008769  0.003241273  0.002711502
[82]  0.001085364

eigen(cor(morgc/morgu))$values

eigen(cor(morgc[,1:44]))$values

avecorr <- function(x,grps)
{
  grptab <- table(grps)
  k <- length(grptab)
  vals <- matrix(NA, length(grptab), length(grptab))
  for (i in 1:k) vals[i,i] <- 
     (sum( x[grps==names(grptab)[i],grps==names(grptab)[i]] ) - grptab[i]) / (grptab[i]^2 - grptab[i])
  for (i in 2:k) for (j in 1:(i-1)) vals[i,j] <- mean(x[grps==names(grptab)[i],grps==names(grptab)[j]])
  dimnames(vals) <- list(names(grptab),names(grptab))
  vals
}

avecorr(cor(morgc[,45:82]),c(rep("1.diacids",12),rep("2.hop-stear",16),rep("3.pahs",9),"4.levogluc"))
            1.diacids 2.hop-stear    3.pahs 4.levogluc
1.diacids   0.5511768          NA        NA         NA
2.hop-stear 0.2359935   0.6643681        NA         NA
3.pahs      0.2462205   0.4460423 0.7373831         NA
4.levogluc  0.3558998   0.6035459 0.6116094        NaN





























































###
###
### Weight organics against metals
###
###

numfacs <- i <- 7

  write(t(morgc),file="morgc.txt",ncol=ncol(morgc))
  write(t(morgu),file="morgu.txt",ncol=ncol(morgu))
  concs <- morgc
  uncs <- morgu
mynames <- dimnames(morgc)[[2]]
mynames
names <- mynames
newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
newdates
# datevec <- newdates
runpmfplain("default_nogkey.ini","morgc.txt","morgu.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0.2,i,nrow(morgc),ncol(morgc))    

lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))
lamUNWT <- lam
facUNWT <- fac
#lamUNWTold <- lamUNWT
#facUNWTold <- facUNWT
lamplots2old(lamUNWT,"morg",names=mynames) 
#lamplotslog(lamUNWT,"morg",names=mynames) 
facplots2(facUNWT,"morgf",newdates)
#weekendplots2(i,"morgw",newdates)
weekendplots3(i,newdates,facname=facUNWT,"morgw")
system(mypaste("multipanel_morg",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
## not working: system(mypaste("copy"," gkeyn6ADJ.jpg"," gkeyn6UNWT.jpg"))

## Try with just metalsplus
  write(t(morgc[,1:44]),file="morgmetc.txt",ncol=ncol(morgc[,1:44]))
  write(t(morgu[,1:44]),file="morgmetu.txt",ncol=ncol(morgu[,1:44]))
mynames <- dimnames(morgc[,1:44])[[2]]
runpmfplain("default_nogkeyMET.ini","morgmetc.txt","morgmetu.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(morgc[,1:44]),ncol(morgc[,1:44]))   
      # had to adjust lims to 10,.3,.1 
lamtemp<-read.table(mypaste("testlam",i,".txt"))
factemp<-read.table(mypaste("testfac",i,".txt"))
#lammet <- lamtemp
#facmet <- factemp
lamplots2old(lammet,"morg",names=mynames) 
#lamplotslog(lammet,"morg",names=mynames) 
facplots2(facmet,"morgf",newdates)
weekendplots3(i,newdates,facname=facmet,"morgw")
system(mypaste("multipanel_morg",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license

cor(facmet,facUNWT)

#AAER <- function(x,true){
#  if(dim(x)[1]<dim(x)[2])
#  {
#    x<-t(x);true<-t(true)
#  }
#  AAEval <- sum(apply(sqrt((x - true)^2)/true, 2, sum)/dim(x)[1])
#  return(AAEval)
#}
#AAE(lamUNWT[1:44,],lammet)
#AAER(lammet,lamUNWT[1:44,])
#AAERmat <- function(x,true)
#{
#  vals <- matrix(NA,ncol(x),ncol(true))
#  for (i in 1:ncol(x)) for (j in 1:ncol(true)) vals[i,j] <-  mean(abs(x[,i]-true[,j]) /true[,j])
#  vals
#}
#AAERmat(lammet,lamUNWT[1:44,])

stdize <- function(x)
{
  means <- apply(x,2,mean)
  stds <- sqrt(apply(x,2,var))
  vals <- matrix(NA,nrow(x),ncol(x))
  for (i in 1:nrow(x)) vals[i,] <- t((x[i,]-means)/stds)
  vals
}

revcormat <- function(x,trueval)
{
  temp <- cbind(x,trueval)
  temp <- t(stdize(t(temp)))
  cor(temp[,(1:ncol(x))],temp[,-(1:ncol(x))])
}

expcormat <- function(x,y,xfac,yfac)
{
  xbar <- apply(xfac,2,mean)
  ybar <- apply(yfac,2,mean)
  tempx <- t(t(x)*xbar)
  tempy <- t(t(y)*ybar)
  xexplained <- matrix(NA,nrow(tempx),ncol(tempx))
  yexplained <- matrix(NA,nrow(tempy),ncol(tempy))
  for (spec in 1:nrow(tempx)) 
  {
    xexplained[spec,] <- tempx[spec,]/sum(tempx[spec,])
    yexplained[spec,] <- tempy[spec,]/sum(tempy[spec,])
  }
  #return(cor(xexplained,yexplained),cbind(round(xexplained,3),round(yexplained,3)))
  cor(xexplained,yexplained)
}

expcormat(lammet,lamUNWT[-grp2,],facmet,facUNWT)


getrmse <- function(x,y)
{
  temp <- matrix(NA,ncol(x),ncol(y))
  for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- sqrt(mean( (x[,i] - y[,j])^2  ))
  #for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- sqrt(mean( (x[,i] - y[,j])^2/
                                                               #(var(x[,i])*var(y[,j]))  ))
  temp
}

getmse <- function(x,y)
{
  temp <- matrix(NA,ncol(x),ncol(y))
  for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- (mean( (x[,i] - y[,j])^2 ))
  temp
}

getmrae <- function(x,y)
{
  temp <- matrix(NA,ncol(x),ncol(y))
  for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- median( abs(x[,i] - y[,j])/((x[,i]+y[,j])/2) )
  temp
}

getaae <- function(x,y)
{
  temp <- matrix(NA,ncol(x),ncol(y))
  for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- mean( abs(x[,i] - y[,j]) )
  temp
}

getmccor <- function(x,y)
{
  temp <- matrix(NA,ncol(x),ncol(y))
  for (i in 1:ncol(x)) for (j in 1:ncol(y)) temp[i,j] <- cor(x[,i],y[,j])*
                                   min(mean(x[,i]),mean(y[,j]))/max(mean(x[,i]),mean(y[,j]))
  temp
}


getmse(facmet,facUNWT)
getrmse(facmet,facUNWT)
getmrae(facmet,facUNWT)

round(cor(facmet,facUNWT),3)
round(getmccor(facmet,facUNWT),3)

revcormat(lammet,lamUNWT[1:44,])   ## Works well

source("C:/Program Files/Insightful/splus61/library/combinat/COMBINAT.Q")
bestcormat <- function(x,y)
{
  mat <- cor(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  max(meancor)
}

bestcormatcors <- function(x,y)
{
  mat <- cor(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  #max(meancor)
  diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==max(meancor)]]] ])
}

bestrevcormat <- function(x,y)
{
  mat <- revcormat(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  max(meancor)
}

bestrevcormatcors <- function(x,y)
{
  mat <- revcormat(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  #max(meancor)
  diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==max(meancor)]]] ])
}

bestexpcormat <- function(x,y,facx,facy)
{
  mat <- expcormat(x,y,facx,facy)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  max(meancor)
}
bestexpcormat(lammet,lamUNWT[-grp2,],facmet,facUNWT)

bestexpcormatcors <- function(x,y,facx,facy)
{
  mat <- expcormat(x,y,facx,facy)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- mean(diag(mat[1:nrow(mat),perms[[i]]]))
  #max(meancor)
  diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==max(meancor)]]] ])
}

bestmsemat <- function(x,y)
{
  mat <- getmse(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- sum(diag(mat[1:nrow(mat),perms[[i]]]))
  return(min(meancor), mat, perms[[(1:length(perms))[meancor==min(meancor)]]],
         diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ]) )
}

bestaaemat <- function(x,y)
{
  mat <- getaae(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- sum(diag(mat[1:nrow(mat),perms[[i]]]))
  #return(min(meancor), mat, perms[[(1:length(perms))[meancor==min(meancor)]]],
  #       diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ]) )
  min(meancor)
}

bestrmsemat <- function(x,y)
{
  mat <- getrmse(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- sum(diag(mat[1:nrow(mat),perms[[i]]]))
  #return(min(meancor), mat, perms[[(1:length(perms))[meancor==min(meancor)]]],
  #       diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ]) )
  min(meancor)
}

bestrmsemat(facmet,facUNWT)

bestrmsematcors <- function(x,y)
{
  #mat <- getrmse(x,y)
  mat <- getaae(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- sum(diag(mat[1:nrow(mat),perms[[i]]]))
  #return(min(meancor), mat, perms[[(1:length(perms))[meancor==min(meancor)]]],
  #       diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ]) )
  #min(meancor)
  diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ])
}

bestaaematmeans <- function(x,y)
{
  #mat <- getrmse(x,y)
  mat <- getaae(x,y)
  perms <- permn(1:ncol(mat))
  meancor <- rep(NA,length(perms))
  for (i in 1:length(perms)) meancor[i] <- sum(diag(mat[1:nrow(mat),perms[[i]]]))
  #return(min(meancor), mat, perms[[(1:length(perms))[meancor==min(meancor)]]],
  #       diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ]) )
  #min(meancor)
  diag(mat[1:nrow(mat),perms[[(1:length(perms))[meancor==min(meancor)]]] ])
}




## Try with just organics plus Al, Si, OC, EC (#3,4,41,42)
i <- 7
  write(t(morgc[,c(3,4,41,42,45:82)]),file="morgorgc.txt",ncol=ncol(morgc[,c(3,4,41,42,45:82)]))
  write(t(morgu[,c(3,4,41,42,45:82)]),file="morgorgu.txt",ncol=ncol(morgu[,c(3,4,41,42,45:82)]))
#  concs <- morgc[,c(3,4,41,42,45:82)]
#  uncs <- morgu[,c(3,4,41,42,45:82)]
mynames <- dimnames(morgc[,c(3,4,41,42,45:82)])[[2]]
runpmfplain("default_nogkey.ini","morgorgc.txt","morgorgu.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(morgc[,c(3,4,41,42,45:82)]),ncol(morgc[,c(3,4,41,42,45:82)]))    
lamtemp<-read.table(mypaste("testlam",i,".txt"))
factemp<-read.table(mypaste("testfac",i,".txt"))
#cor(factemp,facorg)
#revcormat(lamtemp,lamorg)
#lamorg <- lamtemp
#facorg <- factemp
lamplots2old(lamorg,"morg",names=mynames) 
#lamplotslog(lamorg,"morg",names=mynames) 
facplots2(facorg,"morgf",newdates)
weekendplots3(i,newdates,facname=facorg,"morgw")
system(mypaste("multipanel_morg",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license

cor(facorg,facUNWT)

revcormat(lamorg,lamUNWT[c(3,4,41,42,45:82),])   ## Works well
#cor(lamorg,lamUNWT[c(3,4,41,42,45:82),])   ## Works poorly


lamplots(lams7,"gkey")







K <- .5
totreps <- 15
finetunereps <- 6
numfacs <- 7

Lambda <- lamUNWT
if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

dffitvals <- matrix(NA,nrow(concs),ncol(concs))
for (obs in 1:nrow(concs)) dffitvals[obs,] <- 
   dffits(lm(as.double(concs[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(uncs[obs,]^2) ))

tempold <- apply(abs(dffitvals),2,median)
temp <-  apply(abs(dffitvals),2,median)
dfplot(dffitvals)




#grp1 <- 1:44
#grp2 <- c(3,4,41,42,45:82)
grp1 <- c(1:2,5:40,43:44)
grp2 <- c(45:82)
# K uses grp2/grp1

par(mfrow=c(1,2))
dfplot(dffitvals)
Ktemp <- matrix(NA,totreps+1,finetunereps+1)
Ktemp[1,1] <- Knat <- round(sum(tempold[grp2])/sum(tempold[grp1]),3)
  title(paste("iteration 0-- ratio is",signif(sum(tempold[grp2])/sum(tempold[grp1]),3),
              "(K=",K,")"))
lamarray <- array(NA,dim=c(ncol(concs),numfacs,totreps+1))
facarray <- array(NA,dim=c(nrow(concs),numfacs,totreps+1))

lamarray[,,1] <- as.matrix(lamUNWT)
facarray[,,1] <- as.matrix(facUNWT)

uniquify <- function(x,grp1,grp2)
{
  allvals <- 1:(length(x))
  G1 <- G2 <- rep(F,length(x))
  for (spec in 1:length(grp1)) G1[grp1[spec]] <- T
  for (spec in 1:length(grp2)) G2[grp2[spec]] <- T
  g1 <- x[G1 & !G2] 
  g2 <- x[G2 & !G1]
  s <- x[G1 & G2]
  return(g1,g2,s)
}

wtadj <- rep(1,ncol(concs))
for (jjj in 1:finetunereps)
{
  groups <- uniquify(1:(ncol(concs)),grp1,grp2)

  if (K < Knat) wtadj[grp2] <- wtadj[grp2] * (K*sum(temp[groups$g1])) / sum(temp[groups$g2])
  if (K > Knat) wtadj[grp1] <- wtadj[grp1] * sum(temp[groups$g2]) / (K*sum(temp[groups$g1]))
  dffitvals <- matrix(NA,nrow(concs),ncol(concs))
  for (obs in 1:nrow(concs)) dffitvals[obs,] <- 
     dffits(lm(as.double(concs[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
     weights= 1/c(uncs[obs,]^2/wtadj) ))
  #par(mfrow=c(1,1))
#  dfplot(dffitvals)
  temp <-  apply(abs(dffitvals),2,median)
#  title(paste("iteration 0: fine-tune \#:",jjj,"-- \n ratio is",signif(sum(temp[grp2])/sum(temp[grp1]),3),
#              "(K=",K,")"))
  Ktemp[1,1+jjj] <- round(sum(temp[grp2])/sum(temp[grp1]),3)
}
temp <-  apply(abs(dffitvals),2,median)
oldtemp <- temp
tempunc <- matrix(NA,nrow(uncs),ncol(uncs))
for (obs in 1:nrow(uncs)) tempunc[obs,] <- sqrt( uncs[obs,]^2 / wtadj )
newunc <- tempunc
wtadjpower <- 1

## RUN THE CODE BELOW MULTIPLE TIMES TO FINE TUNE WEIGHT ADJUSTMENTS
for (jj in 1:3)
#for (jj in 4:totreps)
{
  #if (min(newunc) < 0.00000001) newunc <- newunc*0.00000001/min(newunc)
  write(t(newunc),file="newunc.txt",ncol=ncol(newunc))
  write(t(concs),file="concs.txt",ncol=ncol(concs))
  runpmfplain("default_nogkey.ini","concs.txt", "newunc.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0.1,numfacs,nrow(concs),ncol(concs))    
  lam<-read.table(mypaste("testlam",numfacs,".txt"))
  Lambda <- lam
  if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
  if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
  if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
  if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
  if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
  if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
  if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
  if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
  if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

  dffitvals <- matrix(NA,nrow(concs),ncol(concs))
  for (obs in 1:nrow(concs)) dffitvals[obs,] <- 
     dffits(lm(as.double(concs[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(newunc[obs,]^2) ))
  #par(mfrow=c(1,1))
  temp <-  apply(abs(dffitvals),2,median)
  dfplot(dffitvals)
  title(paste("iteration",jj,"-- \n ratio is",signif(sum(temp[grp2])/sum(temp[grp1]),3),
              "(K=",K,")"))
  Ktemp[1+jj,1] <- signif(sum(temp[grp2])/sum(temp[grp1]),3)

  wtadj <- rep(1,ncol(concs))
  for (jjj in 1:finetunereps)
  {
    if (K < Knat) wtadj[grp2] <- wtadj[grp2] * (K*sum(temp[groups$g1])) / sum(temp[groups$g2])
    if (K > Knat) wtadj[grp1] <- wtadj[grp1] * sum(temp[groups$g2]) / (K*sum(temp[groups$g1]))
       #wtadj[grp2] <- wtadj[grp2] * (K*sum(temp[groups$g1])) / 
       #           ( (1-K)*sum(temp[groups$s])+sum(temp[groups$g2]) )
    dffitvals <- matrix(NA,nrow(concs),ncol(concs))
    for (obs in 1:nrow(concs)) dffitvals[obs,] <- 
       dffits(lm(as.double(concs[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
       weights= 1/c(newunc[obs,]^2/wtadj) ))
    #par(mfrow=c(1,1))
    temp <-  apply(abs(dffitvals),2,median)
    #dfplot(dffitvals)
    #title(paste("iteration",jj,": fine-tune \#:",jjj,"-- \n ratio is",signif(sum(temp[grp2])/sum(temp[grp1]),3),
    #          "(K=",K,")"))
    Ktemp[1+jj,1+jjj] <- signif(sum(temp[grp2])/sum(temp[grp1]),3)
  }

  Ktemp
  temp <-  apply(abs(dffitvals),2,median)
  tempunc <- matrix(NA,nrow(uncs),ncol(uncs))
  if (jj > 1) if (sum(Ktemp[jj:(1+jj),1] > K)==1) wtadjpower <- wtadjpower/2
  #if (jj > 1 & abs(Ktemp[1+jj,1]-K)/K > .3) 
    {if (sum(Ktemp[jj:(1+jj),1] > K)==0 | sum(Ktemp[jj:(1+jj),1] > K)==2) wtadjpower <- wtadjpower*1.5}
  #else 
  if (jj > 3 & abs(Ktemp[1+jj,1]-K)/K > .02) 
    {if (sum(Ktemp[(jj-2):(1+jj),1] > K)==0 | sum(Ktemp[(jj-2):(1+jj),1] > K)==4) wtadjpower <- wtadjpower*1.25}
  wtadjpower <- min(wtadjpower,4)
  wtadjpower <- max(wtadjpower,.01)

  for (obs in 1:nrow(newunc)) tempunc[obs,] <- sqrt( newunc[obs,]^2 / wtadj^wtadjpower )
  oldunc <- newunc
  newunc <- tempunc
  lam<-read.table(mypaste("testlam",numfacs,".txt"))
  fac<-read.table(mypaste("testfac",numfacs,".txt"))
  lamjj <- lamarray[,,1+jj] <- as.matrix(lam)
  facjj <- facarray[,,1+jj] <- as.matrix(fac)
  lamplots2(lamjj,"morg",names=mynames)
  #lamplots(lamjj,"gkey") 
  facplots2(facjj,"morgf",newdates)
  weekendplots3(numfacs,newdates,facname=facjj,"morgw")
  system(mypaste("multipanel_morg",numfacs,".bat"),invisible=TRUE) 
       #see License in this directory for the imagemagick license
}
## END OF WEIGHT ADJUSTMENT CODE

lamfinal <- lamjj
facfinal <- facjj
uncfinal <- oldunc
i <- numfacs

  lamplots2(lamfinal,"morg",names=mynames)
  #lamplotslog(lamfinal,"morg",names=mynames)
  facplots2(facfinal,"morgf",newdates)
  weekendplots3(numfacs,newdates,facname=facfinal,"morgw")
  system(mypaste("multipanel_morg",numfacs,".bat"),invisible=TRUE) 



totreps <- 15  # 12 actually used
## Compare facs and lams across time series
faccomps <- matrix(NA,totreps+1,totreps+1)
for (ii in 2:(totreps+1)) for (j in 1:(ii-1)) faccomps[ii,j] <- bestcormat(facarray[,,ii],facarray[,,j])
round(faccomps,3)

lamcomps <- matrix(NA,totreps+1,totreps+1)
for (ii in 2:(totreps+1)) for (j in 1:(ii-1)) lamcomps[ii,j] <- bestrevcormat(lamarray[,,ii],lamarray[,,j])
round(lamcomps,3)

#lamcompsdumb <- matrix(NA,totreps+1,totreps+1)
#for (i in 2:(totreps+1)) for (j in 1:(ii-1)) lamcompsdumb[i,j] <- bestcormat(lamarray[,,i],lamarray[,,j])
#round(lamcompsdumb,3)

lamfinal.1 <- lamfinal
facfinal.1 <- facfinal
uncfinal.1 <- uncfinal
lamarray.1 <- lamarray
facarray.1 <- facarray
Ktemp.1 <- Ktemp
faccomps.1 <- faccomps
lamcomps.1 <- lamcomps

cor(facfinal,facUNWT)
bestcormat(facfinal,facUNWT)
revcormat(lamfinal,lamUNWT)
bestrevcormat(lamfinal,lamUNWT)

cor(facfinal,facorg)
bestcor(facfinal,facorg)
revcormat(lamfinal[-grp1,],lamorg)
bestrevcormat(lamfinal[-grp1,],lamorg)

cor(facfinal,facmet)
bestcor(facfinal,facmet)
revcormat(lamfinal[-grp2,],lammet)
bestrevcormat(lamfinal[-grp2,],lammet)

cor(facfinal.1,facmet)
bestcor(facfinal.1,facmet)
revcormat(lamfinal.1[-grp2,],lammet)
bestrevcormat(lamfinal.1[-grp2,],lammet)

cor(facfinal.01,facmet)
bestcor(facfinal.01,facmet)
revcormat(lamfinal.01[-grp2,],lammet)
bestrevcormat(lamfinal.01[-grp2,],lammet)

cor(facfinal.0000001,facmet)
bestcor(facfinal.0000001,facmet)
revcormat(lamfinal.0000001[-grp2,],lammet)
bestrevcormat(lamfinal.0000001[-grp2,],lammet)

cor(facfinal.01,facfinal.0000001)
bestcor(facfinal.01,facfinal.0000001)
revcormat(lamfinal.01,lamfinal.0000001)
bestrevcormat(lamfinal.01,lamfinal.0000001)
revcormat(lamfinal.01[-grp2,],lamfinal.0000001[-grp2,])
bestrevcormat(lamfinal.01[-grp2,],lamfinal.0000001[-grp2,])

cor(facfinal100000,facorg)
bestcormat(facfinal100000,facorg)
revcormat(lamfinal100000[-grp1,],lamorg)
bestrevcormat(lamfinal100000[-grp1,],lamorg)

cor(facfinal100,facorg)
bestcormat(facfinal100,facorg)
revcormat(lamfinal100[-grp1,],lamorg)
bestrevcormat(lamfinal100[-grp1,],lamorg)

cor(facfinal10,facfinal100000)
revcormat(lamfinal10,lamfinal100000)
revcormat(lamfinal10[-grp1,],lamfinal100000[-grp1,])


junk <- facfinal100000[,2]
names(junk) <- newdates
junk








# For presentation
round(cor(facmet,facorg),3)
round(bestcor(facmet,facorg),3)

round(cor(facmet,facUNWT),3)
round(bestcor(facmet,facUNWT),3)
round(revcormat(lammet,lamUNWT[-grp2,]),3)
round(bestrevcormat(lammet,lamUNWT[-grp2,]),3)

round(cor(facorg,facUNWT),3)
round(bestcor(facorg,facUNWT),3)
round(revcormat(lamorg,lamUNWT[-grp1,]),3)
round(bestrevcormat(lamorg,lamUNWT[-grp1,]),3)


#Corrs with facmet
corrswithfacmet <- c(
1,
bestcor(facfinal.0000001,facmet),
bestcor(facfinal.01,facmet),
bestcor(facUNWT,facmet),  # unwt yields influence of 0.0358
bestcor(facfinal.1,facmet),
bestcor(facfinal1,facmet),
bestcor(facfinal10,facmet),
bestcor(facfinal100,facmet),
bestcor(facfinal100000,facmet),
bestcor(facorg,facmet))

#Corrs with lammet
corrswithlammet <- c(
1,
bestrevcormat(lamfinal.0000001[-grp2,],lammet),
bestrevcormat(lamfinal.01[-grp2,],lammet),
bestrevcormat(lamUNWT[-grp2,],lammet),  # unwt yields influence of 0.0358
bestrevcormat(lamfinal.1[-grp2,],lammet),
bestrevcormat(lamfinal1[-grp2,],lammet),
bestrevcormat(lamfinal10[-grp2,],lammet),
bestrevcormat(lamfinal100[-grp2,],lammet),
bestrevcormat(lamfinal100000[-grp2,],lammet),
NA)


#ExpCorrs with lammet
expcorrswithlammet <- c(
1,
bestexpcormat(lamfinal.0000001[-grp2,],lammet,facfinal.0000001,facmet),
bestexpcormat(lamfinal.01[-grp2,],lammet,facfinal.01,facmet),
bestexpcormat(lamUNWT[-grp2,],lammet,facUNWT,facmet),  # unwt yields influence of 0.0358
bestexpcormat(lamfinal.1[-grp2,],lammet,facfinal.1,facmet),
bestexpcormat(lamfinal1[-grp2,],lammet,facfinal1,facmet),
bestexpcormat(lamfinal10[-grp2,],lammet,facfinal10,facmet),
bestexpcormat(lamfinal100[-grp2,],lammet,facfinal100,facmet),
bestexpcormat(lamfinal100000[-grp2,],lammet,facfinal100000,facmet),
NA)




#Corrs with facUNWT
corrswithfacUNWT <- c(
bestcor(facmet,facUNWT),
bestcor(facfinal.0000001,facUNWT),
bestcor(facfinal.01,facUNWT),
bestcor(facUNWT,facUNWT),  # unwt yields influence of 0.0358
bestcor(facfinal.1,facUNWT),
bestcor(facfinal1,facUNWT),
bestcor(facfinal10,facUNWT),
bestcor(facfinal100,facUNWT),
bestcor(facfinal100000,facUNWT),
bestcor(facorg,facUNWT)
)

#Corrs with lamUNWT
corrswithlamUNWT <- c(
bestrevcormat(lammet,lamUNWT[-grp2,]),
bestrevcormat(lamfinal.0000001,lamUNWT),
bestrevcormat(lamfinal.01,lamUNWT),
bestrevcormat(lamUNWT,lamUNWT),  # unwt yields influence of 0.0358
bestrevcormat(lamfinal.1,lamUNWT),
bestrevcormat(lamfinal1,lamUNWT),
bestrevcormat(lamfinal10,lamUNWT),
bestrevcormat(lamfinal100,lamUNWT),
bestrevcormat(lamfinal100000,lamUNWT),
bestrevcormat(lamorg,lamUNWT[-grp1,]),
)


#ExpCorrs with lamUNWT
expcorrswithlamUNWT <- c(
bestexpcormat(lammet,lamUNWT[-grp2,],facmet,facUNWT),
bestexpcormat(lamfinal.0000001,lamUNWT,facfinal.0000001,facUNWT),
bestexpcormat(lamfinal.01,lamUNWT,facfinal.01,facUNWT),
bestexpcormat(lamUNWT,lamUNWT,facUNWT,facUNWT),  # unwt yields influence of 0.0358
bestexpcormat(lamfinal.1,lamUNWT,lamfinal.1,lamUNWT),
bestexpcormat(lamfinal1,lamUNWT,facfinal1,facUNWT),
bestexpcormat(lamfinal10,lamUNWT,facfinal10,facUNWT),
bestexpcormat(lamfinal100,lamUNWT,facfinal100,facUNWT),
bestexpcormat(lamfinal100000,lamUNWT,facfinal100000,facUNWT),
bestexpcormat(lamorg,lamUNWT[-grp1,],facorg,facUNWT),
)




#Corrs with facorg
corrswithfacorg <- c(
bestcor(facmet,facorg),
bestcor(facfinal.0000001,facorg),
bestcor(facfinal.01,facorg),
bestcor(facUNWT,facorg),  # unwt yields influence of 0.0358
bestcor(facfinal.1,facorg),
bestcor(facfinal1,facorg),
bestcor(facfinal10,facorg),
bestcor(facfinal100,facorg),
bestcor(facfinal100000,facorg),
1)

#Corrs with lamorg
corrswithlamorg <- c(
NA,
bestrevcormat(lamfinal.0000001[-grp1,],lamorg),
bestrevcormat(lamfinal.01[-grp1,],lamorg),
bestrevcormat(lamUNWT[-grp1,],lamorg),  # unwt yields influence of 0.0358
bestrevcormat(lamfinal.1[-grp1,],lamorg),
bestrevcormat(lamfinal1[-grp1,],lamorg),
bestrevcormat(lamfinal10[-grp1,],lamorg),
bestrevcormat(lamfinal100[-grp1,],lamorg),
bestrevcormat(lamfinal100000[-grp1,],lamorg),
1)

#ExpCorrs with lamorg
expcorrswithlamorg <- c(
NA,
bestexpcormat(lamfinal.0000001[-grp1,],lamorg,facmet,facorg),
bestexpcormat(lamfinal.01[-grp1,],lamorg,facfinal.01,facorg),
bestexpcormat(lamUNWT[-grp1,],lamorg,facUNWT,facorg),  # unwt yields influence of 0.0358
bestexpcormat(lamfinal.1[-grp1,],lamorg,facfinal.1,facorg),
bestexpcormat(lamfinal1[-grp1,],lamorg,facfinal1,facorg),
bestexpcormat(lamfinal10[-grp1,],lamorg,facfinal10,facorg),
bestexpcormat(lamfinal100[-grp1,],lamorg,facfinal100,facorg),
bestexpcormat(lamfinal100000[-grp1,],lamorg,facfinal100000,facorg),
1)


# K's
Ks <- c(.00000001,.0000001, .01, .0358, .1, 1, 10, 100, 100000, 1000000)
plot(Ks,corrswithfacUNWT,ylim=c(0,1),ylab="Similarity",
   xlab="K" ,xaxt="n",log="x",type="l")
lines(Ks,corrswithfacmet,col="blue")
lines(Ks,corrswithfacorg,col="red")
axis(1,at=Ks,labels=c("MET",Ks[2:(length(Ks)-1)],"ORG"),las=2,cex=.3)
title("Contribution similarities with UNWT (black), MET (blue), and ORG (red)") 
abline(v=0.0358,lty=2)
 
plot(Ks,corrswithlamUNWT,ylim=c(-.2,1),ylab="Similarity",
   xlab="K" ,xaxt="n",log="x",type="l")
lines(Ks,corrswithlammet,col="blue")
lines(Ks,corrswithlamorg,col="red")
axis(1,at=Ks,labels=c("MET",Ks[2:(length(Ks)-1)],"ORG"),las=2,cex=.3)
title("Profile similarities with UNWT (black), MET (blue), and ORG (red)") 
abline(v=0.0358,lty=2)
abline(h=0,lty=2)

plot(Ks,expcorrswithlamUNWT,ylim=c(0,1),ylab="Similarity",
   xlab="K" ,xaxt="n",log="x",type="l")
lines(Ks,expcorrswithlammet,col="blue")
lines(Ks,expcorrswithlamorg,col="red")
axis(1,at=Ks,labels=c("MET",Ks[2:(length(Ks)-1)],"ORG"),las=2,cex=.3)
title("Profile similarities with UNWT (black), MET (blue), and ORG (red)") 
abline(v=0.0358,lty=2)
#abline(h=0,lty=2)

cbind(Ks,expcorrswithlamUNWT,expcorrswithlammet,expcorrswithlamorg)











###
### Consider influence of a species after removing species in similar group
###
Lambda <- lamUNWT
if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

dffitvals <- matrix(NA,nrow(concs),ncol(concs))
for (obs in 1:nrow(concs)) dffitvals[obs,] <- 
   dffits(lm(as.double(concs[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(uncs[obs,]^2) ))

tempold <- apply(abs(dffitvals),2,median)
temp <-  apply(abs(dffitvals),2,median)
dfplot(dffitvals)

avecorr(cor(morgc[,45:82]),c(rep("1.diacids",12),rep("2.hop-stear",16),rep("3.pahs",9),"4.levogluc"))
 
elimgroups <- c(rep("01.metals",40),rep("02.carbon",2),rep("03.ions",2),rep("1.diacids",12),
    rep("2.hop-stear",16),rep("3.pahs",9),"4.levogluc")
avecorr(cor(morgc),elimgroups)
 

dffitvalsNOGROUP <- matrix(NA,nrow(concs),ncol(concs))

#spec <- 50

for (spec in 1:ncol(concs))
for (obs in 1:nrow(concs)) 
{   
  nogos <- uniquify(1:ncol(concs),spec,
                  (1:length(elimgroups))[elimgroups==elimgroups[spec]])$g2 
  dffitvalsNOGROUP[obs,spec] <- 
     dffits(lm(as.double(concs[obs,-nogos]) ~ 
          l1[-nogos] + l2[-nogos] + l3[-nogos] + l4[-nogos] + l5[-nogos] + l6[-nogos] + l7[-nogos] - 1, 
          weights= 1/ c(as.double(uncs[obs,-nogos])^2) ))[ (1:82)[-nogos] == spec ]
}
}

dffitvalsNOGROUP[,82] <- 0
tempold <- apply(abs(dffitvalsNOGROUP),2,median)
temp <-  apply(abs(dffitvals),2,median)
dfplot(dffitvalsNOGROUP)




 
elimgroups2 <- c(rep("01.metals",40),rep("02.carbon",2),rep("03.ions",2),rep("04.organics",38))
avecorr(cor(morgc),elimgroups2)
 

dffitvalsNOGROUP2 <- matrix(NA,nrow(concs),ncol(concs))

#spec <- 50

for (spec in 1:ncol(concs))
for (obs in 1:nrow(concs)) 
{   
  nogos <- uniquify(1:ncol(concs),spec,
                  (1:length(elimgroups2))[elimgroups2==elimgroups2[spec]])$g2 
  dffitvalsNOGROUP2[obs,spec] <- 
     dffits(lm(as.double(concs[obs,-nogos]) ~ 
          l1[-nogos] + l2[-nogos] + l3[-nogos] + l4[-nogos] + l5[-nogos] + l6[-nogos] + l7[-nogos] - 1, 
          weights= 1/ c(as.double(uncs[obs,-nogos])^2) ))[ (1:82)[-nogos] == spec ]
}
}

#dffitvalsNOGROUP[,82] <- 0
tempold2 <- apply(abs(dffitvalsNOGROUP2),2,median)
tempold <- apply(abs(dffitvalsNOGROUP),2,median)
temp <-  apply(abs(dffitvals),2,median)

par(mfrow=c(3,1))
dfplot(dffitvals)
dfplot(dffitvalsNOGROUP)
dfplot(dffitvalsNOGROUP2)

## NOT MAKING A BIG DIFFERENCE








###
### Consider influence after actually dropping specie from list
###

####### USING BIGMET ##########
i <- 7
  write(t(slunc2),file="slunc2.txt",ncol=ncol(slunc2))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
metnames <- dimnames(morgc)[[2]][1:44]

runpmfplain("default_nogkey.ini","sldat2.txt","slunc2.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))
lammetall <- lam
facmetall <- fac

Lambda <- lam
if (ncol(Lambda) >= 1) l1 <- as.double(Lambda[,1])
if (ncol(Lambda) >= 2) l2 <- as.double(Lambda[,2])
if (ncol(Lambda) >= 3) l3 <- as.double(Lambda[,3])
if (ncol(Lambda) >= 4) l4 <- as.double(Lambda[,4])
if (ncol(Lambda) >= 5) l5 <- as.double(Lambda[,5])
if (ncol(Lambda) >= 6) l6 <- as.double(Lambda[,6])
if (ncol(Lambda) >= 7) l7 <- as.double(Lambda[,7])
if (ncol(Lambda) >= 8) l8 <- as.double(Lambda[,8])
if (ncol(Lambda) >= 9) l9 <- as.double(Lambda[,9])

dffitvals <- matrix(NA,nrow(sldat2),ncol(sldat2))
for (obs in 1:nrow(sldat2)) dffitvals[obs,] <- 
   dffits(lm(as.double(sldat2[obs,]) ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 - 1, 
          weights= 1/ c(slunc2[obs,]^2) ))

tempold <- apply(abs(dffitvals),2,median)

#par(mfrow=c(3,4))
dfplot(dffitvals,namespec=metnames)
#  title(paste("iteration 0-- ratio is",round(sum(tempold[41:44])/sum(tempold[1:40]),3),
#              "(K=",K,")"))



#  write(t(morgc),file="morgcjj.txt",ncol=ncol(morgc))
#  write(t(morgu),file="morgujj.txt",ncol=ncol(morgu))
#  concs <- morgc
#  uncs <- morgu
#mynamesjj <- dimnames(morgc)[[2]]
#mynamesjj
#names <- mynamesjj
#jj <- 17

lamjjarray <- array(NA,dim=c(43,7,44))
facjjarray <- array(NA,dim=c(661,7,44))
#lamjjarray2[,,1:40] <- lamjjarray[,,2:41]
#facjjarray2[,,1:40] <- facjjarray[,,2:41]
lamjjarray <- lamjjarray2
facjjarray <- facjjarray2

for (jj in 43:44)   ### DECIDE WHAT TO DO WITH ESTIMATES...SAVE lamjj[,,] and facjj[,,]...
{

#jj <- jj + 1
numfacs <- i <- 7

  write(t(sldat2[,-jj]),file="metalscjj.txt",ncol=ncol(sldat2[,-jj]))
  write(t(slunc2[,-jj]),file="metalsujj.txt",ncol=ncol(slunc2[,-jj]))
  concs <- sldat2[,-jj]
  uncs <- slunc2[,-jj]
mynamesjj <- metnames[-jj]
mynamesjj
names <- mynamesjj
#newdates <- as.Date(dimnames(sldat2)[[1]], "%d%b%Y")
newdatesall

# datevec <- newdates
runpmfplain("default_nogkey.ini","metalscjj.txt","metalsujj.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0.1,i,nrow(concs),ncol(concs))    

lam <-as.matrix(read.table(mypaste("testlam",i,".txt")))
lamjjarray[,,jj] <- lam
fac <- as.matrix(read.table(mypaste("testfac",i,".txt")))
facjjarray[,,jj] <- fac
print(paste("@@@@@@@@@@@@@@@@ SPECIES #",jj," out of 44 @@@@@@@@@@@@@@@@@@@@@@@"))
write(jj,"1jj.txt")
}












## JUNK
org <- as.matrix(morgc[,c(3,4,41,42,45:82)])
dim(org)
dim(lamorg)
dim(facorg)

plot(c(as.matrix(org)),c(as.matrix(facorg)%*%t(as.matrix(lamorg))))
plot(c( (as.matrix(facorg)%*%t(as.matrix(lamorg)))[,-3] ),c(as.matrix(facorg)%*%t(as.matrix(lamorg[-3,]))))
abline(a=0,b=1)
lamorgi <- as.matrix(lamorg[-3,])
apply(lamorgi,2,sum)
facorginew <- org[,-3] %*% solve(lamorgi %*% t(lamorgi)) %*% lamorgi
names

plot(c( (as.matrix(facUNWT)%*%t(as.matrix(lamUNWT)))[,-jj] ),c(as.matrix(fac)%*%t(as.matrix(lam))))
plot(c( as.matrix(facUNWT) ),c(as.matrix(fac)))
plot(c( as.matrix(lamUNWT[-jj,]) ),c(as.matrix(lam)))

par(mfrow=c(2,4))
facnew <- fac
lamnew <- lam
lamUNWTnew <- lamUNWT
facnew <- fac[,c(1,2,4,3,5,6,7)]
lamnew <- lam[,c(1,2,4,3,5,6,7)]
lamUNWTnew <- lamUNWT[,c(1,2,4,3,5,6,7)]
facnewadj <- facnew
#facUNWTjj[,ii] <- lamnew[,ii]/sum(lamnew[,ii])
junk <- rep(NA,7)
for (ii in 1:7) junk[ii] <- (1 - lamUNWT[jj,ii])
for (ii in 1:7) facnewadj[,ii] <- facnew[,ii] / (1 - lamUNWT[jj,ii])

par(mfrow=c(2,4))
for (ii in 1:7)
#ii <- 0
#ii <- ii+1
{
  plot(facUNWT[,ii],type="l",col="blue",lwd=3)
  lines(facnew[,ii],col="red")
  #lines(facnewadj[,ii],col="green")
  title(ii)
}

par(mfrow=c(3,1))
  plot( apply(morgc,1,sum) , ylim=c(0,32)) 
  lines( apply(as.matrix(facUNWT)%*%t(as.matrix(lamUNWT)) ,1,sum) , col="red")
#  lines(facnew[,ii],col="red")
  plot( apply(morgc[,-jj],1,sum) , col="blue", ylim=c(0,32))
  lines( apply(as.matrix(fac)%*%t(as.matrix(lam)) ,1,sum) , col="green")
  plot( apply(facUNWT,1,sum) , type="l", col="black")
  points( apply(fac,1,sum) , type="l", col="red")

  plot( apply(morgc,1,sum) , ylim=c(0,32), type="l", col="grey", lwd=3) 
  lines( apply(morgc[,-jj],1,sum) , col="lightblue", ylim=c(0,32), type="l", lwd=3)
  plot( apply(facUNWT,1,sum) , type="l", col="green")
  points( apply(fac,1,sum) , type="l", col="red")




facorgi <- 


######### USING METORG ###############
  write(t(morgc),file="morgcjj.txt",ncol=ncol(morgc))
  write(t(morgu),file="morgujj.txt",ncol=ncol(morgu))
  concs <- morgc
  uncs <- morgu
mynamesjj <- dimnames(morgc)[[2]]
mynamesjj
names <- mynamesjj
jj <- 17

lamjjarray <- array(NA,dim=c(81,7,82))
facjjarray <- array(NA,dim=c(102,7,82))

#for (jj in 1:82)   ### DECIDE WHAT TO DO WITH ESTIMATES...SAVE lamjj[,,] and facjj[,,]...
#{

jj <- jj + 1
numfacs <- i <- 7

  write(t(morgc[,-jj]),file="morgcjj.txt",ncol=ncol(morgc[,-jj]))
  write(t(morgu[,-jj]),file="morgujj.txt",ncol=ncol(morgu[,-jj]))
  concs <- morgc[,-jj]
  uncs <- morgu[,-jj]
mynamesjj <- dimnames(morgc)[[2]][-jj]
mynamesjj
names <- mynamesjj
newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
newdates

# datevec <- newdates
runpmfplain("default_nogkey.ini","morgcjj.txt","morgujj.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0.5,i,nrow(concs),ncol(concs))    

lam <-read.table(mypaste("testlam",i,".txt"))
lamjjarray[,,jj] <-read.table(mypaste("testlam",i,".txt"))
fac <- read.table(mypaste("testfac",i,".txt"))
facjjarray[,,jj] <- read.table(mypaste("testfac",i,".txt"))
print(paste("@@@@@@@@@@@@@@@@ SPECIES #",jj," out of 82 @@@@@@@@@@@@@@@@@@@@@@@"))

#}












## JUNK
org <- as.matrix(morgc[,c(3,4,41,42,45:82)])
dim(org)
dim(lamorg)
dim(facorg)

plot(c(as.matrix(org)),c(as.matrix(facorg)%*%t(as.matrix(lamorg))))
plot(c( (as.matrix(facorg)%*%t(as.matrix(lamorg)))[,-3] ),c(as.matrix(facorg)%*%t(as.matrix(lamorg[-3,]))))
abline(a=0,b=1)
lamorgi <- as.matrix(lamorg[-3,])
apply(lamorgi,2,sum)
facorginew <- org[,-3] %*% solve(lamorgi %*% t(lamorgi)) %*% lamorgi
names

plot(c( (as.matrix(facUNWT)%*%t(as.matrix(lamUNWT)))[,-jj] ),c(as.matrix(fac)%*%t(as.matrix(lam))))
plot(c( as.matrix(facUNWT) ),c(as.matrix(fac)))
plot(c( as.matrix(lamUNWT[-jj,]) ),c(as.matrix(lam)))

par(mfrow=c(2,4))
facnew <- fac
lamnew <- lam
lamUNWTnew <- lamUNWT
facnew <- fac[,c(1,2,4,3,5,6,7)]
lamnew <- lam[,c(1,2,4,3,5,6,7)]
lamUNWTnew <- lamUNWT[,c(1,2,4,3,5,6,7)]
facnewadj <- facnew
#facUNWTjj[,ii] <- lamnew[,ii]/sum(lamnew[,ii])
junk <- rep(NA,7)
for (ii in 1:7) junk[ii] <- (1 - lamUNWT[jj,ii])
for (ii in 1:7) facnewadj[,ii] <- facnew[,ii] / (1 - lamUNWT[jj,ii])

par(mfrow=c(2,4))
for (ii in 1:7)
#ii <- 0
#ii <- ii+1
{
  plot(facUNWT[,ii],type="l",col="blue",lwd=3)
  lines(facnew[,ii],col="red")
  #lines(facnewadj[,ii],col="green")
  title(ii)
}

par(mfrow=c(3,1))
  plot( apply(morgc,1,sum) , ylim=c(0,32)) 
  lines( apply(as.matrix(facUNWT)%*%t(as.matrix(lamUNWT)) ,1,sum) , col="red")
#  lines(facnew[,ii],col="red")
  plot( apply(morgc[,-jj],1,sum) , col="blue", ylim=c(0,32))
  lines( apply(as.matrix(fac)%*%t(as.matrix(lam)) ,1,sum) , col="green")
  plot( apply(facUNWT,1,sum) , type="l", col="black")
  points( apply(fac,1,sum) , type="l", col="red")

  plot( apply(morgc,1,sum) , ylim=c(0,32), type="l", col="grey", lwd=3) 
  lines( apply(morgc[,-jj],1,sum) , col="lightblue", ylim=c(0,32), type="l", lwd=3)
  plot( apply(facUNWT,1,sum) , type="l", col="green")
  points( apply(fac,1,sum) , type="l", col="red")




facorgi <- 




## Look at performance of PMF in neighborhood of final uncertainties

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

# Check: (1) What is variability of outputs as uncertainties changed slightly?  If small within a neighborhood
#  of uncertainties, we can ignore slight changes in weights as long as variability of outputs in the time 
#  series of attempts is similar.  Check convergence by looking at variability of outputs in a range similar
#  to variability for uncertainty neighborhood of any particular solution (try uncertainties +/- sd equal to
#  1%or5% of uncertainty) 



##
##  LOOK AT VARYING AMOUNTS OF ERRORS ON UNCERTAINTIES
##

## OLD MORG ANALYSIS ##
numfacs <- i <- 7

  write(t(morgc),file="morgc.txt",ncol=ncol(morgc))
  write(t(morgu),file="morgu.txt",ncol=ncol(morgu))
  concs <- morgc
  uncs <- morgu
mynames <- dimnames(morgc)[[2]]
mynames
names <- mynames
newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
newdates

cvvals <- c(.01,.05,.15)
reperrorvals <- 20

lamarrayerrs <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrs <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))
#lamarrayerrsold <- lamarrayerrs
#facarrayerrsold <- facarrayerrs
for (cvi in 1:length(cvvals))
{
  for (reperri in 16:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    #uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","morgc.txt", "uncerrs.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0.1,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrs[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrs[,,cvi,reperri] <- as.matrix(fac)
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


cvi <- 1
faccompsall1 <- matrix(NA,reperri,numfacs)
lamcompsall1 <- matrix(NA,reperri,numfacs)
faccompsallrmse1 <- matrix(NA,reperri,numfacs)
lamcompsallexp1 <- matrix(NA,reperri,numfacs)

for (i in 1:(reperri)) 
{
  #faccompsall1[i,] <- bestcormatcors(facarrayerrs[,,cvi,i],facUNWT)
  #lamcompsall1[i,] <- bestrevcormatcors(lamarrayerrs[,,cvi,i],lamUNWT)
  faccompsallrmse1[i,] <- bestrmsematcors(facarrayerrs[,,cvi,i],facUNWT) ## now does aae
  #lamcompsallexp1[i,] <- bestexpcormatcors(lamarrayerrs[,,cvi,i],lamUNWT,facarrayerrs[,,cvi,i],facUNWT)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall1))
boxplot(data.frame(lamcompsall1))





cvi <- 2
#faccompsall2 <- matrix(NA,reperri,numfacs)
#lamcompsall2 <- matrix(NA,reperri,numfacs)
faccompsallrmse2 <- matrix(NA,reperri,numfacs)
lamcompsallexp2 <- matrix(NA,reperri,numfacs)

for (i in 1:(reperri)) 
{
  #faccompsall2[i,] <- bestcormatcors(facarrayerrs[,,cvi,i],facUNWT)
  #lamcompsall2[i,] <- bestrevcormatcors(lamarrayerrs[,,cvi,i],lamUNWT)
  faccompsallrmse2[i,] <- bestrmsematcors(facarrayerrs[,,cvi,i],facUNWT)## now does aae
  #lamcompsallexp2[i,] <- bestexpcormatcors(lamarrayerrs[,,cvi,i],lamUNWT,facarrayerrs[,,cvi,i],facUNWT)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall2))
boxplot(data.frame(lamcompsall2))




cvi <- 3
#faccompsall3 <- matrix(NA,reperri,numfacs)
#lamcompsall3 <- matrix(NA,reperri,numfacs)
faccompsallrmse3 <- matrix(NA,reperri,numfacs)
lamcompsallexp3 <- matrix(NA,reperri,numfacs)

for (i in 1:(reperri)) 
{
  #faccompsall3[i,] <- bestcormatcors(facarrayerrs[,,cvi,i],facUNWT)
  #lamcompsall3[i,] <- bestrevcormatcors(lamarrayerrs[,,cvi,i],lamUNWT)
  faccompsallrmse3[i,] <- bestrmsematcors(facarrayerrs[,,cvi,i],facUNWT)## now does aae
  #lamcompsallexp3[i,] <- bestexpcormatcors(lamarrayerrs[,,cvi,i],lamUNWT,facarrayerrs[,,cvi,i],facUNWT)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall3))
boxplot(data.frame(lamcompsall3))



par(mfcol=c(2,3))
boxplot(data.frame(faccompsall2),ylim=c(.5,1),xlab="Contributions")
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsall2),ylim=c(.3,1),xlab="Profiles")

boxplot(data.frame(faccompsall3),ylim=c(.5,1),xlab="Contributions")
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsall3),ylim=c(.3,1),xlab="Profiles")

boxplot(data.frame(faccompsall1),ylim=c(.5,1),xlab="Contributions")
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsall1),ylim=c(.3,1),xlab="Profiles")



par(mfcol=c(2,3))
#boxplot(data.frame(faccompsallrmse2),ylim=c(0,2.5),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse2)/apply(facUNWT,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%","1000%"),at=c(.005,.01,.05,.1,.5,1,5,10),las=2)
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)
title("Uncertainties from N(true, .05*true):\n IQR is (.97,1.03)")
boxplot(data.frame(lamcompsallexp2),ylim=c(.3,1),xlab="Profiles",xaxt="n")
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)

#boxplot(data.frame(faccompsallrmse3),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse3)/apply(facUNWT,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse3)/apply(facUNWT,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%","1000%"),at=c(.005,.01,.05,.1,.5,1,5,10),las=2)
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)
title("Uncertainties from N(true, .15*true):\n IQR is (.90,1.10)")
boxplot(data.frame(lamcompsallexp3),ylim=c(.3,1),xlab="Profiles",xaxt="n")
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)

#boxplot(data.frame(faccompsallrmse1),ylim=c(0,2.5),xlab="Contributions")
#boxplot(data.frame(t(t(faccompsallrmse1)/apply(facUNWT,2,mean)) ),ylim=c(0,2),xlab="Contributions")
boxplot(data.frame(t(t(faccompsallrmse1)/apply(facUNWT,2,mean)) ),
      log="y",ylim=c(0.0025,8.1),xlab="Contributions",yaxt="n",xaxt="n")
axis(2,labels=c("0.5%","1%","5%","10%","50%","100%","500%","1000%"),at=c(.005,.01,.05,.1,.5,1,5,10),las=2)
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)
title("Uncertainties from Unif(.5*true, 1.5*true):\n IQR is (.75,1.25)")
boxplot(data.frame(lamcompsallexp1),ylim=c(.3,1),xlab="Profiles",xaxt="n")
axis(1,labels=c("Zn","Cu","Wood","Summer","Mobile","EC","Winter"),at=1:7)





















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


##
##
##
## More levels of cvi
##
##
##
cvvals <- c(.001,.01)
reperrorvals <- 50

lamarrayerrsMORGb <- array(NA,dim=c(ncol(concs),numfacs,length(cvvals),reperrorvals))
facarrayerrsMORGb <- array(NA,dim=c(nrow(concs),numfacs,length(cvvals),reperrorvals))

cvi <- 2
for (cvi in 1:length(cvvals))
{
  for (reperri in 1:reperrorvals)
  #for (reperri in 43:reperrorvals)
  {
    #reperri <- 0
    #reperri <- reperri + 1
    uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,cvvals[cvi]),nrow(uncs),ncol(uncs))
    #if (cvi==1) uncerrs <- uncs * matrix( runif(nrow(uncs)*ncol(uncs),.5,1.5) ,
    #                                      nrow(uncs),ncol(uncs)) #cvi=1
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","morgc.txt", "uncerrs.txt" ,mypaste("testfac",numfacs,".txt"),
     mypaste("testlam",numfacs,".txt"),
     mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrsMORGb[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrsMORGb[,,cvi,reperri] <- as.matrix(fac)
     reperri
  }
}

date()

cvi <- 1
faccompsall4MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsall4MORG <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse4MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp4MORG <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall4MORG[i,] <- bestcormatcors(facarrayerrsMORGb[,,cvi,i],facORIGMORG)
  lamcompsall4MORG[i,] <- bestrevcormatcors(lamarrayerrsMORGb[,,cvi,i],lamORIGMORG)
  faccompsallrmse4MORG[i,] <- bestrmsematcors(facarrayerrsMORGb[,,cvi,i],facORIGMORG) ## now does aae
  lamcompsallexp4MORG[i,] <- bestexpcormatcors(lamarrayerrsMORGb[,,cvi,i],lamORIGMORG,
                                               facarrayerrsMORGb[,,cvi,i],facORIGMORG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall4MORG))
boxplot(data.frame(lamcompsall4MORG))

date()


cvi <- 2
faccompsall5MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsall5MORG <- matrix(NA,reperrorvals,numfacs)
faccompsallrmse5MORG <- matrix(NA,reperrorvals,numfacs)
lamcompsallexp5MORG <- matrix(NA,reperrorvals,numfacs)

for (i in 1:(reperrorvals)) 
{
  faccompsall5MORG[i,] <- bestcormatcors(facarrayerrsMORGb[,,cvi,i],facORIGMORG)
  lamcompsall5MORG[i,] <- bestrevcormatcors(lamarrayerrsMORGb[,,cvi,i],lamORIGMORG)
  faccompsallrmse5MORG[i,] <- bestrmsematcors(facarrayerrsMORGb[,,cvi,i],facORIGMORG) ## now does aae
  lamcompsallexp5MORG[i,] <- bestexpcormatcors(lamarrayerrsMORGb[,,cvi,i],lamORIGMORG,
                                               facarrayerrsMORGb[,,cvi,i],facORIGMORG)
}
par(mfrow=c(1,2))
boxplot(data.frame(faccompsall5MORG))
boxplot(data.frame(lamcompsall5MORG))

date()
##
##
##
## END OF: More levels of cvi
##
##
##





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










pdf("MetalsContributions.pdf",width=12,height=5)
par(mfrow=c(1,3))
junk <- apply(facarrayerrsM[,,,],c(2,3,4),mean)
ymaxcont <- max(junk)
boxplot(data.frame(t( apply(facarrayerrsM[,,2,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .05*true)"))
points(1:9,apply(facORIG,2,mean),pch=10,cex=3)

boxplot(data.frame(t( apply(facarrayerrsM[,,3,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from N(true, .15*true)"))
points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)

boxplot(data.frame(t( apply(facarrayerrsM[,,1,],c(2,3),mean) )), xaxt="n", ylim=c(0,ymaxcont))  # N , k, cvi, rep 
axis(1,labels=facnames9,at=1:9,cex.axis=set.cex)
title(paste("Contributions\nUncertainties from Unif(.5*true, 1.5*true)"))
points(1:9,apply(facORIG,2,mean),col="red",pch="-",cex=3)
dev.off()










### OLD
#cvi <- 3
#faccompsall3 <- matrix(NA,reperri*(reperri-1)/2,numfacs)
#lamcompsall3 <- matrix(NA,reperri*(reperri-1)/2,numfacs)
#
#rownum <- 0
#for (i in 2:(reperri)) for (j in 1:(i-1)) 
#{
#  rownum <- rownum+1
#  faccompsall3[rownum,] <- bestcormatcors(facarrayerrs[,,cvi,i],facarrayerrs[,,cvi,j])
#  lamcompsall3[rownum,] <- bestrevcormatcors(lamarrayerrs[,,cvi,i],lamarrayerrs[,,cvi,j])
#}
#par(mfrow=c(1,2))
#boxplot(data.frame(faccompsall3))
#boxplot(data.frame(lamcompsall3))
#






##
##  LOOK AT VARYING LEVELS OF FPEAK
##


numfacs <- i <- 7

  write(t(morgc),file="morgc.txt",ncol=ncol(morgc))
  write(t(morgu),file="morgu.txt",ncol=ncol(morgu))
  concs <- morgc
  uncs <- morgu
mynames <- dimnames(morgc)[[2]]
mynames
names <- mynames
newdates <- as.Date(dimnames(morgc)[[1]], "%d%b%Y")
newdates

fpeakvals <- seq(-.5,.5,.1)
numfs <- length(fpeakvals)

lamarrayfs <- array(NA,dim=c(ncol(concs),numfacs,numfs))
facarrayfs <- array(NA,dim=c(nrow(concs),numfacs,numfs))

for (cvi in 1:length(fpeakvals))
{
     runpmfplain("default_nogkey.ini","morgc.txt", "morgu.txt" ,mypaste("testfac",numfacs,".txt"),
       mypaste("testlam",numfacs,".txt"),
       mypaste("testq",numfacs,".txt"),fpeak= fpeakvals[cvi],numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrs[,,cvi] <- as.matrix(lam)
     facarrayerrs[,,cvi] <- as.matrix(fac)
  }
  cvi
}



## Randomly generated data
bestcor(matrix(exp(rnorm(100*7)),100,7), matrix(exp(rnorm(100*7)),100,7))
bestcor(matrix((rnorm(100*7)),100,7), matrix((rnorm(100*7)),100,7))


