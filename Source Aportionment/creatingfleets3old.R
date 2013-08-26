setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("SI")
load("SIunc")
load("CI")
load("CIunc")
library(gtools)  ###This library contains functions for Dirichlet Distribution


profilegen<-function()
{
  ### Creating gas profiles ###
  gas<-SI[,4:7]
  gas
  gasunc<-SIunc[,4:7]
  is.matrix(gas)
  gas<-gas[-(121:124),]
  gasdir<-rdirichlet(1,c(1,1,1,1))
  for (i in 1:4)         ###converting profiles from factors to numeric data
  {
    gas[,i]<-as.numeric(as.character(gas[,i]))
    gasunc[,i]<-as.numeric(as.character(gasunc[,i]))
    gas[,i]<-gas[,i]*gasdir[1,i]  
  }
  
  
  gasp<-matrix(NA,120,1)
  for(i in 1:120)
  {
    gasp[i]<-sum(gas[i,])
  }
  
  ### Creating smoker profiles###
  smokers<-SI[,38:50]     
  smokers
  smokerunc<-SIunc[,38:50]
  is.matrix(smokers)
  smokdir<-rdirichlet(1,c(1,1,1,1,1,1,1,1,1,1,1,1,1))
  for (i in 1:13)               ###converting profiles from factors to numeric data
  {
    smokers[,i]<-as.numeric(as.character(smokers[,i]))
    smokerunc[,i]<-as.numeric(as.character(smokerunc[,i]))
    smokers[,i]<-smokers[,i]*smokdir[1,i]
  }
  smokerp<-matrix(,120,1)
  for(i in 1:120)
  {
    smokerp[i]<-sum(smokers[i,])
  }
  
  
  ### Creating diesel profiles###
  diesel<-CI[,8:22]     ###creating profile for diesel vehicles
  diesel
  dieselunc<-CIunc[,8:22]
  is.matrix(diesel)
  diesdir<-rdirichlet(1,c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
  for (i in 1:15)               ###converting profiles from factors to numeric data
  {
    diesel[,i]<-as.numeric(as.character(diesel[,i]))
    dieselunc[,i]<-as.numeric(as.character(dieselunc[,i]))
    diesel[,i]<-diesel[,i]*diesdir[1,i]
  }
  dieselp<-matrix(,120,1)
  for(i in 1:120)
  {
    dieselp[i]<-sum(diesel[i,])
  }
  
  
  yfleet1profiles<-matrix(NA,120,3)
  yfleet1profiles[,1]<-gasp
  yfleet1profiles[,2]<-smokerp
  yfleet1profiles[,3]<-dieselp
  #species<-as.matrix(SI$X)
  #species<-species[-(121:124),]
  
  q=1
  while(q <= dim(yfleet1profiles)[1])
  {                   
    if (sum(yfleet1profiles[q,])==0)
    {
      yfleet1profiles<-yfleet1profiles[-q,];
      #species<-species[-q,];
      q=q-1;
    }
    q=q+1
  }
  #yfleet1profiles
  out<-NULL
  out[[1]]<-yfleet1profiles
  out[[2]]<-gas
  out[[3]]<-smokers
  out[[4]]<-diesel
  return(out)
}

#### Producing Data###


datagen<- function(cv,cvlog,gas,smokers,diesel)
{
  fgas1<-vector(mode="numeric",4)
  fgas<-matrix(NA,120,4)
  for(g in 1:4)
  {
    fgas1[g]<-rnorm(1,15,sd=15*cv)
  }
  for(k in 1:120)
  {
    for(g in 1:4)
    {                     ###All fgas(ji) are drawn from a normal(15,var=9) distribution
      fgas[k,g]<-fgas1[g]
    }
  }
  
  gasmatrix<-matrix(NA,120,4)
  for(k in 1:120)
  {
    for(g in 1:4)
    {
      gasmatrix[k,g]<-gas[k,g]*fgas[k,g]
    }
  }
  
  gasprofile<-vector(mode='numeric',120)
  gasVA<-vector(mode='numeric',120)
  for(k in 1:120)
  {
    gasprofile[k]<-sum(gasmatrix[k,])
    gasVA[k]<-(1/16)*sum((gasunc[k,])^2)
  }
  
  
  fsmokers1<-vector(mode="numeric",13)
  fsmokers<-matrix(NA,120,13)
  for(s in 1:13)
  {
    fsmokers1[s]<-rnorm(1,.937,sd=.937*cv)
  }
  for(k in 1:120)
  {
    for(s in 1:13)
    {                     ###All fsmokers(ji) are drawn from a normal(.937,var=.03515625) distribution
      fsmokers[k,s]<-fsmokers1[s]
    }
  }
  smokersmatrix<-matrix(NA,120,13)
  for(k in 1:120)
  {
    for(s in 1:13)
    {
      smokersmatrix[k,s]<-smokers[k,s]*fsmokers[k,s]
    }
  }
  sum(smokersmatrix[1,1:13])
  smokerprofile<-vector(mode='numeric',120)
  smokerVA<-vector(mode='numeric',120)
  for(s in 1:120)
  {
    smokerprofile[s]<-sum(smokersmatrix[s,])
    smokerVA[s]<-((1/13)^2)*sum((smokerunc[s,])^2)
  } 
  
  
  fdiesel1<-vector(mode="numeric",15)
  fdiesel<-matrix(NA,120,15)
  for(d in 1:15)
  {
    fdiesel1[d]<-rnorm(1,2.8125,sd=2.8125*cv)
  }
  for(k in 1:120)
  {
    for(d in 1:15)
    {                     ###All fdiesel(ji) are drawn from a normal(2.8125,var=.3164) distribution
      fdiesel[k,d]<-fdiesel1[d]
    }
  }
  dieselmatrix<-matrix(NA,120,15)
  for(k in 1:120)
  {
    for(d in 1:15)
    {
      dieselmatrix[k,d]<-diesel[k,d]*fdiesel[k,d]
    }
  }
  dieselprofile<-vector(mode='numeric',120)
  dieselVA<-vector(mode='numeric',120)
  for(d in 1:120)
  {
    dieselprofile[d]<-sum(dieselmatrix[d,])
    dieselVA[d]<-((1/15)^2)*sum((dieselunc[d,])^2)
  }
  
  yfleet1VA<-matrix(NA,120,3)
  yfleet1VA[,1]<-gasVA
  yfleet1VA[,2]<-smokerVA
  yfleet1VA[,3]<-dieselVA
  yfleet1<-gasprofile+smokerprofile+dieselprofile
  q=1
  while(q <= length(yfleet1))
  {                   
    if (yfleet1[q]==0)
    {
      yfleet1VA<-yfleet1VA[-q,];   
      yfleet1<-yfleet1[-q];
      q=q-1
    }
    q=q+1
  }
  ysims1<-matrix(NA,93,1)
  for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
  {
    ysims1[i,1]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cvlog)^2)+1)),sd=sqrt(log(((cvlog)^2)+1))))
  }
  fpre1<-matrix(NA,1,3)
  fpre1[1,1]<-mean(fgas[1,])
  fpre1[1,2]<-mean(fsmokers[1,])
  fpre1[1,3]<-mean(fdiesel[1,])
  
  out<-NULL
  out[[1]]<-ysims1
  out[[2]]<-fpre1
  return(out)
}

#datagen(cv)



days<-50
ysims<-matrix(NA,93,days)
fpre<-matrix(NA,days,3)
cv<-.01
cvlog<-.01
fleet1profilesarr<-array(NA,dim=c(93,3,days))
for(i in 1:days)
{
  pgenerated<-profilegen()
  fleet1profilesarr[,,days]<-pgenerated[[1]]
  dataday<-datagen(cv,cvlog,pgenerated[[2]],pgenerated[[3]],pgenerated[[4]])
  fpre[i,]<-dataday[[2]][1,]
  ysims[,i]<-dataday[[1]][,1]
}
fpre
ysims