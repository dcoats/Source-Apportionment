setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("SI")
load("SIunc")
load("CI")
load("CIunc")
library(gtools)  ###This library contains functions for Dirichlet Distribution
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment")
source("profiles.R")  ### runs code to create single profile and VA
#singleprofile


datagen<-function(cv,cvlog)
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
  dieseldir<-rdirichlet(1,c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
  for (i in 1:15)               ###converting profiles from factors to numeric data
  {
    diesel[,i]<-as.numeric(as.character(diesel[,i]))
    dieselunc[,i]<-as.numeric(as.character(dieselunc[,i]))
    diesel[,i]<-diesel[,i]*dieseldir[1,i]
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
  
  #### Producing Data###
  
  
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
  gasVAmat<-matrix(NA,4,4)
  for(i in 1:4)
  {
    for(j in 1:4)
    {
      if(i != j)
      {
        gasVAmat[i,j]<-gasdir[1,i]*gasdir[1,j]*cov(gas[,i],gas[,j])
      }
      else
      {
        gasVAmat[i,j]<-0
      }
    }
  }
  gasprofile<-vector(mode='numeric',120)
  gasVA<-vector(mode='numeric',120)
  for(k in 1:120)
  {
    gasprofile[k]<-sum(gasmatrix[k,])
    #gasVA[k]<-(1/16)*sum((gasunc[k,])^2)
    #gasVA[k]<-((gasdir[1,1]^2)*gasunc[k,1])+((gasdir[1,2]^2)*gasunc[k,2])+((gasdir[1,3]^2)*gasunc[k,3])+((gasdir[1,4]^2)*gasunc[k,4])+2*gasdir[1,1]*gasdir[1,2]*cov(gas[,1],gas[,2])+2*gasdir[1,1]*gasdir[1,3]*cov(gas[,1],gas[,3])+2*gasdir[1,1]*gasdir[1,4]*cov(gas[,1],gas[,4])+2*gasdir[1,2]*gasdir[1,3]*cov(gas[,2],gas[,3])+2*gasdir[1,2]*gasdir[1,4]*cov(gas[,2],gas[,4])+2*gasdir[1,3]*gasdir[1,4]*cov(gas[,3],gas[,4])
    gasVA[k]<-((gasdir[1,1]^2)*gasunc[k,1])+((gasdir[1,2]^2)*gasunc[k,2])+((gasdir[1,3]^2)*gasunc[k,3])+((gasdir[1,4]^2)*gasunc[k,4])+sum(gasVAmat)
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
  smokerVAmat<-matrix(NA,13,13)
  for(i in 1:13)
  {
    for(j in 1:13)
    {
      if(i != j)
      {
        smokerVAmat[i,j]<-smokdir[1,i]*smokdir[1,j]*cov(smokers[,i],smokers[,j])
      }
      else
      {
        smokerVAmat[i,j]<-0
      }
    }
  }
  smokerprofile<-vector(mode='numeric',120)
  smokerVA<-vector(mode='numeric',120)
  for(k in 1:120)
  {
    smokerprofile[k]<-sum(smokersmatrix[k,])
    #smokerVA[k]<-((1/13)^2)*sum((smokerunc[k,])^2)
    smokerVA[k]<-((smokdir[1,1]^2)*smokerunc[k,1]) + ((smokdir[1,2]^2)*smokerunc[k,2]) + ((smokdir[1,3]^2)*smokerunc[k,3]) + ((smokdir[1,4]^2)*smokerunc[k,4]) + ((smokdir[1,5]^2)*smokerunc[k,5]) + ((smokdir[1,6]^2)*smokerunc[k,6])
    + ((smokdir[1,7]^2)*smokerunc[k,7]) + ((smokdir[1,8]^2)*smokerunc[k,8]) + ((smokdir[1,9]^2)*smokerunc[k,9]) + ((smokdir[1,10]^2)*smokerunc[k,10]) + ((smokdir[1,11]^2)*smokerunc[k,11]) + ((smokdir[1,12]^2)*smokerunc[k,12]) + ((smokdir[1,13]^2)*smokerunc[k,13])
    + sum(smokerVAmat)
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
  
  dieselVAmat<-matrix(NA,15,15)
  for(i in 1:15)
  {
    for(j in 1:15)
    {
      if(i != j)
      {
        dieselVAmat[i,j]<-dieseldir[1,i]*dieseldir[1,j]*cov(diesel[,i],diesel[,j])
      }
      else
      {
        dieselVAmat[i,j]<-0
      }
    }
  }
  dieselsquares<-matrix(NA,120,15)
  for(i in 1:15)
  {
    for(j in 1:120)
    {
      dieselsquares[j,i]<-(dieseldir[1,i]^2)*dieselunc[j,i]
    }
  }
  dieselvar<-apply(dieselsquares,1,sum)
  dieselprofile<-vector(mode='numeric',120)
  dieselVA<-vector(mode='numeric',120)
  for(k in 1:120)
  {
    dieselprofile[k]<-sum(dieselmatrix[k,])
    #dieselVA[k]<-((1/15)^2)*sum((dieselunc[k,])^2)
  }
  dieselVA<-dieselvar+sum(dieselVAmat)
  
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
  out[[3]]<-yfleet1profiles
  out[[4]]<-yfleet1VA
  return(out)
}

#datagen(cv)



days<-50
ysims<-matrix(NA,93,days)
fpre<-matrix(NA,days,3)
cv<-.05
cvlog<-.05
fleet1profilesarr<-array(NA,dim=c(93,3,days))
fleet1VAarr<-array(NA,dim=c(93,3,days))
for(i in 1:days)
{
  dataday<-datagen(cv,cvlog)
  fleet1profilesarr[,,i]<-dataday[[3]]
  fleet1VAarr[,,i]<-dataday[[4]]
  fpre[i,]<-dataday[[2]][1,]
  ysims[,i]<-dataday[[1]][,1]
}
fpre
ysims
fleet1profilesarr
fleet1VAarr

singleprofile