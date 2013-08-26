##### Producing profiles with low volatility   (see profiles5high.R for high volatility)



setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("SI")
load("SIunc")
load("CI")
load("CIunc")
library(gtools)

############## Deciding appropriate Dirichlet distributions for differtnt levels of volatility

### gas
##### High Volatility
(1/4)^.4
mean(apply(rdirichlet(1000,c(.68,.68,.68,.68)),1,max))
##### Low volatility
(1/4)^.8
mean(apply(rdirichlet(1000,c(12,12,12,12)),1,max))

### smoker
##### High Volatility
(1/13)^.4
mean(apply(rdirichlet(1000,rep(.4,13)),1,max))
##### Low Volatility
(1/13)^.8
mean(apply(rdirichlet(1000,rep(8,13)),1,max))

### diesel
##### High Volatility
(1/15)^.4
mean(apply(rdirichlet(1000,rep(.38,15)),1,max))
#####Low Volatilitiy
(1/15)^.8
mean(apply(rdirichlet(1000,rep(9,15)),1,max))



#########################################################################################
################################# Creating Profiles #######################################
#########################################################################################



### Creating gas profiles ###
### Creating gas profiles ###
gas<-SI[,4:7]
gas
gasunc<-SIunc[,4:7]
is.matrix(gas)
gas<-gas[-(121:124),]
gasdir<-rdirichlet(1,c(12,12,12,12))
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
smokdir<-rdirichlet(1,rep(8,13))
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
dieseldir<-rdirichlet(1,rep(9,15))
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

gasVAmat<-matrix(NA,4,4)
for(i in 1:4)
{
  for(j in 1:4)
  {
    if(i != j)
    {
      gasVAmat[i,j]<-.25*.25*cov(gas[,i],gas[,j])
    }
    else
    {
      gasVAmat[i,j]<-0
    }
  }
}
gasVA<-vector(mode='numeric',120)
for(k in 1:120)
{
  gasVA[k]<-(((1/4)^2)*gasunc[k,1])+(((1/4)^2)*gasunc[k,2])+(((1/4)^2)*gasunc[k,3])+(((1/4)^2)*gasunc[k,4])+sum(gasVAmat)
}


smokerVAmat<-matrix(NA,13,13)
for(i in 1:13)
{
  for(j in 1:13)
  {
    if(i != j)
    {
      smokerVAmat[i,j]<-(1/13)*(1/13)*cov(smokers[,i],smokers[,j])
    }
    else
    {
      smokerVAmat[i,j]<-0
    }
  }
}

smokerVA<-vector(mode='numeric',120)
for(k in 1:120)
{
  smokerVA[k]<-(((1/13)^2)*smokerunc[k,1]) + (((1/13)^2)*smokerunc[k,2]) + (((1/13)^2)*smokerunc[k,3]) + (((1/13)^2)*smokerunc[k,4]) + (((1/13)^2)*smokerunc[k,5]) + (((1/13)^2)*smokerunc[k,6])
  + (((1/13)^2)*smokerunc[k,7]) + (((1/13)^2)*smokerunc[k,8]) + (((1/13)^2)*smokerunc[k,9]) + (((1/13)^2)*smokerunc[k,10]) + (((1/13)^2)*smokerunc[k,11]) + (((1/13)^2)*smokerunc[k,12]) + (((1/13)^2)*smokerunc[k,13])
  + sum(smokerVAmat)
} 


dieselVAmat<-matrix(NA,15,15)
for(i in 1:15)
{
  for(j in 1:15)
  {
    if(i != j)
    {
      dieselVAmat[i,j]<-(1/15)*(1/15)*cov(diesel[,i],diesel[,j])
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
    dieselsquares[j,i]<-((1/15)^2)*dieselunc[j,i]
  }
}
dieselvar<-apply(dieselsquares,1,sum)
dieselVA<-vector(mode='numeric',120)
dieselVA<-dieselvar+sum(dieselVAmat)

yfleet1VA<-matrix(NA,120,3)
yfleet1VA[,1]<-gasVA
yfleet1VA[,2]<-smokerVA
yfleet1VA[,3]<-dieselVA


q=1
while(q <= dim(yfleet1profiles)[1])
{                   
  if (sum(yfleet1profiles[q,])==0)
  {
    yfleet1profiles<-yfleet1profiles[-q,];
    #species<-species[-q,];
    yfleet1VA<-yfleet1VA[-q,]
    q=q-1;
  }
  q=q+1
}
singleprofile<-yfleet1profiles