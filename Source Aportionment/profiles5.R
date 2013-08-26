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

singleprofile<-yfleet1profiles