##############################
###### gas vehicle data ######
##############################
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
SI<-read.csv("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\AVERAGINGSIweighted.csv",skip=2,header=T)
head(SI)
SI<-SI[-(1:3),]

SIunc<-SI
for(i in 1:dim(SI)[1])      ### forming data set SIunc containing only uncertainties
{
  if (SIunc$X[i]=="unc"||SIunc$X[i]=="UNC"||SIunc$X[i]=="")
  {
    SIunc<-SIunc[-(i-1),];
  }
  return
}
head(SIunc)
SIunc<-SIunc[-(121:188),]

for(i in 1:dim(SI)[1])      ### deleting rows containing uncertainty from source profiles data frame
{
  if (SI$X[i]=="unc"||SI$X[i]=="UNC"||SI$X[i]=="")
  {
    SI<-SI[-i,];
    i=i+1;
  }
  return
}

###SI$X<-NULL
head(SI)
SIwBadSmokers<-SI[,1:18]
head(SIwBadSmokers)
SInoBadSmokers<-SI[,21:37]
head(SInoBadSmokers)
SIonlyGoodSmokers<-SI[,38:44]
head(SIonlyGoodSmokers)
SIonlyBadSmokers<-SI[,45:48]
head(SIonlyBadSmokers)
SInewDefinedSmokers<-SI[,49:50]
head(SInewDefinedSmokers)
SIbyYear<-SI[,52:65]
head(SIbyYear)
SI<-SI[,1:65]
SIunc<-SIunc[,1:65]
#SI<-SI[-(1:4),]
SI<-SI[-(121:184),]
SI$X.1<-NULL
SI$X.2<-NULL
head(SI)     ###trimmed data frame with all estimated weighted profiles for gas vehicles
#save(SI,file="SI") 
#save(SIunc,file="SIunc")

###Diesel vehicle data###

CI<-read.csv("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\AVERAGINGCIweighted.csv",skip=2,header=T)
CIunc<-CI
for(i in 1:dim(CIunc)[1])      ### forming data set CIunc containing only uncertainties
{
  if (CIunc$X[i]=="unc"||CIunc$X[i]=="UNC"||CIunc$X[i]=="")
  {
    CIunc<-CIunc[-(i-1),];
  }
  return
}
CIunc<-CIunc[1:22]
head(CIunc)
CIunc<-CIunc[-(121:188),]
for(i in 1:dim(CI)[1])      ### deleting rows containing uncertainty from source profiles data frame
{
  if (CI$X[i]=="unc"||CI$X[i]=="UNC"||CI$X[i]=="")
  {
    CI<-CI[-i,];
    i=i+1;
  }
  return
}
###CI$X<-NULL
CI<-CI[,1:22]
#CI<-CI[-(1:4),]
head(CI)        ###trimmed data frame with all estimated weighted profiles for diesel vehicles
CIunc<-CIunc[,1:22]
#save(CI,file="CI")
#save(CIunc,file="CIunc")

#########################################################################################################################################
###  Fleet 1 for BOTH WARM AND COLD STARTS with 80% gas vehicles, 5% smokers(including bad smokers), and 15% of vehicle being diesel  ###
#########################################################################################################################################




cv<-.5
gas<-SI[,4:7]
gas
gasunc<-SIunc[,4:7]
is.matrix(gas)
gas<-gas[-(121:124),]
for (i in 1:4)               ###converting profiles from factors to numeric data
{
  gas[,i]<-as.numeric(as.character(gas[,i]))
  gasunc[,i]<-as.numeric(as.character(gasunc[,i]))
  gas[,i]<-gas[,i]/4
}
gasp<-matrix(NA,120,1)
for(i in 1:120)
{
gasp[i]<-sum(gas[i,])
}
fgas1<-vector(mode="numeric",4)
fgas<-matrix(NA,120,4)
for(i in 1:4)
{
  fgas1[i]<-rnorm(1,15,sd=15*cv)
}
for(j in 1:120)
{
  for(i in 1:4)
  {                     ###All fgas(ji) are drawn from a normal(15,var=9) distribution
    fgas[j,i]<-fgas1[i]
  }
}
fgas
gasmatrix<-matrix(NA,120,4)
for(j in 1:120)
{
  for(i in 1:4)
  {
    gasmatrix[j,i]<-gas[j,i]*fgas[j,i]
  }
}
gasmatrix
gasprofile<-vector(mode='numeric',120)
gasVA<-vector(mode='numeric',120)
for(i in 1:120)
{
  gasprofile[i]<-sum(gasmatrix[i,])
  gasVA[i]<-(1/16)*sum((gasunc[i,])^2)
}
gasprofile
gasVA


smokers<-SI[,38:50]     ###creating profile for smoker cars
smokers
smokerunc<-SIunc[,38:50]
is.matrix(smokers)
for (i in 1:13)               ###converting profiles from factors to numeric data
{
  smokers[,i]<-as.numeric(as.character(smokers[,i]))
  smokerunc[,i]<-as.numeric(as.character(smokerunc[,i]))
  smokers[,i]<-smokers[,i]/13
}
smokerp<-matrix(,120,1)
for(i in 1:120)
{
  smokerp[i]<-sum(smokers[i,])
}
fsmokers1<-vector(mode="numeric",13)
fsmokers<-matrix(NA,120,13)
for(i in 1:13)
{
  fsmokers1[i]<-rnorm(1,.937,sd=.937*cv)
}
for(j in 1:120)
{
  for(i in 1:13)
  {                     ###All fsmokers(ji) are drawn from a normal(.937,var=.03515625) distribution
    fsmokers[j,i]<-fsmokers1[i]
  }
}
fsmokers
smokersmatrix<-matrix(NA,120,13)
for(j in 1:120)
{
  for(i in 1:13)
  {
    smokersmatrix[j,i]<-smokers[j,i]*fsmokers[j,i]
  }
}
smokersmatrix
sum(smokersmatrix[1,1:13])
smokerprofile<-vector(mode='numeric',120)
smokerVA<-vector(mode='numeric',120)
for(i in 1:120)
{
  smokerprofile[i]<-sum(smokersmatrix[i,])
  smokerVA[i]<-((1/13)^2)*sum((smokerunc[i,])^2)
}
smokerprofile
smokerVA


diesel<-CI[,8:22]     ###creating profile for diesel vehicles
diesel
dieselunc<-CIunc[,8:22]
is.matrix(diesel)
for (i in 1:15)               ###converting profiles from factors to numeric data
{
  diesel[,i]<-as.numeric(as.character(diesel[,i]))
  dieselunc[,i]<-as.numeric(as.character(dieselunc[,i]))
  diesel[,i]<-diesel[,i]/15
}
dieselp<-matrix(,120,1)
for(i in 1:120)
{
  dieselp[i]<-sum(diesel[i,])
}
fdiesel1<-vector(mode="numeric",15)
fdiesel<-matrix(NA,120,15)
for(i in 1:15)
{
  fdiesel1[i]<-rnorm(1,2.8125,sd=2.8125*cv)
}
for(j in 1:120)
{
  for(i in 1:15)
  {                     ###All fdiesel(ji) are drawn from a normal(2.8125,var=.3164) distribution
    fdiesel[j,i]<-fdiesel1[i]
  }
}
fdiesel
dieselmatrix<-matrix(NA,120,15)
for(j in 1:120)
{
  for(i in 1:15)
  {
    dieselmatrix[j,i]<-diesel[j,i]*fdiesel[j,i]
  }
}
dieselmatrix
dieselprofile<-vector(mode='numeric',120)
dieselVA<-vector(mode='numeric',120)
for(i in 1:120)
{
  dieselprofile[i]<-sum(dieselmatrix[i,])
  dieselVA[i]<-((1/15)^2)*sum((dieselunc[i,])^2)
}
dieselprofile
dieselVA

yfleet1profiles<-matrix(NA,120,3)
yfleet1VA<-matrix(NA,120,3)
yfleet1profiles[,1]<-gasp
yfleet1profiles[,2]<-smokerp
yfleet1profiles[,3]<-dieselp
yfleet1VA[,1]<-gasVA
yfleet1VA[,2]<-smokerVA
yfleet1VA[,3]<-dieselVA
yfleet1profiles
yfleet1VA
yfleet1<-gasprofile+smokerprofile+dieselprofile
species<-as.matrix(SI$X)
species<-species[-(121:124),]
yfleet1
for(j in 1:5)
{
for (i in 1:116)    ####For some reason I usually have to run this a couple times for it to work, and running it in another for loop doesn't fix it,
{                   ####Run until there are 93 species
  if (yfleet1[i]==0)
  {
    species<-species[-i];
    yfleet1profiles<-yfleet1profiles[-i,];
    yfleet1VA<-yfleet1VA[-i,];   
    yfleet1<-yfleet1[-i];
  }   
}
}
yfleet1
species
yfleet1profiles
yfleet1VA


####Setting up for EV simulations###

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
#save(yfleet1profiles,file="yfleet1profiles")
#save(yfleet1VA,file="yfleet1VA")
#save(yfleet1,file="yfleet1")
load("yfleet1profiles")
load("yfleet1")

days<-50     ###number of days that will be simulated
ysims<-matrix(NA,93,days)
cv<-.5
  for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
  {
    for(j in 1:days)
    {
      ysims[i,j]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
    }
  }

dim(ysims)



#####################################################################
################### With 2 profiles instead of 3 ####################
                   ##############################
yfleet1profiles2<-matrix(NA,120,2)
yfleet1VA2<-matrix(NA,120,2)
yfleet1profiles2[,1]<-gasp
yfleet1profiles2[,2]<-smokerp
yfleet1VA2[,1]<-gasVA
yfleet1VA2[,2]<-smokerVA
yfleet1profiles2
yfleet1VA2
yfleet12<-gasprofile+smokerprofile
species<-as.matrix(SI$X)
species<-species[-(121:124),]
yfleet12
for(j in 1:5)
{
  for (i in 1:120)    ####For some reason I usually have to run this a couple times for it to work, and running it in another for loop doesn't fix it,
  {                   ####Run until there are 93 species
    if (yfleet12[i]==0)
    {
      species<-species[-i];
      yfleet1profiles2<-yfleet1profiles2[-i,];
      yfleet1VA2<-yfleet1VA2[-i,];   
      yfleet12<-yfleet12[-i];
    }   
  }
}
yfleet12
species
yfleet1profiles2
yfleet1VA2


####Setting up for EV simulations###

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
#save(yfleet1profiles2,file="yfleet1profiles2")
#save(yfleet1VA2,file="yfleet1VA2")
#save(yfleet12,file="yfleet12")
load("yfleet1profiles2")
load("yfleet1VA2")
load("yfleet12")

days<-50     ###number of days that will be simulated
ysims2<-matrix(NA,81,days)
cv<-.2
for (i in 1:81)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysims2[i,j]<-exp(rnorm(1,mean=log(yfleet12[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}

dim(ysims2)






#########################################################################################################################################
###  Fleet Cold for ONLY COLD STARTS with 80% gas vehicles, 5% smokers(including bad smokers), and 15% of vehicle being diesel  ###
#########################################################################################################################################

gascold<-SI[,c(4,6)]
gascold
is.matrix(gascold)
for (i in 1:2)               ###converting profiles from factors to numeric data
{
  gascold[,i]<-as.numeric(as.character(gascold[,i]))
}

fgascold1<-vector(mode="numeric",2)
fgascold<-matrix(NA,120,2)
for(i in 1:2)
{
  fgascold1[i]<-rnorm(1,15,sd=3)
}
for(j in 1:120)
{
  for(i in 1:2)
  {                     ###All fgas(ji) are drawn from a normal(15,var=9) distribution
    fgascold[j,i]<-fgascold1[i]
  }
}
fgascold
gascoldmatrix<-matrix(NA,120,2)
for(j in 1:120)
{
  for(i in 1:2)
  {
    gascoldmatrix[j,i]<-gascold[j,i]*fgascold[j,i]
  }
}
gascoldmatrix
sum(gascoldmatrix[1,1:2])
gascoldprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  gascoldprofile[i]<-sum(gascoldmatrix[i,])
}
gascoldprofile



smokersc<-SI[,c(38,40,42,44,46,48,50)]     ###creating profile for smoker cars
smokersc
is.matrix(smokersc)
for (i in 1:7)               ###converting profiles from factors to numeric data
{
  smokersc[,i]<-as.numeric(as.character(smokersc[,i]))
}

fsmokersc1<-vector(mode="numeric",7)
fsmokersc<-matrix(NA,120,7)
for(i in 1:7)
{
  fsmokersc1[i]<-rnorm(1,.937,sd=.1875)
}
for(j in 1:120)
{
  for(i in 1:7)
  {                     ###All fsmokers(ji) are drawn from a normal(.937,var=.03515625) distribution
    fsmokersc[j,i]<-fsmokersc1[i]
  }
}
fsmokersc
smokerscmatrix<-matrix(NA,120,7)
for(j in 1:120)
{
  for(i in 1:7)
  {
    smokerscmatrix[j,i]<-smokersc[j,i]*fsmokersc[j,i]
  }
}
smokerscmatrix
sum(smokerscmatrix[1,1:7])
smokercprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  smokercprofile[i]<-sum(smokerscmatrix[i,])
}
smokercprofile



dieselc<-CI[,16]     ###creating profile for diesel vehicles
dieselc
is.matrix(dieselc)
for (i in 1:1)               ###converting profiles from factors to numeric data
{
  dieselc[,i]<-as.numeric(as.character(dieselc[,i]))
}

fdieselc1<-vector(mode="numeric",1)
fdieselc<-matrix(NA,120,1)
for(i in 1:1)
{
  fdieselc1[i]<-rnorm(1,2.8125,sd=.5625)
}
for(j in 1:120)
{
  for(i in 1:1)
  {                     ###All fdiesel(ji) are drawn from a normal(2.8125,var=.3164) distribution
    fdieselc[j,i]<-fdieselc1[i]
  }
}
fdieselc
dieselcmatrix<-matrix(NA,120,1)
for(j in 1:120)
{
    dieselcmatrix[j]<-dieselc[j]*fdieselc[j]
}
dieselcmatrix
dieselcprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  dieselcprofile[i]<-sum(dieselcmatrix[i,])
}
dieselcprofile


yfleetcold<-gascoldprofile+smokercprofile+dieselcprofile
yfleetcold




ysimscold<-matrix(NA,120,days)
cv<-.2
for (i in 1:120)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysimscold[i,j]<-exp(rnorm(1,mean=log(yfleetcold[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}

ysimscold




#########################################################################################################################################
###  Fleet Warm for ONLY WARM STARTS with 80% gas vehicles, 5% smokers(including bad smokers), and 15% of vehicle being diesel  ###
#########################################################################################################################################


gasw<-SI[,c(5,7)]
gasw
is.matrix(gasw)
for (i in 1:2)               ###converting profiles from factors to numeric data
{
  gasw[,i]<-as.numeric(as.character(gasw[,i]))
}

fgasw1<-vector(mode="numeric",2)
fgasw<-matrix(NA,120,2)
for(i in 1:2)
{
  fgasw1[i]<-rnorm(1,15,sd=3)
}
for(j in 1:120)
{
  for(i in 1:2)
  {                     ###All fgas(ji) are drawn from a normal(15,var=9) distribution
    fgasw[j,i]<-fgasw1[i]
  }
}
fgasw
gaswmatrix<-matrix(NA,120,2)
for(j in 1:120)
{
  for(i in 1:2)
  {
    gaswmatrix[j,i]<-gasw[j,i]*fgasw[j,i]
  }
}
gaswmatrix
sum(gaswmatrix[1,1:2])
gaswprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  gaswprofile[i]<-sum(gaswmatrix[i,])
}
gaswprofile



smokersw<-SI[,c(39,41,43,45,47,49)]     ###creating profile for smoker cars
smokersw
is.matrix(smokersw)
for (i in 1:6)               ###converting profiles from factors to numeric data
{
  smokersw[,i]<-as.numeric(as.character(smokersw[,i]))
}

fsmokersw1<-vector(mode="numeric",6)
fsmokersw<-matrix(NA,120,6)
for(i in 1:6)
{
  fsmokersw1[i]<-rnorm(1,.937,sd=.1875)
}
for(j in 1:120)
{
  for(i in 1:6)
  {                     ###All fsmokers(ji) are drawn from a normal(.937,var=.03515625) distribution
    fsmokersw[j,i]<-fsmokersw1[i]
  }
}
fsmokersw
smokerswmatrix<-matrix(NA,120,6)
for(j in 1:120)
{
  for(i in 1:6)
  {
    smokerswmatrix[j,i]<-smokersw[j,i]*fsmokersw[j,i]
  }
}
smokerswmatrix
sum(smokerswmatrix[1,1:6])
smokerwprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  smokerwprofile[i]<-sum(smokerswmatrix[i,])
}
smokerwprofile



dieselw<-CI[,14]     ###creating profile for diesel vehicles
dieselw<-as.matrix(dieselw)
for (i in 1:1)               ###converting profiles from factors to numeric data
{
  dieselw[,i]<-as.numeric(as.character(dieselw[,i]))
}

fdieselw1<-vector(mode="numeric",1)
fdieselw<-matrix(NA,120,1)
for(i in 1:1)
{
  fdieselw1[i]<-rnorm(1,2.8125,sd=.5625)
}
for(j in 1:120)
{
  for(i in 1:1)
  {                     ###All diesels(ji) are drawn from a normal(2.8125,var=.3164) distribution
    fdieselw[j,i]<-fdieselw1[i]
  }
}
fdieselw
dieselwmatrix<-matrix(NA,120,1)
for(j in 1:120)
{
  dieselwmatrix[j]<-dieselw[j]*fdieselw[j]
}
dieselwmatrix
dieselwprofile<-vector(mode='numeric',120)
for(i in 1:120)
{
  dieselwprofile[i]<-sum(dieselwmatrix[i,])
}
dieselwprofile


yfleetw<-gaswprofile+smokerwprofile+dieselwprofile
yfleetw




ysimswarm<-matrix(NA,120,days)
cv<-.2
for (i in 1:120)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysimswarm[i,j]<-exp(rnorm(1,mean=log(yfleetw[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}

ysimswarm



