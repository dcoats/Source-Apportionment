library(mail)
library(R2WinBUGS)
library(arm)
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")

days<-10     ###number of days that will be simulated
ysims<-matrix(NA,93,days)
cv<-.2
for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysims[i,j]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}

dim(ysims)

############################################################################
############################# linear model with no priors on lambda ########
############################################################################
y<-ysims
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles
profiles<-3 #number of profiles
cv<-.2
days<-10
plot(density(rgamma(10000,shape=.001,rate=1)))
mdl<-"
model
{

for(p in 1:93)
{
for(t in 1:days)
{
for(k in 1:profiles)
{
arr[p,t,k]<-yfleet1profiles[p,k]*f[k,t]
}
}
}

for(p in 1:93)
{
for(t in 1:days)
{
multiplied[p,t]<-sum(arr[p,t,])
mu[p,t]<-log(multiplied[p,t])-.5*vari
logy[p,t] ~ dnorm(mu[p,t],prec)
}
}

for(k in 1:profiles)    
{
for(t in 1:days)  
{
f[k,t]~dgamma(2.5,.41666)    
}
}
vari<-log(((cv*cv)+1))
prec<-1/vari


}"
writeLines(mdl,'mvn.txt')

data<-c('logy','days','profiles','cv','yfleet1profiles')
parameters<-c('f')
mvn.sim<-bugs(data,inits=NULL,parameters,bugs.directory="C:/Users/David Coats/Downloads/winbugs14/WinBUGS14",model.file='mvn.txt',
                 n.iter=10000,n.burnin=2000,
                 n.chains=1,n.thin=1,debug=TRUE,codaPkg=T)
baysolution.out<-read.bugs(mvn.sim)
summary(baysolution.out)
plot(baysolution.out)

###############################################################

gammaparma<-function(mu,s2)
{
  beta<-mu/s2
  alpha<-mu*beta
  return(alpha)
}

gammaparmb<-function(mu,s2)
{
  beta<-mu/s2
  alpha<-mu*beta
  return(beta)
}

gammaparma(yfleet1profiles[7,2],yfleet1VA[7,2])
gammaparmb(yfleet1profiles[7,2],yfleet1VA[7,2])

alpha<-matrix(NA,93,3)
beta<-matrix(NA,93,3)
for(p in 1:93)
{
  for(k in 1:3)
  {
    alpha[p,k]<-gammaparma(yfleet1profiles[p,k],yfleet1VA[p,k])
    beta[p,k]<-gammaparmb(yfleet1profiles[p,k],yfleet1VA[p,k])
    if(is.na(alpha[p,k]))
    {
      alpha[p,k]=1
    }
    if(is.na(beta[p,k]))
    {
      beta[p,k]=1000
    }
    if(alpha[p,k]==Inf)
    {
      alpha[p,k]<-yfleet1profiles[p,k]
    }
    if(beta[p,k]==Inf)
    {
      beta[p,k]<-1
    }
  }
}
alpha
beta
alpha[46,2]
beta[46,2]
yfleet1profiles[48,2]
yfleet1VA[48,2]
                      ###WinBUGS didn't like it when the alphas were smaller than .1
                       ##This wasn't a problem in Jags, but for this model I had to modify quite a few of the parameters.
                       ##Small values for the shape parameter produce "cannot slice bracket" errors in WinBUGS
alpha[46,2]<-2.958
beta[46,2]<-100
alpha[50,2]<-1.617567
beta[50,2]<-100
alpha[48,2]<-2.116346
beta[48,2]<-100
alpha[7,1]<-2.555441
beta[7,1]<-61.75879
alpha[23,2]<-.9251652
beta[23,2]<-100
alpha[7,2]<-1
beta[7,2]<-1/.9115566
alpha[11,2]<-1
beta[11,2]<-1000
alpha[30,2]<-1
beta[30,2]<-1/.3910008
alpha[52,2]<-1
beta[52,2]<-1000
alpha[71,2]<-1
beta[71,2]<-1000
alpha[52,2]<-1
beta[52,2]<-1000
alpha[25,2]<-1
beta[25,2]<-1000
alpha[53,2]<-yfleet1profiles[53,2]*100
beta[53,2]<-100

library(R2WinBUGS)
library(arm)
y<-ysims
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles
profiles<-3 #number of profiles
cv<-.2
days<-10
plot(density(rgamma(10000,shape=.16,rate=10)))
mdl<-"
model
{

for(p in 1:93)
{
for(t in 1:days)
{
for(k in 1:profiles)
{
arr[p,t,k]<-lambda[p,k]*f[k,t]
}
}
}

for(p in 1:93)
{
for(t in 1:days)
{
multiplied[p,t]<-sum(arr[p,t,])
mu[p,t]<-log(multiplied[p,t])-.5*vari
logy[p,t] ~ dnorm(mu[p,t],prec)
}
}

for(p in 1:93)
{
  for(k in 1:profiles)
  {
    lambda[p,k]~dgamma(alpha[p,k],beta[p,k])
  }
}

for(k in 1:profiles)    
{
for(t in 1:days)  
{
f[k,t]~dgamma(2.5,.41666)    
}
}
vari<-log(((cv*cv)+1))
prec<-1/vari


}"
writeLines(mdl,'mvn.txt')

data<-c('logy','days','profiles','cv','alpha','beta')
parameters<-c('f')
mvn.sim<-bugs(data,inits=NULL,parameters,bugs.directory="C:/Users/David Coats/Downloads/winbugs14/WinBUGS14",model.file='mvn.txt',
                 n.iter=10000,n.burnin=2000,
                 n.chains=1,n.thin=1,debug=TRUE,codaPkg=T)
baysolution.out<-read.bugs(mvn.sim)
###sendmail("davidwilliamcoats@gmail.com","R notification","Simulations finished")
summary(baysolution.out)
autocorr.plot(baysolution.out)
baydat<-as.mcmc(baysolution.out)
variable.names(baydat)
f1sim<-matrix(NA,8000,1)
f2sim<-matrix(NA,8000,1)
f3sim<-matrix(NA,8000,1)
for (i in 1:8000)
{
  f1sim<-sum(baydat[i,2:51])/50
  f2sim<-sum(baydat[i,52:101])/50
  f3sim<-sum(baydat[i,102:151])/50
}
f1<-mean(f1sim)
f1
f2<-mean(f2sim)
f2
f3<-mean(f3sim)
f3