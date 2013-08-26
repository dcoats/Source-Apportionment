
library(R2jags)
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")

days<-50     ###number of days that will be simulated
#days<-1
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
days<-50
plot(density(rgamma(10000,shape=2.5,rate=.41666)))
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
bjags.sim<-jags(data,inits=NULL,parameters,model.file='mvn.txt',
                 n.iter=10000,n.burnin=2000,
                 n.chains=1,n.thin=1)
#baysolution.out<-read.bugs(mvn.sim)
bjags.sim
bjags2.out<-as.mcmc(bjags.sim)
sim2<-matrix(NA,days,9)
sim2means<-matrix(NA,1,9)
colnames(sim2)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
colnames(sim2means)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
m<-summary(bjags2.out)[[1]][-1,]
#day<-row.names(summary(bjags2.out)[[1]])[-1]
#pos<-regexpr("[[(.*?),]",day)
#day<-substr(day,pos+1,pos+attr(pos,"match.length"))
#day<-as.numeric(as.character(day))
#day
c<-0
for(i in 1:days)
{
for(j in 1:profiles)
{
  c<-c+1
  sim2[i,j]<-m[c,1]
  sim2[i,j+3]<-m[c,2]
  sim2[i,7]<-(sim2[i,1]-15)^2
  sim2[i,8]<-(sim2[i,2]-.937)^2
  sim2[i,9]<-(sim2[i,3]-2.8215)^2
}
}
sim2means[1,1]<-mean(sim2[,1])
sim2means[1,2]<-mean(sim2[,2])
sim2means[1,3]<-mean(sim2[,3])
sim2means[1,4]<-mean(sim2[,4])
sim2means[1,5]<-mean(sim2[,5])
sim2means[1,6]<-mean(sim2[,6])
sim2means[1,7]<-mean(sim2[,7])
sim2means[1,8]<-mean(sim2[,8])
sim2means[1,9]<-mean(sim2[,9])

sim2means
###########################################################################################
######################### Simulations with priors for lambda  #############################
###########################################################################################
y<-ysims
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles
profiles<-3 #number of profiles
cv<-.2
days<-10

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
      alpha[p,k]=.1
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
    if((alpha[p,k]==0) & (beta[p,k]==0))
    {
      alpha[p,k]<-.1
      beta[p,k]<-1000
    }
    if(alpha[p,k]<.005)
    {
      alpha[p,k]<-.01
    }
    if(beta[p,k]<.005)
    {
      beta[p,k]<-.01
    }
  }
}



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
bjags2.sim<-jags(data,inits=NULL,parameters,model.file='mvn.txt',
                n.iter=10000,n.burnin=2000,
                n.chains=1,n.thin=1)

traceplot(bjags2.sim)

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Bayesian simulations")
load("jags100Kiter")   ###Model run with 100000 iterations,10000 burn in###
bjags2.sim
summary(bjags2.sim)
bjags2.out<-as.mcmc(bjags2.sim)
variable.names(bjags2.out)
head(bjags2.out[,2])
summary(bjags2.out)[[1]]

head(bjags2.out)   ####shows 1st iteration for all of the parameters


#############################################################################
################### Simulations  of multiple iterations of 1 day #############################################
############################################################################
##################
library(R2jags)
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")


y<-ysims
logy1<-log(y)
y<-as.matrix(y)
#yfleet1profiles
profiles<-3 #number of profiles
cv<-.05


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
      alpha[p,k]=.1
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
    if((alpha[p,k]==0) & (beta[p,k]==0))
    {
      alpha[p,k]<-.1
      beta[p,k]<-1000
    }
    if(alpha[p,k]<.005)
    {
      alpha[p,k]<-.01
    }
    if(beta[p,k]<.005)
    {
      beta[p,k]<-.01
    }
  }
}

days<-1

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

iter<-50
sim1<-matrix(NA,iter,9)
sim1means<-matrix(NA,1,9)
colnames(sim1)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
colnames(sim1means)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")

for(i in 1:iter)
{
  logy<-as.matrix(logy1[,iter])
  bjags.sim<-jags(data,inits=NULL,parameters,model.file='mvn.txt',
                   n.iter=10000,n.burnin=2000,
                   n.chains=1,n.thin=1)
  bjags.out<-as.mcmc(bjags.sim)
  sim1[i,1]<-summary(bjags.out)[[1]][2,1]
  sim1[i,2]<-summary(bjags.out)[[1]][3,1]
  sim1[i,3]<-summary(bjags.out)[[1]][4,1]
  sim1[i,4]<-summary(bjags.out)[[1]][2,2]
  sim1[i,5]<-summary(bjags.out)[[1]][3,2]
  sim1[i,6]<-summary(bjags.out)[[1]][4,2]
  sim1[i,7]<-(sim1[i,1]-fpre[i,1])^2
  sim1[i,8]<-(sim1[i,2]-fpre[i,2])^2
  sim1[i,9]<-(sim1[i,3]-fpre[i,3])^2
}

sim1means[1,1]<-mean(sim1[,1])
sim1means[1,2]<-mean(sim1[,2])
sim1means[1,3]<-mean(sim1[,3])
sim1means[1,4]<-mean(sim1[,4])
sim1means[1,5]<-mean(sim1[,5])
sim1means[1,6]<-mean(sim1[,6])
sim1means[1,7]<-mean(sim1[,7])
sim1means[1,8]<-mean(sim1[,8])
sim1means[1,9]<-mean(sim1[,9])

sim1means