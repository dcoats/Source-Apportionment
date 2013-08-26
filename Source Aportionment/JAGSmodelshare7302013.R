###########################################################################################
######################### Simulations with priors for lambda  #############################
###########################################################################################


                                ########################
                                ###### Preparation #####
                                ########################
library(R2jags)
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")
days<-50     ###number of days n that will be simulated
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

#### Data forJAGS model
dim(ysims)
y<-ysims   ### y is the data
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles  ### matrix with p species and k profiles
profiles<-3 #number of profiles
cv<-.2

                             ####These 2 functions solve for alpha and beta given the desired mean and variance.
                             ####Use these functions to get alpha and beta when you want informative priors for lambda
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

                            #######Creating informative gamma priors for every species of the 3 profiles
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
  }
}

                    ####These priors were causing problems with Jags so they were given generic values
alpha[11,2]<-.1
beta[11,2]<-1000
alpha[30,2]<-.01
alpha[52,2]<-.1
beta[52,2]<-1000
alpha[71,2]<-.1
beta[71,2]<-1000
alpha
beta

                    ###########################################
                    #### Script of model to be run in JAGS ####
                    ###########################################

### arr[p,t,k] is an array of the values resulting from matrix multiplication before summing each column
### multiplied[p,t] is the result of lambda[p,k]*f[k,t]
### In this model there is just a single flat gamma prior for F, but you could have different priors for every source

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

data<-c('logy','days','profiles','cv','alpha','beta')  #### data must include all variables of actual data that is to be used in the simulations
parameters<-c('f')   #### parameters are all parameters you wish to observe as far as results 
bjags2.sim<-jags(data,inits=NULL,parameters,model.file='mvn.txt',
                 n.iter=10000,n.burnin=2000,
                 n.chains=1,n.thin=1)                           ###jags() runs jags with options for number of iterations, burnins, chains, and thinning

traceplot(bjags2.sim)   ### Traceplots for all parameters listed in "parameters"
bjags2.sim           #### Shows summary with means and sd for all parameters and deviance as well as DIC information
bjags2.out<-as.mcmc(bjags2.sim)  ###bjags2.out is a matrix with the values of all the iterations for all the parameters
variable.names(bjags2.out)   ### Shows the order of the parameters
head(bjags2.out[,2])
summary(bjags2.out)[[1]]  ###basically gives same information as bjags2.sim, has means and sd's for all the parameters and deviance

bjags2.out[1,]   ####shows 1st iteration for all of the parameters


###################################################################################################################
################### Simulations  of multiple iterations of 1 day ##################################################
###################################################################################################################

####The following code uses the same model except it runs jags multiple times to simulate n iterations instead of using n days in the data

library(R2jags)
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")

     ###number of days that will be simulated
days<-1
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

y<-ysims
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles
profiles<-3 #number of profiles
cv<-.2


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
  }
}
alpha[11,2]<-.1
beta[11,2]<-1000
alpha[30,2]<-.01
alpha[52,2]<-.1
beta[52,2]<-1000
alpha[71,2]<-.1
beta[71,2]<-1000
alpha
beta




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

iter<-100      #### of times to run JAGS with one day
sim1<-matrix(NA,iter,9)
sim1means<-matrix(NA,1,9)
colnames(sim1)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
colnames(sim1means)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")



                 #### The following loop runs jags iter times and puts the means for F and the mse's for every iteration in the matrix sim1
for(i in 1:iter)
{
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
  sim1[i,7]<-(sim1[i,1]-15)^2
  sim1[i,8]<-(sim1[i,2]-.937)^2
  sim1[i,9]<-(sim1[i,3]-2.8215)^2
}
                              ### sim1means is a row of values for the three estimated f's, their std deviations, and mse's
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