#############################################################################
########### Multiple days with 1 profile for prior information ##############
#############################################################################

library(R2jags)

y<-ysims
logy<-log(y)
y<-as.matrix(y)
yfleet1profiles<-singleprofile
profiles<-3 #number of profiles
cv<-.05
days<-50

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

#traceplot(bjags2.sim)
summary(bjags2.sim)
bjags2.out<-as.mcmc(bjags2.sim)
variable.names(bjags2.out)
#summary(bjags2.out)[[1]]

sim2<-matrix(NA,days,9)
sim2means<-matrix(NA,1,9)
colnames(sim2)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
colnames(sim2means)<-c("fgas","fsmoker","fdiesel","sdgas","sdsmoker","sddiesel","msegas","msesmoker","msediesel")
m<-summary(bjags2.out)[[1]][-1,]
c<-0
for(i in 1:days)
{
  for(j in 1:profiles)
  {
    c<-c+1
    sim2[i,j]<-m[c,1]
    sim2[i,j+3]<-m[c,2]
    sim2[i,7]<-(sim2[i,1]-fpre[i,1])^2
    sim2[i,8]<-(sim2[i,2]-fpre[i,2])^2
    sim2[i,9]<-(sim2[i,3]-fpre[i,3])^2
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