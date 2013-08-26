setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\PMFstuff\\pmf_r_interface_old")
source("jeff's functions.R")
#load(".Rdata")



pmf3<-function(cv)
{
  days<-1
  VC<-cv*ysimspmf
  VCt<-t(VC)
  ysimst<-t(ysimspmf)
  write.table(ysimst,file="ysim1dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
  write.table(VCt,file="yVC1dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
 
  n<-1
  
  uncertainty<-t(yfleet1VA1)
  for(i in 1:3)
  {
    for(j in 1:93)
    {
      if (uncertainty[i,j]==0)
      {
        uncertainty[i,j]=.00001
      }
    }
  }
  
  for(i in 1:93)
  {
    for(j in 1:3)
    {
      if (yfleet1profiles1[i,j]==0)
      {
        yfleet1profiles1[i,j]=.0001
      }
    }
  }
  yfleet1profiles1[34,1]<-.0005
  #for(i in 1:93)
  #{
  #if(yfleet1profiles[i,3]==0)
  #{
  # yfleet1profiles[i,3]=.00001
  # }
  #}
  i<-min(dim(yfleet1profiles1))
  k<-max(dim(yfleet1profiles1))
  n<-1
  
  #runpmfplain("default_nogkey.ini",ysimst,VCt,mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
  #mypaste("testq",i,".txt"),fpeak= 0.1,i,nrow(ysims),ncol(ysims))    
  #   if no gkeying is used, change "keydat" to "sldat"
  
  
  #lam<-read.table(mypaste("testlam",i,".txt"))
  #fac<-read.table(mypaste("testfac",i,".txt"))
  
  
  gkeyprep(yfleet1profiles1,uncertainty,"yVC1dat.txt","ysim1dat.txt")     #gets ready to gkey
  
  runpmf("default_gkey_3source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
         mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
  #   if no gkeying is used, change "keydat" to "sldat"
  #runpmf("default_gkey_3source.ini","ysim50dat.txt",mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
  # mypaste("testq",i,".txt"),i,n+i,k)      
  
  #runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
  #mypaste("testq",i,".txt"),fpeak= -0.1,i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"
  
  lam<-read.table(mypaste("testlam",i,".txt"))
  fac<-read.table(mypaste("testfac",i,".txt"))
  fac<-fac[-c(1:i),]    #if gkeying is used
  ##lamplots(lam,"gkey")  #the imagemagick batch files use these names, so if they are changed the .bat files must be changed
  ##facplots(fac,"gkeyf","2001/6/1","2003/6/6")
  ##weekendplots(i,"gkeyw")
  fac<-as.matrix(fac)
  apply(fac/83.59759,2,mean)
  apply(fac/135.25885,2,mean)
  apply(fac/932.02049,2,mean)
  fac1<-fac
  fac1[,1]<-fac[,1]/83.59759
  fac1[,2]<-fac[,2]/135.25885
  fac1[,3]<-fac[,3]/932.02049
  ###for cv=.5
  #fac1[,1]<-fac[,2]/83.59759
  #fac1[,2]<-fac[,1]/135.25885
  #fac1[,3]<-fac[,3]/932.02049
  ###
  facmse<-matrix(NA,days,3)
  for(i in 1:days)
  {
    facmse[i,1]<-(fac1[i,1]-fpre1[1])^2
    facmse[i,2]<-(fac1[i,2]-fpre1[2])^2
    facmse[i,3]<-(fac1[i,3]-fpre1[3])^2
  }
  
  out<-NULL
  out[[1]]<-fac1
  out[[2]]<-facmse
  return(out)
}


iter<-100
cv<-.05
#iter<-10
fac1<-matrix(NA,iter,3)
for(t in 1:iter)
{
  fpre1<-fpre[t,]
  yfleet1profiles1<-fleet1profilesarr[,,t]
  yfleet1VA1<-fleet1VAarr[,,t]
  ysimspmf<-ysims[,t]
  
  pmf3run<-pmf3(cv)
  fac1<-pmf3run[[1]]
  facmse<-pmf3run[[2]] 
}

facmeans<-apply(fac1,2,mean)
facmsemeans<-apply(facmse,2,mean)
pmfsim<-matrix(NA,1,6)
colnames(pmfsim)<-c("fgas","smoker","fdiesel","msegas","msesmoker","msediesel")

for(i in 1:3)
{
  pmfsim[1,i]<-facmeans[[i]]
} 
pmfsim[,4]<-facmsemeans[1]
pmfsim[,5]<-facmsemeans[2]
pmfsim[,6]<-facmsemeans[3]
pmfsim

lam<-as.matrix(lam)
ysims/fac%*%t(lam)
yfleet1profiles/lam




fpre1<-fpre[9,]
yfleet1profiles1<-fleet1profilesarr[,,9]
yfleet1VA1<-fleet1VAarr[,,9]
ysimspmf<-ysims[,9]

yfleet1profiles1[20,3]<-.00001

pmf3run<-pmf3(cv)
fac1<-pmf3run[[1]]
facmse<-pmf3run[[2]] 
