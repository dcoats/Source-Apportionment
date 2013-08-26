setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\PMFstuff\\pmf_r_interface")
source("jeff's functions.R")
days<-500    ###number of days that will be simulated
ysims<-matrix(NA,93,days)
cv<-.000001
for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysims[i,j]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}





dim(ysims)
data1<-t(ysims)
dim(data1)
dim(yfleet1profiles)
lambda<-yfleet1profiles
unclam<-t(yfleet1VA)
for(i in 1:3)
{
  for(j in 1:93)
  {
    if (unclam[i,j]==0)
    {
      unclam[i,j]=.00001
    }
  }
}
unclam
uncsims<-(data1*.2)
write.table(data1,file="ysim50dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
write.table(uncsims,file="yunc50dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
i<-min(dim(lambda))
k<-max(dim(lambda))
n<-days

gkeyprep<-function(lambdakey,uncertainty,ambienterror="yVC50dat.txt",datamatrix="ysim50dat.txt")
{
  
  ## with this variant of gkeying I'm using the key as a starting value for lambda, and 
  ## not using a starting value for F (sorting it out would be a pain)
  i<-min(dim(lambdakey))
  meserror<-read.table(ambienterror)
  uncertainty<-rbind(uncertainty,meserror)
  write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)
  
  
  key<-t(lambdakey)
  write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)
  
  data<-read.table(datamatrix)
  keydata<-rbind(as.matrix(key),as.matrix(data))
  write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}

gkeyprep(lambda,unclam,ambienterror="yunc50dat.txt",datamatrix="ysim50dat.txt")     #gets ready to gkey

runpmf("default_gkey_3source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"

lam<-as.matrix(read.table(mypaste("testlam",i,".txt")))
fac<-as.matrix(read.table(mypaste("testfac",i,".txt")))
fac<-fac[-c(1:i),]  

fac
lam
mean(fac[,1]/fac[,2])
mean(fac[,3]/fac[,2])
fac%*%t(lam)/data1

runpmfplain<-function(inifilename,data,prec,F,Lambda,Q,fpeak=0,sources=6,rows=nrow(data),cols=ncol(data)){
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[9] <-paste("     ",fpeak)
  inifile[38] <-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[39] <-paste("    31   T \"OLD    \"  2000  \" ",prec,"                                       \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  #inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  
  mywrite.table(data,file="pmf_data.txt",quote=FALSE,na="-999.9")
  mywrite.table(prec,file="pmf_unc.txt",quote=FALSE,na="1")
  
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

runpmfplain("default_nogkey.ini",data1,uncsims,mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
            mypaste("testq",i,".txt"),fpeak= 0.1,i,nrow(data1),ncol(uncsims))  

################################################################################################
############################# With 2 profiles instead of 3 #####################################
                             ##############################

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles2")
load("yfleet1VA2")
load("yfleet12")
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\PMFstuff\\pmf_r_interface")

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
data2<-t(ysims2)
dim(data2)
dim(yfleet1profiles2)
lambda2<-yfleet1profiles2
unclam2<-t(yfleet1VA2)
for(i in 1:2)
{
  for(j in 1:81)
  {
    if (unclam2[i,j]==0)
    {
      unclam2[i,j]=.0000001
    }
  }
}
unclam2
uncsims2<-(data2*.2)
write.table(data2,file="ysim50dat2.txt",quote=F,row.names=FALSE,col.names=FALSE)
write.table(uncsims2,file="yunc50dat2.txt",quote=F,row.names=FALSE,col.names=FALSE)
i<-min(dim(lambda2))
k<-max(dim(lambda2))
n<-50
#
gkeyprep<-function(lambdakey,uncertainty,ambienterror="yVC50dat.txt",datamatrix="ysim50dat.txt")
{
  
  ## with this variant of gkeying I'm using the key as a starting value for lambda, and 
  ## not using a starting value for F (sorting it out would be a pain)
  i<-min(dim(lambdakey))
  meserror<-read.table(ambienterror)
  uncertainty<-rbind(uncertainty,meserror)
  write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)
  
  
  key<-t(lambdakey)
  write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)
  
  data<-read.table(datamatrix)
  keydata<-rbind(as.matrix(key),as.matrix(data))
  write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}
#

gkeyprep(lambda2,unclam2,ambienterror="yunc50dat2.txt",datamatrix="ysim50dat2.txt")     #gets ready to gkey

runpmf("default_gkey_2source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"

lam2<-as.matrix(read.table(mypaste("testlam",i,".txt")))
fac2<-as.matrix(read.table(mypaste("testfac",i,".txt")))
fac2<-fac2[-c(1:i),]  

fac2
lam2
mean(fac2[,1]/fac2[,2])
fac2%*%t(lam2)/data2

runpmfplain<-function(inifilename,data,prec,F,Lambda,Q,fpeak=0,sources=6,rows=nrow(data),cols=ncol(data)){
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[9] <-paste("     ",fpeak)
  inifile[38] <-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[39] <-paste("    31   T \"OLD    \"  2000  \" ",prec,"                                       \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  #inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  
  mywrite.table(data,file="pmf_data.txt",quote=FALSE,na="-999.9")
  mywrite.table(prec,file="pmf_unc.txt",quote=FALSE,na="1")
  
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

runpmfplain("default_nogkey.ini",data2,uncsims2,mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
            mypaste("testq",i,".txt"),fpeak= 0.1,i,nrow(data1),ncol(uncsims))  



###########################################################################
#######################  Scaled to add up to 1  ###########################
###########################################################################
ysimsprop<-matrix(NA,93,50)
VCprop<-matrix(NA,93,50)
for(i in 1:93)
{
  for(j in 1:50)
  {
    ysimsprop[i,j]<-ysims[i,j]/sum(ysims[,j])
    VCprop[i,j]<-.2*ysimsprop[i,j]
  }
}
VCpropt<-t(VCprop)
ysimsprop<-t(ysimsprop)
write.table(ysimsprop,file="ysim50datp.txt",quote=F,row.names=FALSE,col.names=FALSE)
write.table(VCpropt,file="yVC50datp.txt",quote=F,row.names=FALSE,col.names=FALSE)
lambda<-matrix(NA,93,3)
for(i in 1:93)
{
  for(j in 1:3)
  {
    lambda[i,j]<-yfleet1profiles[i,j]/sum(yfleet1profiles[,j])
  }
}
i<-min(dim(lambda))
k<-max(dim(lambda))
n<-50

uncertainty<-t(yfleet1VA)
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
uncertainty[,34]<-.0001


gkeyprep(lambda,uncertainty,"yVC50datp.txt","ysim50datp.txt")     #gets ready to gkey

runpmf("default_gkey_3source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"


runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
        mypaste("testq",i,".txt"),fpeak= -0.1,i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"

lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))
fac<-fac[-c(1:i),]    #if gkeying is used
##lamplots(lam,"gkey")  #the imagemagick batch files use these names, so if they are changed the .bat files must be changed
##facplots(fac,"gkeyf","2001/6/1","2003/6/6")
##weekendplots(i,"gkeyw")
fac<-as.matrix(fac)
lam<-as.matrix(lam)
data/fac%*%t(lam)




###############################################
######### Without 1st 5 rows ##################
###############################################

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
load("yfleet1profiles")
load("yfleet1VA")
load("yfleet1")
setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\PMFstuff\\pmf_r_interface")
source("jeff's functions.R")
days<-500    ###number of days that will be simulated
ysims<-matrix(NA,93,days)
cv<-.000001
for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysims[i,j]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}





dim(ysims)
data1<-t(ysims[-c(1:5),])
dim(data1)
dim(yfleet1profiles)
lambda<-yfleet1profiles[-c(1:5),]
unclam<-t(yfleet1VA[-c(1:5),])
for(i in 1:3)
{
  for(j in 1:88)
  {
    if (unclam[i,j]==0)
    {
      unclam[i,j]=.00001
    }
  }
}
unclam
uncsims<-(data1*.2)
write.table(data1,file="ysim50dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
write.table(uncsims,file="yunc50dat.txt",quote=F,row.names=FALSE,col.names=FALSE)
i<-min(dim(lambda))
k<-max(dim(lambda))
n<-days

gkeyprep<-function(lambdakey,uncertainty,ambienterror="yVC50dat.txt",datamatrix="ysim50dat.txt")
{
  
  ## with this variant of gkeying I'm using the key as a starting value for lambda, and 
  ## not using a starting value for F (sorting it out would be a pain)
  i<-min(dim(lambdakey))
  meserror<-read.table(ambienterror)
  uncertainty<-rbind(uncertainty,meserror)
  write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)
  
  
  key<-t(lambdakey)
  write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)
  
  data<-read.table(datamatrix)
  keydata<-rbind(as.matrix(key),as.matrix(data))
  write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}
unclam[3,31]<-.00001
unclam[3,81]<-.00001
unclam[3,84]<-.00001
unclam[2,6]<-.00001
gkeyprep(lambda,unclam,ambienterror="yunc50dat.txt",datamatrix="ysim50dat.txt")     #gets ready to gkey

runpmf("default_gkey_3source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"

lam<-as.matrix(read.table(mypaste("testlam",i,".txt")))
fac<-as.matrix(read.table(mypaste("testfac",i,".txt")))
fac<-fac[-c(1:i),]  

fac
lam
mean(fac[,1]/fac[,2])
mean(fac[,3]/fac[,2])
fac%*%t(lam)/data1
