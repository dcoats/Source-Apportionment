

######################################################
################    5 Sources   ######################
######################################################

lambda<-lam5
apply(fac5,2,mean)

i<-min(dim(lambda))
k<-max(dim(lambda))
n<-749


#testing
uncertainty<-array(NA,dim=c(i,k))
#uncertainty[1,]<-.0001         #copper smelter  
uncertainty[1,]<-.01           #fireworks 
uncertainty[2,]<-.01           #zinc smelter 
uncertainty[3,]<-.1           #summer 
uncertainty[4,]<-1           #winter
uncertainty[5,]<-1           #vehicle    
unc<-read.table("slunc.txt")
dim(unc)
gkeyprep(lambda,uncertainty,"slunc.txt","sldat.txt")     #gets ready to gkey

runpmf("default_gkey_5source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"

#runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       # mypaste("testq",i,".txt"),fpeak= -0.1,i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"

lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))

fac<-fac[-c(1:i),]    #if gkeying is used
apply(fac,2,mean)
sum(apply(fac,2,mean))



#######################################################
################    6 Sources   #######################
#######################################################
lams7<-read.table("lam7.txt")  #the one I was working with
smelter<-apply(lams7[,c(1,3)],1,mean)   #smelter

lambda<-cbind(lams7[,c(1:6)])

i<-min(dim(lambda))
k<-max(dim(lambda))
n<-749


#testing
uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-.0001         #copper smelter  
uncertainty[2,]<-.01           #fireworks 
uncertainty[3,]<-.01           #zinc smelter 
uncertainty[4,]<-.1           #summer 
uncertainty[5,]<-1           #winter
uncertainty[6,]<-1           #vehicle    
unc<-read.table("slunc.txt")
dim(unc)
gkeyprep(lambda,uncertainty,"slunc.txt","sldat.txt")     #gets ready to gkey

runpmf("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
       mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
#   if no gkeying is used, change "keydat" to "sldat"

runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
        mypaste("testq",i,".txt"),fpeak= -0.1,i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"

lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))

fac<-fac[-c(1:i),]    #if gkeying is used
apply(fac,2,mean)
sum(apply(fac,2,mean))