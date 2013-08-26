## Try 3 through 10 sources
##
  write(t(slunc2),file="slunc2.txt",ncol=ncol(slunc2))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
datesnew <- getdates()
firstofmonth <- c(
  "2001-06-01","2001-07-01","2001-08-01","2001-09-01","2001-10-01","2001-11-01","2001-12-01",
  "2002-01-01","2002-02-01","2002-03-01","2002-04-01","2002-05-01","2002-06-01",
  "2002-07-01","2002-08-01","2002-09-01","2002-10-01","2002-11-01","2002-12-01",
  "2003-01-01","2003-02-01","2003-03-01","2003-04-01","2003-05-01")


#OLD facplots3 <- function(fac,name="NOgkeyf",datevec)
{
  dates<-datevec
  for(iii in 1:ncol(fac))
  {
    jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=175)
    par(mar=c(2,6,1,1))
    par(cex=1) 
    plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",
          ylab=mypaste(iii," of ",i," factors"))
    par(cex=.7)
    axis.Date(1,at=firstofmonth,
      format="%b-%y",las=2)
    dev.off()
  }
}

facplots3 <- function(fac,name="NOgkeyf",datevec,facnames)
{
  dates<-datevec
  for(iii in 1:ncol(fac))
  {
    jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=125)
    par(mar=c(2.75,4.5,.5,1))
    par(mgp=c(2,.5,0))
    par(cex=1) 
    plot(dates,fac[,iii],type="l",col="black",xlab="",xaxt="n",
          ylab="",las=2,tcl=-.2)
    par(cex=1.25)
    title(ylab=facnames[iii])
    par(cex=.6)
    par(mgp=c(1.5,.4,0))
    axis.Date(1,at=firstofmonth,
      format="%b%y",las=2,tcl=-.2)
    dev.off()
  }
}
#facplots3(fac,name="hey",newdates)
facplots3(fac9gkeynostart,"gkeyf",newdates,facnames)


facplots3ps <- function(fac,name="gkeyf",datevec,facnames)
{
  dates<-datevec
  for(iii in 1:ncol(fac))
  {
    pdf(file=mypaste(name,i,iii,".pdf") ,width=8.5, height=1.5)
    par(mar=c(2.75,4.5,.5,1))
    par(mgp=c(2,.5,0))
    par(cex=1) 
    plot(dates,fac[,iii],type="l",col="black",xlab="",xaxt="n",
          ylab="",las=2,tcl=-.2)
    par(cex=1.25)
    title(ylab=facnames[iii])
    par(cex=.6)
    par(mgp=c(1.5,.4,0))
    axis.Date(1,at=firstofmonth,
      format="%b%y",las=2,tcl=-.2)
    dev.off()
  }
}
#facplots3(fac,name="hey",newdates)
facplots3ps(fac9gkeynostart,"gkeyf",newdates,facnames)

date()

i <- 9
newdates <- datesnew[apply(sldat==-999.9,1,sum)==0]
mynames <- names(sldat2)

runpmfplain("default_nogkey.ini","sldat2.txt","slunc2.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
#runpmfplain("default_nogkeyLIMCHANGE.ini","sldat2.txt","slunc2.txt",mypaste("testfac",i,".txt"),
#   mypaste("testlam",i,".txt"),
#   mypaste("testq",i,".txt"),fpeak= 0,i,nrow(sldat2),ncol(sldat2))    
lam9 <- lam <- read.table(mypaste("testlam",i,".txt"))
fac9 <- fac <- read.table(mypaste("testfac",i,".txt"))
lamplots2old(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3old(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


i   <-    9
lam <- lam9
fac <- fac9
lamplots2old(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3old(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 



i   <-    8
lam <- lam8
fac <- fac8
lamplots2(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 



i   <-    7
lam <- lam7
fac <- fac7
lamplots2(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


i   <-    6
lam <- lam6
fac <- fac6
lamplots2(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


i   <-    5
lam <- lam5
fac <- fac5
lamplots2(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


i   <-    4
lam <- lam4
fac <- fac4
lamplots2(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


###
### Pb smelter
###
profiledata <- cbind(lam9[,1],lam10[,1],lam11[,1],lam12[,1],lam13[,1])
facdata <- cbind(fac9[,1],fac10[,1],fac11[,1],fac12[,1],fac13[,1])
faccors <- cor(facdata)
facrmse <- matrix(NA,ncol(profiledata),ncol(profiledata))
for (i in 2:ncol(profiledata)) for (j in 1:(i-1)) 
{
  facrmse[i,j] <- c(getrmse(as.matrix(facdata[,i]),as.matrix(facdata[,j])))
}
lamcors <- cor(profiledata)
lamcors




#testing


lamplots2old(lamORIGMORG,"gkey",names=mynamesmorg) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3old(facORIGMORG,"gkeyf",newdatesmorg)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdatesmorg,facname=facORIGMORG,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


