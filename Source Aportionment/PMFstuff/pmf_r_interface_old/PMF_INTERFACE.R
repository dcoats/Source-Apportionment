#########################
#
#    Simple PMF-R Interface
#
#            Jeff Lingwall June 05
#
#########################

######
#
#   Pre-reqs
#
######
setwd("c:/my documents/research/pmf_r_interface")
source("jeff's functions.R")
#load(".Rdata")

alldat<-read.table("sldat_naremoved.txt",header=T)

data<-cbind(alldat$NAXC,alldat$MGXC,alldat$ALXC,alldat$SIXC,alldat$PHXC,
		alldat$SUXC,alldat$CLXC,alldat$KPXC,alldat$CAXC,alldat$TIXC,
		alldat$VAXC,alldat$CRXC,alldat$MNXC,alldat$FEXC,alldat$COXC,
		alldat$NIXC,alldat$CUXC,alldat$ZNXC,alldat$GAXC,alldat$ASXC,
		alldat$SEXC,alldat$BRXC,alldat$RBXC,alldat$SRXC,alldat$YTXC,
		alldat$ZRXC,alldat$MOXC,alldat$PDXC,alldat$AGXC,alldat$CDXC,
		alldat$INXC,alldat$SNXC,alldat$SBXC,alldat$BAXC,
		alldat$LAXC,alldat$AUXC,alldat$HGXC,alldat$TLXC,alldat$PBXC,
		alldat$URXC,alldat$OC,alldat$EC,alldat$Sulfate,alldat$Nitrate)

write.table(data,file="sldat.txt",quote=F,na="-999.9",row.names=FALSE,col.names=FALSE)

names<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")


########
#
#   Ladies and Gentlemen . . .
#
########


##create a target profile matrix, lambda
#note-- I assume lambda is p by k

lams7<-read.table("lam7.txt")  #the one I was working with
smelter<-apply(lams7[,c(1,3)],1,mean)   #smelter

lambda<-cbind(lams7[,c(2,4,6)],smelter,lams7[,c(5,7)])

i<-min(dim(lambda))
k<-max(dim(lambda))
n<-749


#testing
uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-.0001         #fireworks  
uncertainty[2,]<-.01           #summer secondary 
uncertainty[3,]<-.01           #winter secondary 
uncertainty[4,]<-.1            #smelter (artificial) 
uncertainty[5,]<-1           #vehicle 
uncertainty[6,]<-1           #vehicle    

gkeyprep(lambda,uncertainty,"slunc.txt","sldat.txt")     #gets ready to gkey

runpmf("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),i,n+i,k)                  # ini file must be set up for gkeying, 
                                                       #   if no gkeying is used, change "keydat" to "sldat"

runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= -0.1,i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"

lam<-read.table(mypaste("testlam",i,".txt"))
fac<-read.table(mypaste("testfac",i,".txt"))

fac<-fac[-c(1:i),]    #if gkeying is used
lamplots(lam,"gkey")  #the imagemagick batch files use these names, so if they are changed the .bat files must be changed
facplots(fac,"gkeyf","2001/6/1","2003/6/6")
weekendplots(i,"gkeyw")

#the following command doesn't want to work on your system.  Maybe R 2.02 has a bug--
#run the batch file from the command line (I have batch files for 5,6,7 sources)
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) #see License in this directory for the imagemagick license
#system(mypaste("multipanel_gkey6EXAMPLE.bat"),invisible=TRUE) 
Q<-getQ(mypaste("testq",i,".txt"))
Q


alldat<-read.table("dataforweekends_2.txt",header=T)
library(survival)
#dates<-seq(1,nrow(data),27)
dates<-alldat[,1]
datess<-alldat[,1]
dates<-as.character(dates)
dates<-as.date(dates)
days<-date.mdy(dates,weekday=TRUE)$weekday




### Compare Fpeak and Q ###
fpeakvals <- c(-.1,0,.1,.2,.3,.4)
qvals <- rep(NA,length(fpeakvals))
for (jj in 1:length(fpeakvals))
{
   runpmf2("default_gkey_6source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
      mypaste("testq",i,".txt"),fpeak=fpeakvals[jj],i,n+i,k)        #   if no gkeying is used, change "keydat" to "sldat"
   qvals[jj] <- getQ(mypaste("testq",i,".txt"))
   cat(paste("VALUE",jj,"COMPLETED"))
}
## rbind(fpeakvals,qvals)
#              [,1]       [,2]       [,3]     [,4]       [,5]       [,6]     [,7]
#fpeakvals     -1.0     -0.666     -0.333      0.0      0.333      0.666      1.0
#qvals     168788.9 169001.750 168734.484 168708.0 168707.969 171602.094 172083.8

save.image()


################################################NOTES
#GKEY INI FILES
# default_gkey_5source.ini
# default_gkey_6source.ini
# default_gkey_7source.ini
#
#NO GKEY FILE WITH SAME SETTINGS OTHERWISE
# default_nogkey.ini
#
#The functions below aren't generic-- if data other than the st louis data with the files set up how I have them
#are run, they will need to be changed (especially the functions concerned with the dates)
#
#
#################################################

#################################################
# 
#    this line comes in handy because creating the jpegs kills the internal graphics device
#    so if you want to plot something this calls a new one
#
get(getOption("device"))()
#
#################################################


#################################################
#
#
#    FUNCTIONS USED ABOVE
#
#
#################################################

runpmf<-function(inifilename,data,F,Lambda,Q,sources=5,rows=788,cols=23){
sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
inifile<- readLines(con=inifilename)
inifile[38]<-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
writeLines(inifile,con=inifilename)
system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

runpmf2 <- function(inifilename,data,F,Lambda,Q,fpeak=0,sources=5,rows=788,cols=23)
{
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[9] <-paste("     ",fpeak)
  inifile[38]<-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

weekendplots<-function(i,name="gkeyw"){
####weekend plots
alldat<-read.table("dataforweekends_2.txt",header=T)
library(survival)
#dates<-seq(1,nrow(data),27)
dates<-alldat[,1]
datess<-alldat[,1]
dates<-as.character(dates)
dates<-as.date(dates)
days<-date.mdy(dates,weekday=TRUE)$weekday

weekdays<-c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")

daycounts<-array(NA,dim=c(i,7,1))

for(ii in 1:i){
for(iii in c(1:7)){
daycounts[ii,iii,1]<-mean(fac[,ii][days==iii])
}
jpeg(filename=mypaste(name,i,ii,".jpg") ,width=400, height=250)
plot(daycounts[ii,,],type="l",xaxt="n",ylab="",xlab=mypaste("Mean by days of week ",i,", ",ii))
axis(1,at=c(1:7),labels=weekdays)
dev.off()}

}

facplots<-function(fac,name="gkeyf",startdate="2001/6/1",enddate="2003/6/6"){
#start on a logical date near the beginning of the data
dates<-getdates()
for(iii in 1:ncol(fac)){
jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
par(cex=1) 
plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=mypaste(iii," of ",i," factors"))
par(cex=.7)
axis.Date(1,at=seq(as.Date(startdate),as.Date(enddate), "months"),
format="%b-%y",las=2)
dev.off()}
}

gkeyprep<-function(lambdakey,uncertainty,ambienterror="slunc.txt",datamatrix="sldat.txt"){

## with this variant of gkeying I'm using the key as a starting value for lambda, and 
## not using a starting value for F (sorting it out would be a pain)
i<-min(dim(lambdakey))
meserror<-read.table(ambienterror)
uncertainty<-rbind(uncertainty,meserror)
write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)

#adds a matrix of zeros on top of the starting array F, obtained from a previous run

#zeros<-array(0,dim=c(i,i))
#gnot<-read.table(paste("testfac",i,".txt",sep=""))
#keygnot<-rbind(zeros,gnot)
#write.table(keygnot,file="keygnot.txt",row.names=FALSE,col.names=FALSE)

key<-t(lambdakey)
write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)

data<-read.table(datamatrix)
keydata<-rbind(as.matrix(key),as.matrix(data))
write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}

lamplots<-function(lam,name="gkey"){
for(iii in 1:ncol(lam)){
	jpeg(filename=mypaste(name,i,iii,".jpg") ,width=500, height=250)
	par(cex=1)
	plot(lam[,iii],ylim=c(0,max(lam[,iii])+.05),ylab="",
	xlab=mypaste("Source ",i,", ",iii),xaxt="n")
	par(cex=.6)
	axis(1,at=1:length(lam[,iii]),label=names,las=2)
	###  makes a box around significant values
	par(cex=1)	
	for(ii in 1:nrow(lam)){
	if(lam[ii,iii]>=.01){
	rect(ii-.5,par("usr")[3],ii+.5,par("usr")[4]-.01,col="grey90",border="transparent")
	}}
	points(lam[,iii],bty="o") #redraws the points on top of the shaded boxes
	#end of box maxing
	#print(par("usr"))
	dev.off()
	
	}
}

getdates<-function(){
#i used excel to get the dates in the right format for R to read
alldat<-read.table("data_for_dates.txt",header=T)
dates<-alldat[,1]
dates<-as.Date(dates)
return(dates)}

getQ<-function(filename){
Q<-readLines(con=filename)
Q<-Q[1]
Q<-unlist(strsplit(Q,"=",fixed=TRUE))[2]
Q<-as.numeric(Q)
return(Q)}



