setwd("c:/work/slouis")
source("c:/work/jeff's functions.R")
save.image()
load(".Rdata")

alldat<-read.table("sldat_naremoved.txt",header=T)
names(alldat)
data<-cbind(alldat$NAXC,alldat$MGXC,alldat$ALXC,alldat$SIXC,alldat$PHXC,
		alldat$SUXC,alldat$CLXC,alldat$KPXC,alldat$CAXC,alldat$TIXC,
		alldat$VAXC,alldat$CRXC,alldat$MNXC,alldat$FEXC,alldat$COXC,
		alldat$NIXC,alldat$CUXC,alldat$ZNXC,alldat$GAXC,alldat$ASXC,
		alldat$SEXC,alldat$BRXC,alldat$RBXC,alldat$SRXC,alldat$YTXC,
		alldat$ZRXC,alldat$MOXC,alldat$PDXC,alldat$AGXC,alldat$CDXC,
		alldat$INXC,alldat$SNXC,alldat$SBXC,alldat$BAXC,
		alldat$LAXC,alldat$AUXC,alldat$HGXC,alldat$TLXC,alldat$PBXC,
		alldat$URXC,alldat$OC,alldat$EC,alldat$Sulfate,alldat$Nitrate)



for(i in 3:10){
assign(mypaste("lams",i),as.matrix(read.table(mypaste("lam",i,".txt"))))
assign(mypaste("facs",i),as.matrix(read.table(mypaste("fac",i,".txt"))))}


runpmf<-function(inifilename,data,Ftest,Ltest,Qtest,sources=5,rows=788,cols=23){
#inifilename<-"default.ini" #CHANGE 
sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
inifile<- readLines(con=inifilename)
#data<-"ydat1.txt"
#Ftest<-"FTest.txt"
#Ltest<-"Ltest.txt"
#Qtest<-"Qtest.txt"
#sources<-5
#rows<-788
#cols<-23
inifile[38]<-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",Ftest,"    \"",sep="")
inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Ltest,"    \"",sep="")
inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Qtest,"    \"",sep="")
inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
writeLines(inifile,con=inifilename)

system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

weekendplots<-function(i){
####weekend plots
alldat<-read.table("dataforweekends_2.txt",header=T)
library(survival)
#dates<-seq(1,nrow(data),27)
dates<-alldat[,1]
dates<-as.character(dates)
dates<-as.date(dates)
days<-date.mdy(dates,weekday=TRUE)$day

weekdays<-c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")

daycounts<-array(NA,dim=c(i,7,1))

for(ii in 1:i){
for(iii in c(1:7)){
daycounts[ii,iii,1]<-mean(fac[,ii][days==iii])
}
jpeg(filename=mypaste("gkeyw",i,ii,".jpg") ,width=400, height=250)
plot(daycounts[ii,,],type="l",xaxt="n",ylab="",xlab=mypaste("Mean by days of week ",i,", ",ii))
axis(1,at=c(1:7),labels=weekdays)
dev.off()}

}

facplots<-function(fac,name="gkeyf"){
dates<-getdates()
for(iii in 1:ncol(fac)){
jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
par(cex=1) 
plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=mypaste(iii," of ",i," factors"))
par(cex=.7)
axis.Date(1,at=seq(as.Date("2001/6/1"),as.Date("2003/6/6"), "months"),
format="%b-%y",las=2)
dev.off()}
}

gkeyprep<-function(uncertainty,i){
meserror<-read.table("slunc.txt")
uncertainty<-rbind(uncertainty,meserror)
write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)

#adds a matrix of zeros on top of the starting array F, obtained from a previous run

#zeros<-array(0,dim=c(i,i))
#gnot<-read.table(paste("testfac",i,".txt",sep=""))
#keygnot<-rbind(zeros,gnot)
#write.table(keygnot,file="keygnot.txt",row.names=FALSE,col.names=FALSE)

key<-t(lambda)
#startlam<-read.table(mypaste("startlam",i,".txt"))
write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)

data<-read.table(paste("sldat",".txt",sep=""))
keydata<-rbind(as.matrix(key),as.matrix(data))
write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}


names<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")

lamplots<-function(lam){
for(iii in 1:ncol(lam)){
	jpeg(filename=mypaste("gkey",i,iii,".jpg") ,width=500, height=250)
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
alldat<-read.table("data_for_dates.txt",header=T)
dates<-alldat[,1]
dates<-as.Date(dates)
return(dates)}

################End of RUNPMF function
########
#
#   Gkeying
#
########
#gets the default values to use as starting values
i<-6

smelter<-apply(lams7[,c(1,3)],1,mean)   #smelter

lambda<-cbind(lams7[,c(2,4,6)],smelter,lams7[,c(5,7)])
i<-6

#####################################   best results
uncertainty<-array(.0001,dim=c(6,44))
uncertainty[1,]<-.0001         #fireworks  
uncertainty[2,]<-.01       #summer secondary 
uncertainty[3,]<-.01       #winter secondary 
uncertainty[4,]<-.1    #smelter (artificial) 
uncertainty[5,]<-100    #vehicle 
uncertainty[6,]<-100    #vehicle    


#testing
uncertainty<-array(.0001,dim=c(6,44))
uncertainty[1,]<-.0001         #fireworks  
uncertainty[2,]<-.01       #summer secondary 
uncertainty[3,]<-.01       #winter secondary 
uncertainty[4,]<-.1    #smelter (artificial) 
uncertainty[5,]<-100    #vehicle 
uncertainty[6,]<-100  #vehicle    


gkeyprep(uncertainty,6)
runpmf("gkey_combine_test6.ini","keydat6.txt",mypaste("testfac",i,".txt"),mypaste("testlam",i,".txt"),
mypaste("testq",i,".txt"),i,dim(data)[1]+i,dim(data)[2])

lam<-read.table("testlam6.txt")
fac<-read.table("testfac6.txt")

fac<-fac[-c(1:i),]    #if gkeying is used
lamplots(lam)
facplots(fac)
weekendplots(i)
system("multipanel_gkey6.bat")

###########################################extra goodies for plots
sources<-c("Fireworks","Summer Secondary","Smelter","Vehicular 1",
	     "Winter Secondary","Vehicular 2")

dates<-getdates()
name<-"gkeyf"
for(iii in 1:ncol(fac)){
jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
par(cex=1) 
plot(dates,fac[,iii],type="l",col="blue",xlab=sources[iii],xaxt="n",ylab=mypaste(iii," of ",i," factors"))
par(cex=.6)
axis.Date(1,at=seq(as.Date("2001/6/1"),as.Date("2003/6/6"), "months"),
format="%b-%y",las=2)
dev.off()}
system("multipanel_gkey6.bat")


#################################################
#
#
#
get(getOption("device"))()
#
#################################################