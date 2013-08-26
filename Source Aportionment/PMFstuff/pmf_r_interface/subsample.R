#########################
#
#    Subsampling PMF Interface
#
#            Jeff Lingwall June 05
#
#########################

######
#
#   Pre-reqs
#
######
setwd("c:/work/slouis/subsampling")
rm(list=ls())
source("c:/work/slouis/subsampling/jeff's functions.R")
load(".Rdata") #changes directory and cleans up
save.image()
#data, names,prec in .Rdata
library(survival)

#CV<-27
#sulfunc<-alldat$Sulfat/CV
#CV<-21
#nitunc<-alldat$Nitrate/CV

prec[prec==0]<-.0001

tlam<-read.table("goldlam.txt");tfac<-read.table("goldfac.txt")  #the 6 source run that was good

dates<-read.table("dataforweekends_2.txt",header=T)[,1]
dates<-as.date(as.character(dates))
days<-date.mdy(dates,weekday=TRUE)$weekday


###create subsampling data sets

#Every sixth day
#results<-list()
#results$every_6<-array(dim=c(6,2))
for(i in 4:6){
	# 1 is sunday, etc . . .
	sampleday<-seq(i,dim(data)[1],6)
	sampledata<-data[sampleday,]
	sampleprec<-prec[sampleday,] #uncertainties
	tfac.i<-tfac[sampleday,]	
		
		runpmf("default_2.ini",sampledata,sampleprec,
		mypaste("one_in_six",i,"fac.txt"),mypaste("one_in_six",i,"lam.txt"),
		"junk.txt",sources=6,rows=nrow(sampledata),cols=44)

fac<-read.table(mypaste("one_in_six",i,"fac.txt"))
lam<-read.table(mypaste("one_in_six",i,"lam.txt"))
fac<-sort(lam,fac,tfac.i)$fac
lam<-sort(lam,fac,tfac.i)$lam
results$every_6[i,1]<-AAE(fac,tfac.i)
results$every_6[i,2]<-AAE(lam,tlam)

#at this point, I need to make the graphs
lamplots.2(lam,"every_6")
facplots.2(fac,"every_6f")
weekendplots.2(6,"every_6w")
batfile<- readLines(con="every_6.bat")
batfile[7]<-mypaste("convert -append every_6a61.jpg every_6a62.jpg every_6a65.jpg every_6a63.jpg every_6a64.jpg every_6a66.jpg every_6_",i,"_n6.jpg")
writeLines(batfile,con="every_6.bat")
system("every_6.bat",show.output.on.console = TRUE,invisible=TRUE)
}


#1 in 6 random
results$random_6<-array(dim=c(12,2))
for(i in 1:6){
if(i!=4){
	interval<-seq(1,dim(data)[1],6) #sampling interval
	sampleday<-numeric(length(interval))
	for(ii in 1:length(interval)) sampleday[ii]<-interval[ii]+ sample(0:5,1)
	while(sampleday[length(sampleday)] > dim(data)[1] ) 
	sampleday[ii]<-interval[ii]+ sample(0:5,1)
	sampledata<-data[sampleday,]
	sampleprec<-prec[sampleday,] #uncertainties
	tfac.i<-tfac[sampleday,]	

	#	runpmf("default_2.ini",sampledata,sampleprec,
	#	mypaste("random_six",i,"fac.txt"),mypaste("random_six",i,"lam.txt"),
	#	"junk.txt",sources=6,rows=nrow(sampledata),cols=44)
fac<-read.table(mypaste("random_six",i,"fac.txt"))
lam<-read.table(mypaste("random_six",i,"lam.txt"))
fac<-sort(lam,fac,tfac.i)$fac
lam<-sort(lam,fac,tfac.i)$lam
results$random_6[i,1]<-AAE(fac,tfac.i)
results$random_6[i,2]<-AAE(lam,tlam)
#lamplots.2(lam,"every_6")
#facplots.2(fac,"every_6f")
#weekendplots.2(6,"every_6w")
#batfile<- readLines(con="every_6.bat")
#batfile[7]<-mypaste("convert -append every_6a61.jpg every_6a62.jpg every_6a65.jpg every_6a63.jpg every_6a64.jpg every_6a66.jpg random_6_",i,"_n6.jpg")
#writeLines(batfile,con="every_6.bat")
#system("every_6.bat",show.output.on.console = TRUE,invisible=TRUE)
	}}

results

#they are basically the same as far as AAE goes, with the random draw
#doing perhaps a little worse on the estimation of F

#Now, let's look at the plots and see how they look
















#################################################
#
#
#    FUNCTIONS USED ABOVE
#
#
#################################################

runpmf<-function(inifilename,data,prec,F,Lambda,Q,sources=5,rows=788,cols=23){
sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
inifile<- readLines(con=inifilename)
inifile[38]<-paste("    30   T \"OLD    \"  2000  \" pmf_data.txt    \"",sep="")
inifile[39]<-paste("    31   T \"OLD    \"  2000  \" pmf_unc.txt    \"",sep="")
inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
writeLines(inifile,con=inifilename)

mywrite.table(data,file="pmf_data.txt",quote=FALSE,na="-999.9")
mywrite.table(prec,file="pmf_unc.txt",quote=FALSE,na="1")

system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

weekendplots.2<-function(i,name="gkeyw"){
####weekend plots
#alldat<-read.table("dataforweekends_2.txt",header=T)
#library(survival)
#dates<-seq(1,nrow(data),27)
#dates<-alldat[,1]
#dates<-as.character(dates)
#dates<-as.date(dates)
#days<-date.mdy(dates,weekday=TRUE)$weekday
days.i<-days[sampleday]

weekdays<-c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")

daycounts<-array(NA,dim=c(i,7,1))

for(ii in 1:i){
for(iii in c(1:7)){
daycounts[ii,iii,1]<-mean(fac[,ii][days.i==iii])}
jpeg(filename=mypaste(name,i,ii,".jpg") ,width=400, height=250)
plot(daycounts[ii,,],type="l",xaxt="n",ylab="",xlab=mypaste("Mean by days of week ",i,", ",ii))
axis(1,at=c(1:7),labels=weekdays)
dev.off()}

}

facplots.2<-function(fac,name="gkeyf",startdate="2001/6/1",enddate="2003/6/6"){
#start on a logical date near the beginning of the data
#dates<-getdates()
i<-ncol(fac)
dates.i<-as.Date(dates[sampleday])
for(iii in 1:ncol(fac)){
jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
par(cex=1) 
plot(dates.i,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=mypaste(iii," of ",i," factors"))
par(cex=.7)
axis.Date(1,at=seq(as.Date(startdate),as.Date(enddate), "months"),
format="%b-%y",las=2)
dev.off()}
}

lamplots.2<-function(lam,name="gkey"){
for(iii in 1:ncol(lam)){
i<-ncol(lam)
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



