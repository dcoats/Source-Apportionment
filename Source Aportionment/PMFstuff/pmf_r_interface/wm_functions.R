facnames <- c("Fireworks","Summer Secondary","Smelter","Vehicle 1",
              "Winter Secondary","Vehicle 2")

weekendplots2<-function(i,name="gkeyw"){
####weekend plots
alldat<-read.table("dataforweekends_2.txt",header=T)
library(survival)
#dates<-seq(1,nrow(data),27)
dates<-alldat[,1]
dates<-as.character(dates)
dates<-as.date(dates)
days<-date.mdy(dates,weekday=TRUE)$weekday

weekdays<-c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")

daycounts<-array(NA,dim=c(i,7,1))

for(ii in 1:i){
for(iii in c(1:7)){
daycounts[ii,iii,1]<-mean(fac[,ii][days==iii])
}
pdf(file=mypaste(name,i,ii,".pdf") )#,width=400, height=250)
plot(daycounts[ii,,],type="l",xaxt="n",ylab="",xlab=mypaste("Mean by days of week ",i,", ",ii))
axis(1,at=c(1:7),labels=weekdays)
dev.off()}

}

facplots2<-function(fac,name="gkeyf",startdate="2001/6/1",enddate="2003/6/6"){
#start on a logical date near the beginning of the data
dates<-getdates()
for(iii in 1:ncol(fac)){
pdf(file=mypaste(name,i,iii,".pdf"))# ,width=700, height=250)
par(cex=1) 
plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=mypaste(iii," of ",i," factors"))
par(cex=.7)
axis.Date(1,at=seq(as.Date(startdate),as.Date(enddate), "months"),
format="%b-%y",las=2)
dev.off()}
}

lamplots2<-function(lam,name="gkey"){
for(iii in 1:ncol(lam)){
	pdf(file=mypaste(name,i,iii,".pdf") )#,width=500, height=250)
	par(cex=1)
	plot(lam[,iii],ylim=c(0,max(lam[,iii])+.05),ylab="",
	xlab=facnames[iii],xaxt="n")
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


lamplots3<-function(lam,name="gkey"){
for(iii in 1:ncol(lam)){
	jpeg(filename=mypaste(name,i,iii,".jpg") ,width=500, height=250)
	par(cex=1)
	plot(lam[,iii],ylim=c(0,max(lam[,iii])+.05),ylab="",
	xlab=facnames[iii],xaxt="n")
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


facplots3<-function(fac,name="gkeyf",startdate="2001/6/1",enddate="2003/6/6"){
#start on a logical date near the beginning of the data
dates<-getdates()
for(iii in 1:ncol(fac)){
jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
par(cex=1) 
plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=facnames[iii])
par(cex=.7)
axis.Date(1,at=seq(as.Date(startdate),as.Date(enddate), "months"),
format="%b-%y",las=2)
dev.off()}
}



lamplots3(lam,"gkey")  #the imagemagick batch files use these names, so if they are changed the .bat files must be changed
facplots3(fac,"gkeyf","2001/6/1","2003/6/6")
weekendplots(i,"gkeyw")

