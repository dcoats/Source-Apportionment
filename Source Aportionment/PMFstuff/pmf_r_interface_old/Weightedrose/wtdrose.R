wtdrose<-function(wind,species,nsector=16,main='sample',
		measurename="cv",radius=3){
	step<-360/nsector
	deg2rad <- 180/pi
	wind <- (wind+step/2)%%360 #Values like 359 go to sector 0
	histwind <- hist(wind,breaks=seq(0,360,by=step),plot=F)
	counts<-histwind$counts
	breaks<-histwind$breaks
	mids <- (histwind$mids-step/2)/deg2rad
	step<-step/deg2rad
	contrib<-c(1:length(counts))
	for (i in 1:length(counts)){
		contrib[i]<-sum(species[breaks[i]<wind
		 & wind<breaks[i+1]])
	}
	secavcont<-contrib/histwind$counts
	totavcont<-sum(species)/length(species)
	height<-secavcont/totavcont
	cv<-sqrt(var(height))/mean(height)
	var<-var(height)
	median<-median(height)
	max<-max(height)
	if (measurename=="cv"){measure<-cv}
	if (measurename=="var"){measure<-var}
	if (measurename=="max/med"){measure<-max/median}
	plot(c(-1.2*radius,1.2*radius),c(-1.2*radius,1.2*radius),,xlab='',
		ylab='',xaxt='n',yaxt='n',pch=' ',asp=1)
	title(paste("Wtdrose:  ",main,"\nsectors=",nsector,"\n",
		measurename,"=",round(measure,3)))
	ang<-seq(0,2*pi,by=2*pi/50)
	if (radius<=.5){
		for (i in seq(0,radius,by=.1)){
			lines(sin(ang)*i,cos(ang)*i,lty="11")
		}
		i<-seq(0,radius,by=.1)
		text(i,0,i,adj=c(1.2,1.2),cex=.7)
	}
	if (radius > .5 & radius <= 1){
		for (i in seq(0,radius,by=.2)){
			lines(sin(ang)*i,cos(ang)*i,lty="11")
		}
		i<-seq(0,radius,by=.2)
		text(i,0,i,adj=c(1.2,1.2),cex=.7)
	}
	if (radius > 1 & radius <= 2){
		for (i in seq(0,radius,by=.5)){
			lines(sin(ang)*i,cos(ang)*i,lty="11")
		}
		i<-seq(0,radius,by=.5)
		text(i,0,i,adj=c(1.2,1.2),cex=.7)
	}
	if (radius > 2) {
		for (i in 1:radius){
			lines(sin(ang)*i,cos(ang)*i,lty="11")
		}
		i<-1:round(radius)
		text(i,0,i,adj=c(1.2,1.2),cex=.7)
	}
	symbols(0,0,radius,add=T,inches=F,lwd=3)
	for (i in 1:length(counts)){
		w1 <- mids[i]-step/2
		arc<-(w1)+(0:50)*step/50
		lines(height[i]*c(0,sin(arc),0),
			height[i]*c(0,cos(arc),0),
			col="navy",lwd=2) #draw sector
		lines(c(0,sin(w1)*radius),c(0,cos(w1)*radius),lty="11")
	}
	if (main=='lead') {
		lines(c(0,sin(206/deg2rad)*radius),c(0,cos(206/deg2rad)*radius),col="red",lwd=1)
		lines(c(0,sin(217/deg2rad)*radius),c(0,cos(217/deg2rad)*radius),col="red",lwd=1)
		text(sin(206/deg2rad)*radius*.9,cos(206/deg2rad)*radius*.9,1,adj=c(-1,1),col="red",cex=.8)
		text(sin(217/deg2rad)*radius*.9,cos(217/deg2rad)*radius*.9,2,adj=c(1,-1),col="red",cex=.8)
	}
	if (main=='copper') {
		lines(c(0,sin(208/deg2rad)*radius),c(0,cos(208/deg2rad)*radius),col="red",lwd=1)
		lines(c(0,sin(217/deg2rad)*radius),c(0,cos(217/deg2rad)*radius),col="red",lwd=1)
		text(sin(208/deg2rad)*radius*.9,cos(208/deg2rad)*radius*.9,1,adj=c(-1,1),col="red",cex=.8)
		text(sin(217/deg2rad)*radius*.9,cos(217/deg2rad)*radius*.9,2,adj=c(1,-1),col="red",cex=.8)
	}
	if (main=='zinc') {
		lines(c(0,sin(219/deg2rad)*radius),c(0,cos(219/deg2rad)*radius),col="red",lwd=1)
		lines(c(0,sin(10/deg2rad)*radius),c(0,cos(10/deg2rad)*radius),col="red",lwd=1)
		text(sin(219/deg2rad)*radius*.9,cos(219/deg2rad)*radius*.9,1,adj=c(-1,1),col="red",cex=.8)
		text(sin(10/deg2rad)*radius*.9,cos(10/deg2rad)*radius*.9,2,adj=c(1.5,.5),col="red",cex=.8)
	}
	if (main=='steel') {
		lines(c(0,sin(10/deg2rad)*radius),c(0,cos(10/deg2rad)*radius),col="red",lwd=1)
		text(sin(10/deg2rad)*radius*.9,cos(10/deg2rad)*radius*.9,1,adj=c(-1,1),col="red",cex=.8)
	}
	text(sin(mids)*1.2*radius,cos(mids)*1.2*radius,
		round(mids*deg2rad,1),font=3,cex=.85)
	direction<-mids*deg2rad
	cbind(direction,contrib,counts,secavcont,totavcont,height)
}

dat <- read.table(file="D:/EPA/Weightedrose/src3.txt")
names(dat)<-c("date","winddir","lead","fireworks","copper","sumsec",
	"zinc","soil","steel","winsec","mobile")
specnames<-c("lead","fireworks","copper","sumsec","zinc","soil","steel",
			"winsec","mobile")
attach(dat)
pdf("D:/EPA/Weightedrose/wtdrose24cvadj.pdf",paper="letter")
par(mar=c(.5,1,3.5,1),mfrow=c(3,3))
for (i in 1:length(specnames)){
wtdrose(winddir,dat[,i+2],24,specnames[i],measurename="cv",radius=4)
}
dev.off()