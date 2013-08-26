directory<-"C:/My Documents/research/pmf_r_interface/Weightedrose/Dr. Christensen/"
throwout<-.4
range<-pi/6
nsector<-24
radiusPb<-4
radiusZn<-4



roseplots<-function(style,wind,species,nsector=16,upper=.25,main='sample',
		measurename="cv",radius=3){
	step<-360/nsector
	deg2rad <- 180/pi
	wind <- (wind+step/2)%%360 #Values like 359 go to sector 0
	histwind <- hist(wind,breaks=seq(0,360,by=step),plot=F)
	counts<-histwind$counts
	breaks<-histwind$breaks
	mids <- (histwind$mids-step/2)/deg2rad
	contrib<-c(1:length(counts))
	for (i in 1:length(counts)){
		contrib[i]<-sum(species[breaks[i]<wind
		 & wind<breaks[i+1]])
	}
	secavcont<-contrib/histwind$counts
	totavcont<-sum(species)/length(species)
	thres<-quantile(species,1-upper)
	highwind<- wind[species > thres]
	histhigh<-hist(highwind,breaks=seq(0,360,by=step),plot=F)
	highcounts<-histhigh$counts
	step<-step/deg2rad
	if (style=='wtrose'){
	height<-secavcont/totavcont
	} else if (style=='hopke'){
	height<-highcounts/counts
	} else if (style=='simple'){
	height<-contrib # /sum(species)
	}
	if (main=='Fe') {
		srcdir<-c(10)
		srcoutput<-c(1)
	} else if (main=='Cu') {
		srcdir<-c(208,217,312,292)
		srcoutput<-c(9453,8797,3708,1900)
	} else if (main=='Zn') {
		srcdir<-c(219,010,044,045)
		srcoutput<-c(79953,35000,26490,18139)
	} else if (main=='Pb') {
		srcdir<-c(206,217,219,216)
		srcoutput<-c(117626,36818,31032,29829)
	}
	srcrad<-srcdir/deg2rad
	src<-cbind(srcrad,srcoutput)
	weight<-1:dim(src)[1]
	for (i in 1:dim(src)[1]){
		weight[i]<-sum(cos(mids-src[i,1])*(height/sum(height))
			*(src[i,2]/sum(src[,2])))
	}
	wtsum<-sum(weight)	
	wtavg<-sum(weight)/length(weight)
	petal<-1:length(mids)
	clsrc<-1:length(srcdir)
	for (i in 1:length(mids)){
		srcnum<-clsrc[cos(mids[i]-srcrad)==max(cos(mids[i]-srcrad))]
		petal[i]<-max(cos(mids[i]-srcrad))*(height[i]/max(height))*
			(srcoutput[srcnum]/
			sum(srcoutput))
	}
	closestsrc<-sum(petal)
	
	if (length(srcrad)==1){
		minwt<-1
	} else {
		minwt<-1
		for (i in 1:(length(srcrad)-1)){
			minwt<-minwt*(.5*cos(srcrad[i]-srcrad[i+1])+.5)
		}
		minwt<-minwt*(.5*cos(srcrad[1]-srcrad[length(srcrad)])+.5)
	}
	newweight<-1:dim(src)[1]
	for (i in 1:dim(src)[1]){
		wt<-(1-minwt)*(.5*cos(mids-srcrad[i]))+minwt
		newweight[i]<-sum(cos(mids-src[i,1])*(height/sum(height))
			*(src[i,2]/sum(src[,2]))*wt)
	}
	newwtsum<-sum(newweight)
	newwtavg<-sum(newweight)/length(newweight)
			
	cv<-sqrt(var(height))/mean(height)
	median<-median(height)
	max<-max(height)
	if (measurename=="cv"){measure<-cv}
	if (measurename=="max/med"){measure<-max/median}
	if (measurename=="wtsum"){measure<-wtsum}
	if (measurename=="wtavg"){measure<-wtavg}
	if (measurename=="closestsrc"){measure<-closestsrc}
	if (measurename=="newwtsum"){measure<-newwtsum}
	if (measurename=="newwtavg"){measure<-newwtavg}
	plot(c(-1.2*radius,1.2*radius),c(-1.2*radius,1.2*radius),,xlab='',
		ylab='',xaxt='n',yaxt='n',pch=' ',asp=1)
	plot(c(-1.2*radius,1.2*radius),c(-1.2*radius,1.2*radius),xlab='',
		ylab='',xaxt='n',yaxt='n',pch=' ',asp=1)
	#title(paste(style,":  ",main,"\nsectors=",nsector,"\n",
	#	measurename,"=",round(measure,3)))
	ang<-seq(0,2*pi,by=2*pi/50)
	off<-c(1.7,2.3)
	for (i in 1:(radius-1)){
		lines(sin(ang)*i,cos(ang)*i,lty="11")
		text(i,0,i,adj=off,cex=.7)	
	}
	text(radius,0,radius,adj=off,cex=.7)
	symbols(0,0,radius,add=T,inches=F,lwd=3)
	for (i in 1:length(counts)){
		w1 <- mids[i]-step/2
		arc<-(w1)+(0:50)*step/50
		lines(height[i]*c(0,sin(arc),0),
			height[i]*c(0,cos(arc),0),
			col="navy",lwd=5) #draw sector
#			col="navy",lwd=2) #draw sector
		lines(c(0,sin(w1)*radius),c(0,cos(w1)*radius),lty="11")
	}
	off<-5/deg2rad
	for (i in 1:dim(src)[1]){
		lines(c(0,sin(src[i,1])*radius*src[i,2]/max(src[,2])),
			c(0,cos(src[i,1])*radius*src[i,2]/max(src[,2])),
			col="red",lwd=5)
	}
	text(sin(mids)*1.1*radius,cos(mids)*1.1*radius,
		round(mids*deg2rad,1),font=3,cex=.85)
	if (style=='wtrose'){
	wmeas<<-measure
	} else if (style=='hopke'){
	hmeas<<-measure
	} else if (style=='simple'){
	smeas<<-measure
	}
}

elements <- read.table(file=paste(directory,"sldat2.txt",sep=""))
wind <- read.table(file=paste(directory,"src3.txt",sep=""))[,1:2]
dat1<-cbind(wind,elements)
names(dat1)<-c("date","winddir",
		"Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")
specnames<-names(dat1)[c(-1,-2)]
dat1$date<-as.Date(dat1$date,"%Y-%m-%d")

pdf(paste(directory,"Pb.raw.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",dat1$winddir,dat1$Pb,nsector,.15,"Pb","newwtsum",radius=radiusPb)
dev.off()

pdf(paste(directory,"Znraw.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",dat1$winddir,dat1$Zn,nsector,.15,"Zn","newwtsum",radius=radiusZn)
dev.off()

hr<-read.csv(file=paste(directory,"hourlywind.csv",sep=""),as.is=T)
hr<-(na.omit(hr))[,c(1,3)]
hr$date<-as.Date(hr$date,"%m/%d/%Y")
hr$dir<-hr$dir*pi/180

day <- read.table(file=paste(directory,"src3.txt",sep=""))[,1:2]
names(day)<-c("date","dir")
day$date<-as.Date(day$date,"%Y-%m-%d")
day$dir<-day$dir*pi/180

sdmat<-NULL
for (i in day$date){
	dat2<-hr[hr$date==i,]
	xcomp<-sin(dat2$dir)
	ycomp<-cos(dat2$dir)
	sdsum<-sd(xcomp)+sd(ycomp)
	sdmat<-rbind(sdmat,c(i,sdsum))
}
sdmat<-as.data.frame(sdmat)
names(sdmat)<-c("date","sdsum")


step<-2*pi/nsector
cutoff<-range*nsector/(2*pi)

hrplus<-as.data.frame(cbind(hr$date,(hr$dir+step/2)%%(2*pi))) #Values like 359 go to sector 0
names(hrplus)<-c("date","dir")
dayplus<-as.data.frame(cbind(day$date,(day$dir+step/2)%%(2*pi))) #Values like 359 go to sector 0
names(dayplus)<-c("date","dir")
histday <- hist(dayplus$dir,breaks=seq(0,2*pi,by=step),plot=F)

percmat<-NULL
for (j in 1:nsector){
	datesec<-dayplus$date[histday$breaks[j]<dayplus$dir & dayplus$dir<histday$breaks[j+1]]
	for (i in 1:length(datesec)){
		dat3<-hrplus[hrplus$date==datesec[i],]
		if (cutoff<= i & i <= nsector-cutoff){
			mids<-histday$mids
		} else if (i<cutoff) {
			dat3$dir<-dat3$dir+range
			mids<-histday$mids+range
		} else {
			dat3$dir<-dat3$dir-range
			mids<-histday$mids-range
		}
		close<-dat3[mids[j]-range<dat3[,2] & dat3[,2]<mids[j]+range,]
		percmat<-rbind(percmat,c(datesec[i],dim(close)[1]/dim(dat3)[1],(histday$mids[j]-step/2)%%(2*pi)))
	}
}
percmat<-percmat[order(percmat[,1]),]
percmat<-as.data.frame(percmat)
names(percmat)<-c("date","perc","secang")	

pdf(paste(directory,"sdsumhist.pdf",sep=""))
hist(sdmat$sdsum,main="SD Sum",col="lemonchiffon",border="red3")	
dev.off()

pdf(paste(directory,"perchist.pdf",sep=""))
hist(percmat$perc,main="Percentage in Range of Sector Mid",col="lemonchiffon",border="red3")
dev.off()

sddates<-sdmat$date[sdmat$sdsum<quantile(sdmat$sdsum,1-throwout)]
datsd<-dat1[dat1$date %in% sddates,]

pdf(paste(directory,"Pb.sdsum.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",datsd$winddir,datsd$Pb,nsector,.15,"Pb","newwtsum",radius=radiusPb)
dev.off()

pdf(paste(directory,"Zn.sdsum.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",datsd$winddir,datsd$Zn,nsector,.15,"Zn","newwtsum",radius=radiusZn)
dev.off()


percdates<-percmat$date[percmat$perc>quantile(percmat$perc,throwout)]
datperc<-dat1[dat1$date %in% percdates,]

pdf(paste(directory,"Pb.perc.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",datperc$winddir,datperc$Pb,nsector,.15,"Pb","newwtsum",radius=radiusPb)
dev.off()

pdf(paste(directory,"Znperc40.pdf",sep=""))
par(mar=c(1,1,1,1),mfrow=c(1,1))
roseplots("wtrose",datperc$winddir,datperc$Zn,nsector,.15,"Zn","newwtsum",radius=radiusZn)
dev.off()


