ewtdrose<-function(wind,species,nsector=16,main='sample',
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
	if (main=='Fe') {
		srcdir<-c(10)
		srcoutput<-c(1)
	}
	if (main=='Cu') {
		srcdir<-c(208,217,312,292)
		srcoutput<-c(9453,8797,3708,1900)
	}
	if (main=='Zn') {
		srcdir<-c(219,010,044,045)
		srcoutput<-c(79953,35000,26490,18139)
	}
	if (main=='Pb') {
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
	title(paste("Wtdrose:  ",main,"\nsectors=",nsector,"\n",
		measurename,"=",round(measure,3)))
	ang<-seq(0,2*pi,by=2*pi/50)
	if (radius <= 1.5){
		label<-round(radius/2,1)
	} else {label<-round(radius/2)}
	lines(sin(ang)*label,cos(ang)*label,lty="11")
	text(label,0,label,adj=c(-.1,1.4),cex=.7)
	symbols(0,0,radius,add=T,inches=F,lwd=3)
	for (i in 1:length(counts)){
		w1 <- mids[i]-step/2
		arc<-(w1)+(0:50)*step/50
		lines(height[i]*c(0,sin(arc),0),
			height[i]*c(0,cos(arc),0),
			col="navy",lwd=1.5) #draw sector
		lines(c(0,sin(w1)*radius),c(0,cos(w1)*radius),lty="11")
	}
	off<-5/deg2rad
	for (i in 1:dim(src)[1]){
		lines(c(0,sin(src[i,1])*radius*src[i,2]/max(src[,2])),
			c(0,cos(src[i,1])*radius*src[i,2]/max(src[,2])),
			col="red",lwd=2)
	}
	text(sin(mids)*1.2*radius,cos(mids)*1.2*radius,
		round(mids*deg2rad,1),font=3,cex=.85)
	measure1<<-measure
}

ehopke<-function(wind,species,nsector=16,upper=.25,main='sample',
		measurename="cv",radius=3){
	step<-360/nsector
	deg2rad <- 180/pi
	wind <- (wind+step/2)%%360 #Values like 359 go to sector 0
	histwind <- hist(wind,breaks=seq(0,360,by=step),plot=F)
	counts<-histwind$counts
	breaks<-histwind$breaks
	mids <- (histwind$mids-step/2)/deg2rad
	thres<-quantile(species,1-upper)
	highwind<- wind[species > thres]
	histhigh<-hist(highwind,breaks=seq(0,360,by=step),plot=F)
	highcounts<-histhigh$counts
	height<-highcounts/counts
	step<-step/deg2rad
	if (main=='Fe') {
		srcdir<-c(10)
		srcoutput<-c(1)
	}
	if (main=='Cu') {
		srcdir<-c(208,217,312,292)
		srcoutput<-c(9453,8797,3708,1900)
	}
	if (main=='Zn') {
		srcdir<-c(219,010,044,045)
		srcoutput<-c(79953,35000,26490,18139)
	}
	if (main=='Pb') {
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
	title(paste("Hopke:  ",main,"\nsectors=",nsector," ",
		"%=",upper,"\n",measurename,"=",round(measure,3)))
	ang<-seq(0,2*pi,by=2*pi/50)
	if (radius <= 1.5){
		label<-round(radius/2,1)
	} else {label<-round(radius/2)}
	lines(sin(ang)*label,cos(ang)*label,lty="11")
	text(label,0,label,adj=c(-.1,1.4),cex=.7)
	symbols(0,0,radius,add=T,inches=F,lwd=3)
	for (i in 1:length(counts)){
		w1 <- mids[i]-step/2
		arc<-(w1)+(0:50)*step/50
		lines(height[i]*c(0,sin(arc),0),
			height[i]*c(0,cos(arc),0),
			col="navy",lwd=1.5) #draw sector
		lines(c(0,sin(w1)*radius),c(0,cos(w1)*radius),lty="11")
	}
	off<-5/deg2rad
	for (i in 1:dim(src)[1]){
		lines(c(0,sin(src[i,1])*radius*src[i,2]/max(src[,2])),
			c(0,cos(src[i,1])*radius*src[i,2]/max(src[,2])),
			col="red",lwd=2)
	}
	text(sin(mids)*1.2*radius,cos(mids)*1.2*radius,
		round(mids*deg2rad,1),font=3,cex=.85)
	measure2<<-measure
}

elements <- read.table(file="D:/EPA/Weightedrose/elements/sldat2.txt")
wind <- read.table(file="D:/EPA/Weightedrose/src3.txt")[,1:2]
dat<-cbind(wind,elements)
names(dat)<-c("date","winddir",
		"Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")
specnames<-names(dat)[c(-1,-2)]
attach(dat)
elements<-c(14,17,18,39)
meas<-rep(NA,144)
perc<-rep(NA,144)
sec<-rep(NA,144)
elem<-rep(NA,144)
style<-rep(NA,144)
compare<-cbind(meas,perc,sec,elem,style)
b<-0
for (m in c("wtsum","closestsrc","newwtsum")){
	for (p in c(.1,.15,.25)){
		for (s in c(12,16,20,24)){
			pdf(paste("D:/EPA/Weightedrose/elements/",m,"/",m,"-",p,"-",s,".pdf",sep=""),paper="letter")
			par(mar=c(.5,1,3.5,1),mfrow=c(3,3))
			for (i in 1:length(elements)){
				ewtdrose(winddir,dat[,elements[i]+2],s,specnames[elements[i]],measurename=m,radius=3)
				ehopke(winddir,dat[,elements[i]+2],s,p,specnames[elements[i]],measurename=m,radius=.6)
				b<-b+1
				compare[b,1]<-m
				compare[b,2]<-p
				compare[b,3]<-s
				compare[b,4]<-specnames[elements[i]]
				if (measure1>measure2){
				compare[b,5]<-"W"
				}else{
				compare[b,5]<-"H"
				}
			}
			dev.off()
		}
	}
}
compare<-as.data.frame(compare)
compare
write.csv(compare,file="D:/EPA/Weightedrose/elements/plotcomparison.csv")
