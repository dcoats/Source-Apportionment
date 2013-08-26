weightedrose<-function(wind,species,step=30,main='sample'){
	deg2rad <- 180/pi
	wind <- (wind+step/2)%%360 #Values like 359 go to sector 0
	histwind <- hist(wind,breaks=seq(0,360,by=step),plot=F)#use hist for breaks, counting
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
	maxheight<-max(height)
	par(mar=c(2,2,2,2))
	plot(c(-3.5,3.5),c(-3.5,3.5),,xlab='',ylab='',
		main=main,xaxt='n',yaxt='n',pch=' ',asp=1)
	ang<-seq(0,2*pi,by=2*pi/50)
	for (i in 1:3){
		lines(sin(ang)*i,cos(ang)*i,lty=2)
	}
	for (i in 1:length(counts)){
		w1 <- mids[i]-step/2
		arc<-(w1)+(0:50)*step/50
		lines(height[i]*c(0,sin(arc),0),
			height[i]*c(0,cos(arc),0),
			col="navy",lwd=3) #draw sector
		lines(c(0,sin(w1)*3),c(0,cos(w1)*3),lty=2)
	}
	text(sin(mids)*3.3,cos(mids)*3.3,
		(mids)*deg2rad,font=3)
	text(c(.85,1.85,2.85),c(-.15,-.15,-.15),c('1.0','2.0','3.0'),cex=.85)
	direction<-mids*deg2rad
	cbind(contrib,counts,secavcont,totavcont,height)
}

dat <- read.table(file="src3.txt")
names(dat)<-c("date","winddir","lead","fireworks","copper","sumsec","zinc","soil","steel","mobile","winsec")
newfacnames <- scan("fac9ChemometricsNames.txt",what="character")
newfac <- read.table("fac9Chemometrics.txt",col.names=newfacnames)

specnames<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")
sldat2 <- read.table("../sldat2.txt",col.names=specnames)

weightedrose(dat$winddir,newfac$LeadSmelter,360/24,'LeadSmelter')
weightedrose(dat$winddir,newfac$ZincSmelter,360/24,'ZincSmelter')
weightedrose(dat$winddir,newfac$CopperSmelter,360/24,'CopperSmelter')
weightedrose(dat$winddir,newfac$SteelMill,360/24,'Steel Mill')
weightedrose(dat$winddir,newfac$Soil,360/24,'Soil')
weightedrose(dat$winddir,newfac$Mobile,360/24,'Mobile')
weightedrose(dat$winddir,newfac$Fireworks,360/24,'Fireworks')
weightedrose(dat$winddir,sldat2$Pb,360/24,'Pb')
weightedrose(dat$winddir,sldat2$Fe,360/24,'Fe')
weightedrose(dat$winddir,sldat2$Cu,360/24,'Cu')
weightedrose(dat$winddir,sldat2$Zn,360/24,'Zn')
weightedrose(dat$winddir,sldat2$Si,360/24,'Si')

cuts <- 12
par(mfcol=c(2,4))
weightedrose(dat$winddir,newfac$LeadSmelter,360/cuts,'LeadSmelter')
weightedrose(dat$winddir,sldat2$Pb,360/cuts,'Pb')

weightedrose(dat$winddir,newfac$ZincSmelter,360/cuts,'ZincSmelter')
weightedrose(dat$winddir,sldat2$Zn,360/cuts,'Zn')

weightedrose(dat$winddir,newfac$CopperSmelter,360/cuts,'CopperSmelter')
weightedrose(dat$winddir,sldat2$Cu,360/cuts,'Cu')

weightedrose(dat$winddir,newfac$SteelMill,360/cuts,'Steel Mill')
weightedrose(dat$winddir,sldat2$Fe,360/cuts,'Fe')


wind<-dat$winddir(
species<-dat$lead
step<-360/16
main<-'lead'

plot(c(-3.5,3.5),c(-3.5,3.5),,xlab='',ylab='',
	main=main,xaxt='n',yaxt='n',pch=' ',asp=1)
text(sin(mids[i])*3.3,cos(mids[i])*3.3,mids[i]*deg2rad,font=3)
symbols(0,0,2,add=T,lwd=5,lty="dashed")

