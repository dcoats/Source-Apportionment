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
thres<-quantile(Pb,.9)
dat[322.5<winddir & winddir<352.5 & Pb>thres,c(1,2,39+2)]
dat[67.5<winddir & winddir<97.5 & Pb>thres,c(1,2,39+2)]
dat[202.5<winddir & winddir<217.5 & Pb>thres,c(1,2,39+2)]
dat[37.5<winddir & winddir<52.5 & Pb>thres,c(1,2,39+2)]