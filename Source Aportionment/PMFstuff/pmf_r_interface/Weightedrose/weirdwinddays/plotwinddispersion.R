windplot<-function(wind=file,days,style="cloud",start=1/12,end=4/6,
	xpar=c(-60,60),ypar=c(-150,100),vw=NA,conc=NA){

dat<-wind[wind[,1]>=as.Date(days[1]) & wind[,1]<=as.Date(days[2]),]

dat[,3]<-dat[,3]*pi/180
dat[,2]<-dat[,2]*3.6
xcomp<-4*sin(dat$dir)
ycomp<-4*cos(dat$dir)
varsum<-sd(xcomp)+sd(ycomp)
inlength<-dim(dat)[1]+1
x<-c(43.5*sin(206*pi/180))
y<-c(43.25*cos(206*pi/180))
s<-c(NA)
plot(xpar,ypar,,xlab='x',ylab='y',pch=' ',asp=1)
title(paste(days[1],",  dir=",vw,
	"\n svar=",round(varsum,3),", c=",round(conc,3),sep=""))
xcoord<-cbind(rep(NA,dim(dat)[1]+1))
ycoord<-cbind(rep(NA,dim(dat)[1]+1))
sizematrix<-cbind(rep(NA,dim(dat)[1]+1))
size<-seq(1,(dim(dat)[1]+1)/6,length=dim(dat)[1]+1)
col<-rainbow(dim(dat)[1]+1,start=start,end=end)
for (j in 1:dim(dat)[1]){
	x[j]<-c(43.25*sin(206*pi/180))
	y[j]<-c(43.25*cos(206*pi/180))
	s[j]<-size[1]
	for (i in j:dim(dat)[1]){

		x[i+1]<-x[i]+dat$sp[i]*sin(dat$dir[i]+pi)
		y[i+1]<-y[i]+dat$sp[i]*cos(dat$dir[i]+pi)
		s[i+1]<-size[i+2-j]
	}
	xcoord<-cbind(xcoord,x)
	ycoord<-cbind(ycoord,y)
	sizematrix<-cbind(sizematrix,s)
	x<-rep(NA,length(x))
	y<-rep(NA,length(y))
	s<-rep(NA,length(s))
}
x[dim(dat)[1]+1]<-c(43.25*sin(206*pi/180))
y[dim(dat)[1]+1]<-c(43.25*cos(206*pi/180))
s[dim(dat)[1]+1]<-size[1]

xcoord<-cbind(xcoord,x)[,-1]
ycoord<-cbind(ycoord,y)[,-1]
sizematrix<-cbind(sizematrix,s)[,-1]
xcoord<-rbind(xcoord,col)
ycoord<-rbind(ycoord,col)

points(43.5*sin(206*pi/180),43.25*cos(206*pi/180),pch=23,col="orange",
	bg="orange",cex=1)
points(0,0,pch=23,col="red",bg="red",cex=1)


if (style=="cloud"){
	for (i in 1:dim(dat)[1]+1){
		points(xcoord[i,],ycoord[i,],pch=21,col=col[i],cex=sizematrix[i,])
		Sys.sleep(.1)
	}
} else if (style=="path"){
	for (i in 1:dim(dat)[1]+1){
		points(xcoord[i,1],ycoord[i,1],pch=20,col="blue",cex=1)
		Sys.sleep(.1)
	}
	lines(xcoord[,1],ycoord[,1],col="red",lwd=1)
}
points(43.5*sin(206*pi/180),43.25*cos(206*pi/180),pch=23,col="orange",
	bg="orange",cex=1)
points(0,0,pch=23,col="red",bg="red",cex=1)
sd(xcomp)+sd(ycomp)
}

elements <- read.table(file="F:/EPA/Weightedrose/elements/sldat2.txt")
wind <- read.table(file="F:/EPA/Weightedrose/src3.txt",as.is=T)[,1:2]
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
dat1<-dat[322.5<winddir & winddir<352.5 & Pb>thres,c(1,2,39+2)]
dat2<-dat[67.5<winddir & winddir<97.5 & Pb>thres,c(1,2,39+2)]
dat3<-dat[202.5<winddir & winddir<217.5 & Pb>thres,c(1,2,39+2)]
dat4<-dat[37.5<winddir & winddir<52.5 & Pb>thres,c(1,2,39+2)]

set<-dat4

file<-read.csv("F:/EPA/Weightedrose/weirdwinddays/hourlywind.csv",as.is=T)
file[,1]<-as.Date(file[,1],"%m/%d/%Y")

pdf("F:/EPA/Weightedrose/weirdwinddays/dat4path.pdf",paper="letter")
par(mar=c(2,2,3,2),mfrow=c(3,3))
for (i in 1:dim(set)[1]){
	windplot(file,c(set[i,1],set[i,1]),"path",1/12,4/6,
		c(-60,60),c(-150,50),set[i,2],set[i,3])
}
dev.off()

pdf("F:/EPA/Weightedrose/weirdwinddays/highdays7.pdf",paper="letter")
par(mar=c(2,2,3,2),mfrow=c(2,2))
windplot(file,c("2002-07-01","2002-07-05"),"cloud",1/12,4/6,
	c(-100,60),c(-150,50),52,.105)
dev.off()