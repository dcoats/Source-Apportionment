wind<-read.csv("D:/EPA/Weightedrose/weirdwinddays/hourlywind.csv")
days<-c("08/28/2001","08/28/2001","08/28/2001")
dat<-wind[wind[,1]==days[1]|wind[,1]==days[2]|wind[,1]==days[3],]
dat[,3]<-dat[,3]*pi/180
dat[,2]<-dat[,2]*3.6
xcomp<-4*sin(dat$dir)
ycomp<-4*cos(dat$dir)

inlength<-dim(dat)[1]+1
x<-c(43.5*sin(206*pi/180))
y<-c(43.25*cos(206*pi/180))
plot(c(-100,100),c(-300,200),,xlab='x',ylab='y',pch=' ',asp=1,main=days[1],
	sub=days[2])
nhours<-24
col<-rainbow(dim(dat)[1]+1,start=0,end=4/6)

for (j in 1:nhours){
	for (i in 1:dim(dat)[1]){
	x[i+1]<-x[i]+dat$sp[i]*sin(dat$dir[i]+pi)
	y[i+1]<-y[i]+dat$sp[i]*cos(dat$dir[i]+pi)
	}
	size<-seq(1,5,length=inlength)
	points(x,y,pch=21,col=col,cex=size)
	cbind(x[i+1],y[i+1],col[i+1],size[i+1])
	#lines(x,y,col=col[1],lwd=2)
	cbind
	col<-col[-1]
	dat<-dat[-1,]
	x<-x[-(i+1)]
	y<-y[-(i+1)]
	Sys.sleep(.15)
}
points(43.5*sin(206*pi/180),43.25*cos(206*pi/180),pch=23,col="orange",
	bg="orange",cex=2)
points(0,0,pch=23,col="red",bg="red",cex=2)
sd(xcomp)+sd(ycomp)
