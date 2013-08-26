dat <- read.table(file="F:/EPA/Weightedrose/src3.txt")
names(dat)<-c("date","winddir","lead","fireworks","copper","sumsec",
	"zinc","soil","steel","winsec","mobile")
specnames<-c("lead","fireworks","copper","sumsec","zinc","soil","steel",
			"winsec","mobile")
attach(dat)
par(mar=c(2,2,2,2),mfrow=c(3,3))
i<-3
wind<-winddir
species<-dat[i+2]
nsector<-24
upper<-.25
main<-specnames[i]
radius<-.8