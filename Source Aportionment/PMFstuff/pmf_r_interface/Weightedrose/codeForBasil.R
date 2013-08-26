names<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")



sldat <- read.table("F:/EPA/sldat.txt")
slunc <- read.table("F:/EPA/slunc.txt")

sldat2 <- sldat
slunc2 <- slunc
sldat2[sldat2== -999.9] <- NA
slunc2[sldat2== -999.9] <- NA
sldat2 <- na.omit(sldat2)
slunc2 <- na.omit(slunc2)

alldat<-read.table("F:/EPA/data_for_dates.txt",header=T)
datesnew<-alldat[,1]
datesnew<-as.Date(datesnew)

newdates <- datesnew[apply(sldat==-999.9,1,sum)==0]
src<-read.table("F:/EPA/Weightedrose/fac9gkeynostart.txt")
src2<-cbind(newdates,src)
write.table(src2,file="F:EPA/Weightedrose/src2.txt",row.names=F,col.names=T)
