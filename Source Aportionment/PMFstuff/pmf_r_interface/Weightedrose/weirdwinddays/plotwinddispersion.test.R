file<-read.csv("D:/EPA/Weightedrose/weirdwinddays/hourlywind.csv",as.is=T)
file[,1]<-as.Date(file[,1],"%m/%d/%Y")
wind<-file
days<-c("2001-09-26","2001-09-30")
style<-"cloud"
start<-1/12
end<-4/6
xpar<-c(-60,60)
ypar<-c(-150,100)



