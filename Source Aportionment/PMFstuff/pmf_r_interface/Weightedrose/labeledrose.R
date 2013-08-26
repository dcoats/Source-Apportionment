rose <- function(data,step=30,main='wind rose'){
  deg2rad <- 180/pi
  data <- (data+step/2)%%360# Values like 359 go to Sector 0
  histdata <- hist(data,breaks=seq(0,360,by=step),plot=F) #use hist for counting
  counts <- histdata$counts
  maxcount <- max(counts)
  mids <- (histdata$mids-step/2)/deg2rad
  step <- step/deg2rad
  plot(c(-1,1),c(-1,1),,xlab='',ylab='',
       main=main,xaxt='n',yaxt='n',pch=' ')
  for (i in 1:length(counts)){
    w1 <- mids[i]-step/2
    w2 <- mids[i]+step/2
    lines(counts[i]/maxcount*c(0,sin(w1),sin(w2),0),
          counts[i]/maxcount*c(0,cos(w1),cos(w2),0))#draw sector
    text(sin(mids[i]),cos(mids[i]),mids[i]*deg2rad)
  }
  names(counts) <- round(mids*deg2rad,3)
  counts
}



exdat <- read.csv(file="F:/EPA/winddir.csv",head=TRUE,sep=",")
x<-exdat$dec02feb03
#x <-90-x
#x<-ifelse(x>0,x,x+360)

#Test with 500 values between 0 and 360 (uniform distribution)
#rose(runif(500)*360,360/16) 
rose(x,360/16) 