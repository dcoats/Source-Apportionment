x <- rnorm(5,100,10)
y <- rnorm(5,100,10)
z <- runif(5,2000,30000)
ptnames <- c("Kip","Napoleon","Rico","Gramma","Tina")
lx <- max(x)-min(x)
ly <- max(y) - min(y)
spacemult <- .08
spaceshrink <- .5
plot(x,y,type="p",cex=sqrt(z/max(z))*10,
     xlim=c(min(x)-spacemult*lx,max(x)+spacemult*lx),
     ylim=c(min(y)-spacemult*ly,max(y)+spacemult*ly) ,xlab="Hooey" )
text(x,y+ly*.125*(z/max(z))^spaceshrink,labels=ptnames,cex=.5+z/max(z)*2)
title("Hooey",cex=3)
abline(v=100)
