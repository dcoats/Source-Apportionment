Ysim <- t(read.table("Y_sim.txt",sep=","))
Ysim
Yunc <- read.table("slunc2.txt")

  write(t(Yunc),file="Yunc.txt",ncol=ncol(Yunc))
  write(t(Ysim),file="Ysim.txt",ncol=ncol(Ysim))

trulam <- read.table("lam9PMFout.txt")
trufac <- read.table("fac9PMFout.txt")

compmats <- function(x,trux)
{
  temp <- abs(x-trux)/trux 
  apply(temp,2,median)
}

i <- 9
newdates <- datesnew[apply(Ysim==-999.9,1,sum)==0]
mynames <- names(sldat2)

runpmfplain("default_nogkey.ini","Ysim.txt","Yunc.txt",mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),fpeak= 0,sources=9,nrow(Ysim),ncol(Ysim))   

limsim9 <- lam <- read.table(mypaste("testlam",i,".txt"))
fac9 <- fac <- read.table(mypaste("testfac",i,".txt"))
lamplots2old(lam,"gkey",names=mynames) 
#lamplotslog(lam,"gkey",names=mynames) 
#facplots3old(fac,"gkeyf",newdates)
facplots3old(fac,"gkeyf",1:661)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,1:661,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 

compmats(lam,trulam)
plot(c(as.matrix(lam)), c(as.matrix(trulam)))


