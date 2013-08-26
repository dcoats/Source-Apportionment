# cvi <- 3
#reperri <- 0
    reperri <- reperri + 1
    uncerrs <- uncs * matrix(rnorm(nrow(uncs)*ncol(uncs),0+1,
                      cvvals[cvi]),nrow(uncs),ncol(uncs))
    #uncerrs <- uncs
    # uncerrs <- uncerrs*10
    write(t(uncerrs),file="uncerrs.txt",ncol=ncol(uncerrs))
    runpmfplain("default_nogkey.ini","morgc.txt", "uncerrs.txt" ,
       mypaste("testfac",numfacs,".txt"),
       mypaste("testlam",numfacs,".txt"),
       mypaste("testq",numfacs,".txt"),fpeak= 0.1,numfacs,nrow(concs),ncol(concs))    

     lam<-read.table(mypaste("testlam",numfacs,".txt"))
     fac<-read.table(mypaste("testfac",numfacs,".txt"))
     lamarrayerrs[,,cvi,reperri] <- as.matrix(lam)
     facarrayerrs[,,cvi,reperri] <- as.matrix(fac)

reperri


