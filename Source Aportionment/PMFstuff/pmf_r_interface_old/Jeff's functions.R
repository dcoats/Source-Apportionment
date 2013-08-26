


#################################################
#
#
#    FUNCTIONS USED ABOVE
#
#
#################################################

runpmf<-function(inifilename,data,F,Lambda,Q,sources=5,rows=788,cols=23){
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[38]<-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

runpmf2 <- function(inifilename,data,F,Lambda,Q,fpeak=0,sources=5,rows=788,cols=23)
{
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[9] <-paste("     ",fpeak)
  inifile[38]<-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}

weekendplots<-function(i,name="gkeyw"){
  ####weekend plots
  alldat<-read.table("dataforweekends_2.txt",header=T)
  library(survival)
  #dates<-seq(1,nrow(data),27)
  dates<-alldat[,1]
  datess<-alldat[,1]
  dates<-as.character(dates)
  dates<-as.date(dates)
  days<-date.mdy(dates,weekday=TRUE)$weekday
  
  weekdays<-c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")
  
  daycounts<-array(NA,dim=c(i,7,1))
  
  for(ii in 1:i){
    for(iii in c(1:7)){
      daycounts[ii,iii,1]<-mean(fac[,ii][days==iii])
    }
    jpeg(filename=mypaste(name,i,ii,".jpg") ,width=400, height=250)
    plot(daycounts[ii,,],type="l",xaxt="n",ylab="",xlab=mypaste("Mean by days of week ",i,", ",ii))
    axis(1,at=c(1:7),labels=weekdays)
    dev.off()}
  
}

facplots<-function(fac,name="gkeyf",startdate="2001/6/1",enddate="2003/6/6"){
  #start on a logical date near the beginning of the data
  dates<-getdates()
  for(iii in 1:ncol(fac)){
    jpeg(filename=mypaste(name,i,iii,".jpg") ,width=700, height=250)
    par(cex=1) 
    plot(dates,fac[,iii],type="l",col="blue",xlab="",xaxt="n",ylab=mypaste(iii," of ",i," factors"))
    par(cex=.7)
    axis.Date(1,at=seq(as.Date(startdate),as.Date(enddate), "months"),
              format="%b-%y",las=2)
    dev.off()}
}

gkeyprep<-function(lambdakey,uncertainty,ambienterror="slunc.txt",datamatrix="sldat.txt"){
  
  ## with this variant of gkeying I'm using the key as a starting value for lambda, and 
  ## not using a starting value for F (sorting it out would be a pain)
  i<-min(dim(lambdakey))
  meserror<-read.table(ambienterror)
  uncertainty<-rbind(uncertainty,meserror)
  write.table(uncertainty,file="keyuncertainties.txt",row.names=FALSE,col.names=FALSE)
  
  #adds a matrix of zeros on top of the starting array F, obtained from a previous run
  
  #zeros<-array(0,dim=c(i,i))
  #gnot<-read.table(paste("testfac",i,".txt",sep=""))
  #keygnot<-rbind(zeros,gnot)
  #write.table(keygnot,file="keygnot.txt",row.names=FALSE,col.names=FALSE)
  
  key<-t(lambdakey)
  write.table(key,file=mypaste("startlam",i,".txt"),row.names=FALSE,col.names=FALSE)
  
  data<-read.table(datamatrix)
  keydata<-rbind(as.matrix(key),as.matrix(data))
  write.table(keydata,file=paste("keydat",i,".txt",sep=""),row.names=FALSE,col.names=FALSE)
}

lamplots<-function(lam,name="gkey"){
  for(iii in 1:ncol(lam)){
    jpeg(filename=mypaste(name,i,iii,".jpg") ,width=500, height=250)
    par(cex=1)
    plot(lam[,iii],ylim=c(0,max(lam[,iii])+.05),ylab="",
         xlab=mypaste("Source ",i,", ",iii),xaxt="n")
    par(cex=.6)
    axis(1,at=1:length(lam[,iii]),label=names,las=2)
    ###  makes a box around significant values
    par(cex=1)	
    for(ii in 1:nrow(lam)){
      if(lam[ii,iii]>=.01){
        rect(ii-.5,par("usr")[3],ii+.5,par("usr")[4]-.01,col="grey90",border="transparent")
      }}
    points(lam[,iii],bty="o") #redraws the points on top of the shaded boxes
    #end of box maxing
    #print(par("usr"))
    dev.off()
    
  }
}

getdates<-function(){
  #i used excel to get the dates in the right format for R to read
  alldat<-read.table("data_for_dates.txt",header=T)
  dates<-alldat[,1]
  dates<-as.Date(dates)
  return(dates)}

getQ<-function(filename){
  Q<-readLines(con=filename)
  Q<-Q[1]
  Q<-unlist(strsplit(Q,"=",fixed=TRUE))[2]
  Q<-as.numeric(Q)
  return(Q)}

runpmfplain<-function(inifilename,data,prec,F,Lambda,Q,fpeak=0,sources=6,rows=nrow(data),cols=ncol(data)){
  sname<-unlist(strsplit(inifilename,".",fixed=TRUE))[1]
  inifile<- readLines(con=inifilename)
  inifile[9] <-paste("     ",fpeak)
  inifile[38] <-paste("    30   T \"OLD    \"  2000  \" ",data,"    \"",sep="")
  inifile[39] <-paste("    31   T \"OLD    \"  2000  \" ",prec,"                                       \"",sep="")
  inifile[44]<-paste(" 36   F \"REPLACE\"  2000  \" ",F,"    \"",sep="")
  inifile[45]<-paste(" 37   F \"REPLACE\"  2000  \" ",Lambda,"    \"",sep="")
  inifile[46]<-paste(" 38   F \"REPLACE\"  2000  \" ",Q,"    \"",sep="")
  #inifile[40]<-paste(" 32   T \"OLD\"  2000  \" startlam",sources,".txt    \"",sep="")
  inifile[7]<-paste("          ",rows,"   ",cols,"    ",sources,"    1",sep="")
  writeLines(inifile,con=inifilename)
  
  mywrite.table(data,file="pmf_data.txt",quote=FALSE,na="-999.9")
  mywrite.table(prec,file="pmf_unc.txt",quote=FALSE,na="1")
  
  system(paste("pmf2wtst.exe ",sname,sep=""),show.output.on.console = TRUE,invisible=TRUE)
}



###Handy functions

mypaste<-function (..., sep = "", collapse = NULL) 
{
    args <- list(...)
    if (length(args) == 0) 
        if (length(collapse) == 0) 
            character(0)
        else ""
    else {
        for (i in seq(along = args)) args[[i]] <- as.character(args[[i]])
        .Internal(paste(args, sep, collapse))
    }
}

clip<-function(x){
write.table(x, file = "clipboard", sep = "\t")}

clipclean<-function(x){
write.table(x, file = "clipboard", sep = "\t",row.names=FALSE,col.names=FALSE)}



mywrite.table<-function (x, file = "", append = FALSE, quote = TRUE, sep = " ", 
    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
    qmethod = c("escape", "double")) 
{
    qmethod <- match.arg(qmethod)
    if (!is.data.frame(x)) 
        x <- data.frame(x)
    nocols <- NCOL(x) == 0
    if (nocols) 
        quote <- FALSE
    if (is.logical(quote) && quote) 
        quote <- which(unlist(lapply(x, function(x) is.character(x) || 
            is.factor(x))))
    if (dec != ".") {
        num <- which(unlist(lapply(x, function(x) is.double(x) || 
            is.complex(x))))
        if (length(num)) 
            x[num] <- lapply(x[num], function(z) gsub("\\.", 
                ",", as.character(z)))
    }
    cmplx <- sapply(x, is.complex)
    if (any(cmplx) && !all(cmplx)) 
        x[cmplx] <- lapply(x[cmplx], as.character)
    x <- as.matrix(x)
    if (!nocols) {
        i <- is.na(x)
        if (any(i)) 
            x[i] <- na
    }
    p <- ncol(x)
    d <- dimnames(x)
    if (is.logical(quote)) 
        quote <- if (quote) 
            1:p
        else NULL
    else if (is.numeric(quote)) {
        if (any(quote < 1 | quote > p)) 
            stop(paste("invalid numbers in", sQuote("quote")))
    }
    else stop(paste("invalid", sQuote("quote"), "specification"))
    rn <- FALSE
    if (is.logical(row.names)) {
        if (row.names) {
            x <- cbind(d[[1]], x)
            rn <- TRUE
        }
    }
    else {
        row.names <- as.character(row.names)
        if (length(row.names) == nrow(x)) 
            x <- cbind(row.names, x)
        else stop(paste("invalid", sQuote("row.names"), "specification"))
    }
    if (!is.null(quote) && (p < ncol(x))) 
        quote <- c(0, quote) + 1
    if (is.logical(col.names)) 
        col.names <- if (is.na(col.names) && rn) 
            c("", d[[2]])
        else if (col.names) 
            d[[2]]
        else NULL
    else {
        col.names <- as.character(col.names)
        if (length(col.names) != p) 
            stop(paste("invalid", sQuote("col.names"), "specification"))
    }
    if (file == "") 
        file <- stdout()
    else if (is.character(file)) {
        file <- file(file, ifelse(append, "a", "w"))
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) 
        stop(paste("argument", sQuote("file"), "must be a character string or connection"))
    qstring <- switch(qmethod, escape = "\\\\\"", double = "\"\"")
    if (!is.null(col.names)) {
        if (append) 
            warning("appending column names to file")
        if (!is.null(quote)) 
            col.names <- paste("\"", gsub("\"", qstring, col.names), 
                "\"", sep = "")
        writeLines(paste(col.names, collapse = sep), file, sep = eol)
    }
    if (NROW(x) == 0) 
        return(invisible(x))
    for (i in quote) x[, i] <- paste("\"", gsub("\"", qstring, 
        as.character(x[, i])), "\"", sep = "")
    if (ncol(x)) 
        writeLines(paste(c(t(x)), c(rep.int(sep, ncol(x) - 1), 
            eol), sep = "", collapse = ""), file, sep = "")
    else cat(eol, file = file)
}


sort<-function(lam,fac,truef){    
tfac<-truef
k<-min(dim(fac))      
n<-max(dim(fac))                         
#requires package combinats
library(combinat)
sortmatrix<-matrix(unlist(permn(c(1:k))),nrow=gamma(k+1),ncol=k,byrow=T)
highest<-fac
highestlam<-lam
mse<-sum((tfac - fac)^2)/(n*k)
for(i in 1:nrow(sortmatrix)){               
prm<-sortmatrix[i,]
testmatrix<-fac[,prm]
testlam<-lam[,prm]
msenew<-sum((tfac - testmatrix)^2)/(n*k)
if (msenew<mse)
{
mse<-msenew
highest<-testmatrix
highestlam<-testlam
last<-i
}}
fac<-highest
lam<-highestlam
result<-list(fac,lam)
print(paste(last,"th iteration",sep=""))
names(result)<-c("fac","lam")
return(result)
}


AAE<-function(x,true){
if(dim(x)[1]<dim(x)[2]){
x<-t(x);true<-t(true)
}
AAEval <- sum(apply(sqrt((x - true)^2), 2, sum)/dim(x)[1])
return(AAEval)
}



log.entry<-function(x,filename="c:/work/gkeyinglog.txt"){
write.table(date(),file=filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(x,file=filename,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}




