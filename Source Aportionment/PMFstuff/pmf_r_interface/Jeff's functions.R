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




