EVold <- function(C,VC,A,VA,P,toli=.01,maxit=50,update.A=TRUE,diagmat=FALSE)
## EVold can be used to get solutions to the CMB equations
##    C is the p by n matrix of measured mixtures, where n is the number of different measurements and
##         and p is the number of species/features (n can be 1 here if there is only one p-dimensional 
##         mixture measurement)
##    VC is the p by n matrix of measurement uncertainties (or measurement error standard deviations) for 
##         each element of C
##    A is the p by k matrix of the k profiles, with the first column being the p-dimensional profile for the 
##         first source, etc.
##    VA is the p by k matrix of profile uncertainties (or profile error standard deviations) for each 
##         element of A  
##    P is an n by n correlation matrix.  For the Effective Variance solution, use P = diag(n)
##    toli specifies max tolerance for checking convergence; default is to stop iterating when relative change in size of estimates
##         is changing by less than 1% (i.e., toli = .01)
##    maxit is maximum number of allowed iterations (default = 50)
##    update.A is a boolean argument indicating whether or not the A matrix should be iteratively update.  For the Effective Variance solution
##         use update.A=FALSE (default is TRUE)
##    diagmat is a boolean argument indicating where or not covariance matrix for the data should be treated as a diagonal matrix.  For the 
##         effective variance solution, use diagmat=TRUE (default is FALSE)
## 
##  NOTE: For the standard effective variance solution, use P=diag(nrow(C)), update.A=FALSE, diagmat=TRUE
##  NOTE: If VA is a matrix of zeroes, you get the WLS solution.  If VC and VA are both matrices of zeroes, you get the OLS solution.
##
{
	convergei <- FALSE
	iter <- 0
	p <- nrow(C)
	n <- ncol(C)
	k <- ncol(A)
	cat("\n Iteration  MaxError \n")
	cat(" __________________ \n")
	oldA <- A
        oldS <- matrix(0,k,n)
	if (diagmat==TRUE) oldVe <- diag(diag( diag(c(VC)) + kronecker( t(oldS) , diag(p) ) %*% diag(c(VA)) %*% kronecker( (oldS) , diag(p) ) ))
	if (diagmat==FALSE) oldVe <- diag(c(VC)) + kronecker( t(oldS) , diag(p) ) %*% diag(c(VA)) %*% kronecker( (oldS) , diag(p) ) 

	for (iter in 1:maxit)
	{
		if (convergei == FALSE)
		{
		   newS <- oldS + matrix( solve( kronecker(diag(n),t(oldA)) %*% solve(oldVe) %*% kronecker(diag(n),oldA) ) %*%
		        kronecker(diag(n),t(oldA)) %*% solve(oldVe) %*% c(C - A %*% oldS) , k , n)
		   if (update.A == TRUE) 
		   {
			   newA <- A + matrix( diag(c(VA)) %*% kronecker( oldS , diag(p) ) %*% solve(oldVe) %*% 
			       ( diag(p*n) - kronecker( diag(n) , (oldA) ) %*% 
					solve(kronecker( diag(n) , t(oldA) ) %*% solve(oldVe) %*% kronecker( diag(n) , (oldA) )) %*%
					kronecker( diag(n) , t(oldA) ) %*% solve(oldVe) ) %*% c(  C - A %*% oldS ) , p , k )

			}

	       if (diagmat==TRUE) newVe <- diag(diag( diag(c(VC)) + kronecker( t(newS) , diag(p) ) %*%
                     diag(c(VA)) %*% kronecker( (newS) , diag(p) ) ))
	       if (diagmat==FALSE) newVe <- diag(c(VC)) + kronecker( t(newS) , diag(p) ) %*% diag(c(VA)) %*%
                    kronecker( (newS) , diag(p) ) 
		}

       cat(iter, "    ", round(max(abs( (newS - oldS)/oldS ) ),4) ,"\n" )
       if (max( abs( (newS - oldS)/oldS ) ) < toli) 
		{
			convergei <- TRUE
			if (update.A == TRUE) oldA <- newA
			cat("Iteration converged.\n")
			break;
		}
       else
		{
			oldS <- newS
			if (update.A == TRUE) oldA <- newA
			oldVe <- newVe
		}
       }
   diagVe <- diag(newVe)
   SES <- sqrt( diag( solve(kronecker( diag(n),t(oldA)) %*% solve(newVe) %*% kronecker(diag(n),oldA)) ) )
   newSbar <- apply(newS,1,mean)

        
   SESbar <- sqrt( 1/(n^2) * apply( matrix(SES^2,k,n) ,1,sum ) )
#   if (update.A == TRUE) return(newS,SES,newSbar,SESbar,oldA)
#   if (update.A == FALSE) return(newS,SES,newSbar,SESbar)
   if (update.A == TRUE) return(list(newS,SES,newSbar,SESbar,oldA))
   if (update.A == FALSE) return(list(newS,matrix(SES,k,n,byrow=FALSE),newSbar,SESbar))
}


C <- exp(t(matrix(rnorm(200),20,10)))
VC <- matrix(0,10,20)
A <- matrix(rnorm(50,10,1),10,5)
VA <- matrix(0,10,5)
EVold(C,VC,A,VA,matrix(0,10,10),update.A=F,diagmat=F)
EVold(C[,1:2],VC[,1:2],A,VA,diag(nrow(C)),update.A=FALSE,diagmat=TRUE)

#############################
### Fleet1 Estimation ###
#############################

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\Fleets")
#load("yfleet1profiles")
#load("yfleet1VA")
#memory.limit(20000)

days<-1
ysims<-matrix(NA,93,days)
cv<-.2
for (i in 1:93)  ###generating y's distributed lognormally,  profiles for fleet1 over days days
{
  for(j in 1:days)
  {
    ysims[i,j]<-exp(rnorm(1,mean=log(yfleet1[i])-(.5*log(((cv)^2)+1)),sd=sqrt(log(((cv)^2)+1))))
  }
}


iter<-days
simev<-matrix(NA,iter,9)
simevmeans<-matrix(NA,1,9)
colnames(simev)<-c("fgas","fsmoker","fdiesel","sesgas","sessmoker","sesdiesel","msegas","msesmoker","msediesel")
colnames(simevmeans)<-c("fgas","fsmoker","fdiesel","sesgas","sessmoker","sesdiesel","msegas","msesmoker","msediesel")
for(i in 1:iter)
{
  C<-as.matrix(ysims[,i])
  VC<-C*.2
  A<-yfleet1profiles
  VA<-yfleet1VA
  effvar<-EVold(C,VC,A,VA,diag(nrow(C)),update.A=FALSE,diagmat=TRUE)
  simev[i,1]<-effvar[[1]][1,1]
  simev[i,2]<-effvar[[1]][2,1]
  simev[i,3]<-effvar[[1]][3,1]
  simev[i,4]<-effvar[[2]][1,1]
  simev[i,5]<-effvar[[2]][2,1]
  simev[i,6]<-effvar[[2]][3,1]
  simev[i,7]<-(simev[i,1]-15)^2
  simev[i,8]<-(simev[i,2]-.937)^2
  simev[i,9]<-(simev[i,3]-2.8215)^2
}

simevmeans[1,1]<-mean(simev[,1])
simevmeans[1,2]<-mean(simev[,2])
simevmeans[1,3]<-mean(simev[,3])
simevmeans[1,4]<-mean(simev[,4])
simevmeans[1,5]<-mean(simev[,5])
simevmeans[1,6]<-mean(simev[,6])
simevmeans[1,7]<-mean(simev[,7])
simevmeans[1,8]<-mean(simev[,8])
simevmeans[1,9]<-mean(simev[,9])

simevmeans

setwd("C:\\Users\\David Coats\\Documents\\Research for Dr. Christensen\\Source Aportionment\\EV simulations")

load("sim1100000cv30")
load("sim1cv30")
#save(sim1means,file="sim1cv30")
#save(sim1,file="sim1100000cv30")
sim1means
cvestgas<-sqrt(sim1means[7])/15
cvestgas
cvestsmoke<-sqrt(sim1means[8])/.937
cvestsmoke
cvestdiesel<-sqrt(sim1means[9])/2.8215
cvestdiesel

#load("sim1100000cv20")
load("sim1cv20")
#save(sim1means,file="sim1cv20")
#save(sim1,file="sim1100000cv20")
sim1means
cvestgas<-sqrt(sim1means[7])/15
cvestgas
cvestsmoke<-sqrt(sim1means[8])/.937
cvestsmoke
cvestdiesel<-sqrt(sim1means[9])/2.8215
cvestdiesel





##################################################################################################
###############

iter<-50
simev<-matrix(NA,iter,9)
simevmeans<-matrix(NA,1,9)
colnames(simev)<-c("fgas","fsmoker","fdiesel","sesgas","sessmoker","sesdiesel","msegas","msesmoker","msediesel")
colnames(simevmeans)<-c("fgas","fsmoker","fdiesel","sesgas","sessmoker","sesdiesel","msegas","msesmoker","msediesel")
for(i in 1:iter)
{
  C<-as.matrix(ysims[,i])
  VC<-C*cv
  A<-yfleet1profiles
  VA<-yfleet1VA
  effvar<-EVold(C,VC,A,VA,diag(nrow(C)),update.A=FALSE,diagmat=TRUE)
  simev[i,1]<-effvar[[1]][1,1]
  simev[i,2]<-effvar[[1]][2,1]
  simev[i,3]<-effvar[[1]][3,1]
  simev[i,4]<-effvar[[2]][1,1]
  simev[i,5]<-effvar[[2]][2,1]
  simev[i,6]<-effvar[[2]][3,1]
  simev[i,7]<-(simev[i,1]-fpre[i,1])^2
  simev[i,8]<-(simev[i,2]-fpre[i,2])^2
  simev[i,9]<-(simev[i,3]-fpre[i,3])^2
}

simevmeans[1,1]<-mean(simev[,1])
simevmeans[1,2]<-mean(simev[,2])
simevmeans[1,3]<-mean(simev[,3])
simevmeans[1,4]<-mean(simev[,4])
simevmeans[1,5]<-mean(simev[,5])
simevmeans[1,6]<-mean(simev[,6])
simevmeans[1,7]<-mean(simev[,7])
simevmeans[1,8]<-mean(simev[,8])
simevmeans[1,9]<-mean(simev[,9])

simevmeans