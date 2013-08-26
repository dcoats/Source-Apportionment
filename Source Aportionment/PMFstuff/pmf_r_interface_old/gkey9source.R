#########################
#
#    Simple PMF-R Interface
#
#            Jeff Lingwall June 05
#
#########################

######
#
#   Pre-reqs
#
######
setwd("c:/my documents/research/pmf_r_interface")
source("jeff's functions.R")
#load(".Rdata")

alldat<-read.table("sldat_naremoved.txt",header=T)

data<-cbind(alldat$NAXC,alldat$MGXC,alldat$ALXC,alldat$SIXC,alldat$PHXC,
		alldat$SUXC,alldat$CLXC,alldat$KPXC,alldat$CAXC,alldat$TIXC,
		alldat$VAXC,alldat$CRXC,alldat$MNXC,alldat$FEXC,alldat$COXC,
		alldat$NIXC,alldat$CUXC,alldat$ZNXC,alldat$GAXC,alldat$ASXC,
		alldat$SEXC,alldat$BRXC,alldat$RBXC,alldat$SRXC,alldat$YTXC,
		alldat$ZRXC,alldat$MOXC,alldat$PDXC,alldat$AGXC,alldat$CDXC,
		alldat$INXC,alldat$SNXC,alldat$SBXC,alldat$BAXC,
		alldat$LAXC,alldat$AUXC,alldat$HGXC,alldat$TLXC,alldat$PBXC,
		alldat$URXC,alldat$OC,alldat$EC,alldat$Sulfate,alldat$Nitrate)

write.table(data,file="sldat.txt",quote=F,na="-999.9",row.names=FALSE,col.names=FALSE)

names<-c("Na","Mg","Al","Si","Ph",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO","NO")

names3<-c("Na","Mg","Al","Si","P",
		"S","Cl","K","Ca","Ti",
		"V","Cr","Mn","Fe","Co",
		"Ni","Cu","Zn","Ga","As",
		"Se","Br","Rb","Sr","Y",
		"Zr","Mo","Pd","Ag","Cd",
		"In","Sn","Sb","Ba",
		"La","Au","Hg","Tl","Pb",
		"U","OC","EC","SO4","NO3")


########
#
#   START HERE . . .
#
########


##create a target profile matrix, lambda
#note-- I assume lambda is p by k

#lams<-read.table("lam7.txt")  #the one I was working with
#smelter<-apply(lams7[,c(1,3)],1,mean)   #smelter

#lambda<-cbind(lams7[,c(2,4,6)],smelter,lams7[,c(5,7)])
#lambda<-cbind(lam11[,1],lam9[,2],lam9[,3],lam9[,5],
#              lam6[,3],lam10[,6],lam6[,4],lam10[,10],lam10[,7])

# The following profiles include Cu, Zn, Pb smelters that are measured before the SO are pulled out
lambda<-cbind(
lam11[,1],  # Pb
lam12[,2],  # Cu
lam10[,3],   # fireworks
lam6[,3],   # Summ
lam13[,4],  # Zn
lam10[,7],   # steel
lam6[,4],   # mobile
lam10[,6],  # soil 
lam10[,10] # winter
)

# The following profiles include Cu, Zn, Pb smelters that are measured before the SO are pulled out
lambda<-cbind(
lam13[,1],  # Pb
lam13[,2],  # Cu
lam10[,3],   # fireworks
lam6[,3],   # Summ
lam13[,4],  # Zn
lam10[,7],   # steel
lam6[,4],   # mobile
lam10[,6],  # soil 
lam10[,10] # winter
)

i<-min(dim(lambda))
k<-max(dim(lambda))
n<-nrow(sldat2)

  write(t(slunc2),file="slunc2.txt",ncol=ncol(slunc2))
  write(t(sldat2),file="sldat2.txt",ncol=ncol(sldat2))
datesnew <- getdates()
firstofmonth <- c(
  "2001-06-01","2001-07-01","2001-08-01","2001-09-01","2001-10-01","2001-11-01","2001-12-01",
  "2002-01-01","2002-02-01","2002-03-01","2002-04-01","2002-05-01","2002-06-01",
  "2002-07-01","2002-08-01","2002-09-01","2002-10-01","2002-11-01","2002-12-01",
  "2003-01-01","2003-02-01","2003-03-01","2003-04-01","2003-05-01")


#testing
#uncertainty<-array(NA,dim=c(i,k))
#uncertainty[1,]<-.001           #Pb smelter  
#uncertainty[2,]<-.001           #Cu smelter 
#uncertainty[3,]<-.001           #fireworks 
#uncertainty[4,]<-.001            #Zn smelter 
#uncertainty[5,]<-.1             #summer secondary  
#uncertainty[6,]<-.001             #soil    
#uncertainty[7,]<-1             #mobile    
#uncertainty[8,]<-.1             #winter secondary    
#uncertainty[9,]<-.1             #steel mill

uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-10           #Pb smelter  
uncertainty[2,]<-10           #Cu smelter 
uncertainty[3,]<-.001           #fireworks 
uncertainty[4,]<-10           #Zn smelter 
uncertainty[5,]<-10             #summer secondary  
uncertainty[6,]<-10             #soil    
uncertainty[7,]<-100             #mobile    
uncertainty[8,]<-10             #winter secondary    
uncertainty[9,]<-10             #steel mill

uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-100           #Pb smelter  
uncertainty[2,]<-100          #Cu smelter 
uncertainty[3,]<-.001           #fireworks 
uncertainty[4,]<-100           #Zn smelter 
uncertainty[5,]<-100             #summer secondary  
uncertainty[6,]<-100             #soil    
uncertainty[7,]<-1000             #mobile    
uncertainty[8,]<-100             #winter secondary    
uncertainty[9,]<-100             #steel mill



uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-.01           #Pb smelter   .01
uncertainty[2,]<-.01          #Cu smelter    .01
uncertainty[3,]<-.001           #fireworks    .001
uncertainty[4,]<-1             #summer secondary  1
uncertainty[5,]<-.01           #Zn smelter   .01
uncertainty[6,]<-.01             #steel mill    1
uncertainty[7,]<-10             #mobile        10
uncertainty[8,]<-1             #soil          1
uncertainty[9,]<-1             #winter secondary 1   

gkeyprep(lambda,uncertainty,"slunc2.txt","sldat2.txt")     #gets ready to gkey

runpmf("default_gkey_9source.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),i,n+i,k)       # ini file must be set up for gkeying, 
                                            #   if no gkeying is used, change "keydat" to "sldat"
#lam8gkey <- lam <- read.table(mypaste("testlam",i,".txt"))
#fac <- read.table(mypaste("testfac",i,".txt"))
#fac8gkey <- fac <- fac[-c(1:i),]    #if gkeying is used

lam9gkey <- lam <- read.table(mypaste("testlam",i,".txt"))
fac <- read.table(mypaste("testfac",i,".txt"))
fac9gkey <- fac <- fac[-c(1:i),]    #if gkeying is used
round(apply(fac9gkey,2,mean),2)
startvalfac <- fac9gkeynostart <- fac9gkey
lam9gkeynostart <- lam9gkey

lamplots2(lam9gkeynostart,"gkey",names=mynames,facnames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac9gkeynostart,"gkeyf",newdates,facnames)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac9gkeynostart,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


## PLOTS FOR ENVIRONMETRICS PAPER ##
facnames <- c("Lead\nSmelter","Fireworks","Copper\nSmelter","Summer\nSecondary","Zinc\nSmelter",
              "Soil","Steel\nMill","Winter\nSecondary","Mobile")
reorderfacs <- c(4,8,9,7,6,5,2,3,1)
lamplots2ps(lam9gkeynostart[,reorderfacs],"gkey",names=mynames,facnames[reorderfacs]) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)

lamplotslogps(lam9gkeynostart[,reorderfacs],"gkey",names=mynames,facnames[reorderfacs]) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)

facplots3ps(fac9gkeynostart[,reorderfacs],"gkeyf",newdates,facnames[reorderfacs])
system("multipanel_gkey9contspdf.bat",invisible=TRUE)

#####
#####
#####
#####
facnames3 <- c("LeadSmelter","Fireworks","CopperSmelter","SummerSecondary","ZincSmelter",
              "Soil","SteelMill","WinterSecondary","Mobile")
lamout <- lam9gkeynostart
facout <- fac9gkeynostart
dimnames(facout)[[2]] <- facnames3
facout <- facout[,reorderfacs]
dataout <- sldat2[-c(47,48,364,365),]
dimnames(dataout)[[2]] <- names3
dataout <- dataout[,apply(dataout,2,mean)>.01]

ranset <- sample(1:nrow(dataout),200)
dataout2 <- dataout[ranset[order(ranset)],]
facout2 <- facout[ranset[order(ranset)],-7]

ranset3 <- sample(1:nrow(dataout),50)
dataout3 <- dataout[ranset3[order(ranset3)],]
facout3 <- facout[ranset3[order(ranset3)],-7]


#write(t(lam9gkeynostart[,reorderfacs]),file="ChemoPaperProfiles.txt",ncol=9)
#write(t(fac9gkeynostart[,reorderfacs]),file="ChemoPaperContributions.txt",ncol=9)
write.table(facout,file="PSAContributions.csv",
            col.names=TRUE,row.names=FALSE,sep=",")
write.table(dataout,file="PSAData.csv",
            col.names=TRUE,row.names=FALSE,sep=",")
write.table(facout2,file="C:/My Documents/stat611/datasets/PSAContributions.txt",
            col.names=TRUE,row.names=FALSE)
write.table(dataout2,file="C:/My Documents/stat611/datasets/PSAData.txt",
            col.names=TRUE,row.names=FALSE)

junkc <- read.csv("PSAdata.csv",header=T)



fac9PMFout <- fac9gkeynostart[,reorderfacs]
dimnames(fac9PMFout)[[2]] <- facnames3[reorderfacs]
lam9PMFout <- lam9gkeynostart[,reorderfacs]
dimnames(lam9PMFout) <- list(metnames, facnames3[reorderfacs])

write(t(fac9PMFout),file="fac9PMFout.txt",ncol=ncol(fac9PMFout))
write(t(lam9PMFout),file="lam9PMFout.txt",ncol=ncol(lam9PMFout))

facnamesout <- facnames3[reorderfacs]



#####
#####
#####
#####





# The following profiles include Cu, Zn, Pb smelters that are measured before the SO are pulled out
lambda<-cbind(
lam11[,1],  # Pb
lam9[,3],   # fireworks
lam12[,2],  # Cu
lam13[,4],  # Zn
lam6[,3],   # Summ
lam10[,6],  # soil 
lam10[,7],   # steel
lam6[,4],   # mobile
lam10[,10] # winter
)

uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-1              #Pb smelter  
uncertainty[2,]<-.001          #fireworks 
uncertainty[3,]<-1              #Cu smelter 
uncertainty[4,]<-1              #Zn smelter 
uncertainty[5,]<-100             #summer secondary  
uncertainty[6,]<-1          #soil    
uncertainty[7,]<-1             #steel mill
uncertainty[8,]<-100             #mobile    
uncertainty[9,]<-100             #winter secondary    

uncertainty<-array(NA,dim=c(i,k))
uncertainty[1,]<-1 * c(lambda[,1])              #Pb smelter  
uncertainty[2,]<-.1 * c(lambda[,2])         #fireworks 
uncertainty[3,]<-1 * c(lambda[,3])              #Cu smelter 
uncertainty[4,]<-1 * c(lambda[,4])              #Zn smelter 
uncertainty[5,]<-1 * c(lambda[,5])             #summer secondary  
uncertainty[6,]<-1 * c(lambda[,6])          #soil    
uncertainty[7,]<-1 * c(lambda[,7])             #steel mill
uncertainty[8,]<-10 * c(lambda[,8])             #mobile    
uncertainty[9,]<-1 * c(lambda[,9])             #winter secondary    

startvalfacOLD<-cbind(
fac11[,1],  # Pb
fac9[,3],   # fireworks
fac12[,2],  # Cu
fac13[,4],  # Zn
fac6[,3],   # Summ
fac10[,6],  # soil 
fac10[,7],   # steel
fac6[,4],   # mobile
fac10[,10] # winter
)



#gnot <- rbind(matrix(0,9,9),startvalfac)
gnot <- rbind(matrix(0,9,9),startvalfacOLD)
#write(as.matrix(gnot),file="testfac9start.txt",ncol=nrow(gnot))
write(t(gnot),file="testfac9start.txt",ncol=ncol(gnot))

gkeyprep(lambda,uncertainty,"slunc2.txt","sldat2.txt")     #gets ready to gkey


runpmf("default_gkey_9source06May.ini",mypaste("keydat",i,".txt"),mypaste("testfac",i,".txt"),
   mypaste("testlam",i,".txt"),
   mypaste("testq",i,".txt"),i,n+i,k)       # ini file must be set up for gkeying, 
                                            #   if no gkeying is used, change "keydat" to "sldat"

lam9gkey <- lam <- read.table(mypaste("testlam",i,".txt"))
fac <- read.table(mypaste("testfac",i,".txt"))
fac9gkey <- fac <- fac[-c(1:i),]    #if gkeying is used
apply(fac9gkey,2,mean)

lamplots2(lam,"gkey",names=mynames,facnames) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",newdates)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(i,newdates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",i,".bat"),invisible=TRUE) 


lamfacplots(lam,fac,newdates,names=mynames,facnames)











