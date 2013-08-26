metorg30 <- read.table("../receptor/StLouis/metorg30conc.txt",na.strings=".",
   header=T,row.names=1)
metorg30u <- read.table("../receptor/StLouis/metorg30unc.txt",na.strings=".",
   header=T,row.names=1)

completeobs <- apply(is.na(metorg30),1,sum) == 0
metorg30c <- metorg30[completeobs,]
metorg30uc <- metorg30u[completeobs,]


minvals <- apply(metorg30uc,2,min,na.rm=TRUE)
for (j in 1:ncol(metorg30uc)) metorg30uc[is.na(metorg30uc[,j]),j] <- minvals[j]

### REMOVE 18th and 19th (July 4-5, 2001) ###

metorg30c <- metorg30c[c(-18,-19),]
metorg30uc <- metorg30uc[c(-18,-19),]

dim(metorg30c)
dim(metorg30uc)

write(t(metorg30c),file="metorg30c.txt",ncol=ncol(metorg30c))
write(t(metorg30uc),file="metorg30uc.txt",ncol=ncol(metorg30uc))

metorg30names <- c( "Na","Mg","Al","Si","Ph","S","Cl","K","Ca",
  "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","As","Se",
"Br","Rb","Sr","Y","Zr","Mo","Pd","Ag","Cd","In","Sn","Sb",
"Ba","La","Au","Hg","Tl","Pb","U","OC","EC","SO","NO",
"Fluoranthene","Pyrene","benzoaanthracene","chrysene",
"benzobKfluoranthene","benzoepyrene","Benzoapyrene",
"Indeno123cdpyrene","Benzoghiperylene","Dibenzoahanthracene",
"X222930trisnorneohopane","X17ah21bh29Norhopane",
"X17aH21bhhopane","X22s17aH21bH30Homohopane",
"X22R17aH21bH30Homohopane","X22S17aH21bH30Bishomohopane",
"X22R17aH21bH30bishomohopane","X20Rabbsitostane",
"X20Sabbsitostane","tetracosane","Pentacosane","Hexacosane",
"Heptacosane","Octacosane","Nonacosane","Triacontane",
"Hentriacontane","Dotriacontane","Tritriacontane",
"Tetratriacontane")

metorg30dates <- as.Date(rownames(metorg30c),"%d%b%Y")

numfacs <- 9

runpmfplain("default_nogkey.ini","metorg30c.txt","metorg30uc.txt",mypaste("testfac",numfacs,".txt"),
   mypaste("testlam",numfacs,".txt"),
   mypaste("testq",numfacs,".txt"),fpeak= 0,numfacs,nrow(metorg30c),ncol(metorg30c))
 lam <- read.table(mypaste("testlam",numfacs,".txt"))
 fac <- read.table(mypaste("testfac",numfacs,".txt"))
    
lamplots2(lam,"gkey",names=metorg30names,facnames=1:9) 
#lamplotslog(lam,"gkey",names=mynames) 
facplots3(fac,"gkeyf",metorg30dates,1:9)
#weekendplots2(i,"gkeyw",newdates)
weekendplots3(numfacs,metorg30dates,facname=fac,"gkeyw")
system(mypaste("multipanel_gkey",numfacs,".bat"),invisible=TRUE) 

lamplots2ps(lam,"gkey",names=metorg30names,1:9) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)

lamplotslogps(lam,"gkey",names=metorg30names,1:9) 
system("multipanel_gkey9profspdf.bat",invisible=TRUE)

facplots3ps(fac,"gkeyf",metorg30dates,1:9)
system("multipanel_gkey9contspdf.bat",invisible=TRUE)



lambda <- matrix(NA,



