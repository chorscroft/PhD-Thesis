## gets the middle window of nSl data

midsection <- c(150001,350000) #For analysis
nslNeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslNeutralStand"
nslSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedPStand"


if (substr(nslNeutralFolder,nchar(nslNeutralFolder),nchar(nslNeutralFolder))!="/"){
  nslNeutralFolder <- paste0(nslNeutralFolder,"/")
}
if (substr(nslSelectedFolder,nchar(nslSelectedFolder),nchar(nslSelectedFolder))!="/"){
  nslSelectedFolder <- paste0(nslSelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(nslNeutralFolder))
##Get nSL max matrix
nslMaxMatrixPStand<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.delim(paste0(nslNeutralFolder,"nslNeutral_",sim,".txt"),sep=" ")
  simDF$nSL <- abs(simDF$nSL)
  nslMaxMatrixPStand[sim,1]<-0
  nslMaxMatrixPStand[sim,2]<-max(simDF$nSL[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]])
}
#Selected
for (sim in 1:nSims){
  simDF <- read.delim(paste0(nslSelectedFolder,"nslSelect_",sim,".txt"),sep=" ")
  simDF$nSL <- abs(simDF$nSL)
  nslMaxMatrixPStand[nSims+sim,1]<-1
  nslMaxMatrixPStand[nSims+sim,2]<-max(simDF$nSL[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]])
}
nslPstandroc<- roc(nslMaxMatrixPStand[,1],nslMaxMatrixPStand[,2])
plot(nslPstandroc,main="nsl roc 50%")
auc(nslPstandroc)
auc(nslPstandroc,partial.auc=c(0.95,1))/0.05

write(t(nslMaxMatrixPStand),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslMaxMatrixPStand.txt",ncolumns = 2)

require(pROC)
drawROCCurve("nslPStand","roc curve of nslPStand",nslMaxMatrixPStand)



##Get nSL MEDIAN matrix
nslMaxMatrixPStandmed<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.delim(paste0(nslNeutralFolder,"nslNeutral_",sim,".txt"),sep=" ")
  simDF$nSL <- abs(simDF$nSL)
  nslMaxMatrixPStandmed[sim,1]<-0
  nslMaxMatrixPStandmed[sim,2]<-median(simDF$nSL[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]])
}
#Selected
for (sim in 1:nSims){
  simDF <- read.delim(paste0(nslSelectedFolder,"nslSelect_",sim,".txt"),sep=" ")
  simDF$nSL <- abs(simDF$nSL)
  nslMaxMatrixPStandmed[nSims+sim,1]<-1
  nslMaxMatrixPStandmed[nSims+sim,2]<-median(simDF$nSL[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]])
}
nslPstandmedroc<- roc(nslMaxMatrixPStandmed[,1],nslMaxMatrixPStandmed[,2])
plot(nslPstandmedroc,main="nslmed roc 50%")
auc(nslPstandmedroc)
auc(nslPstandmedroc,partial.auc=c(0.95,1))/0.05
