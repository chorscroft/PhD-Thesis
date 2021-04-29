## gets the middle window of zalpha data

midsection <- c(150001,350000) #For analysis
zalphaNeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaNeutral"
zalphaSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaSelectedP"


if (substr(zalphaNeutralFolder,nchar(zalphaNeutralFolder),nchar(zalphaNeutralFolder))!="/"){
  zalphaNeutralFolder <- paste0(zalphaNeutralFolder,"/")
}
if (substr(zalphaSelectedFolder,nchar(zalphaSelectedFolder),nchar(zalphaSelectedFolder))!="/"){
  zalphaSelectedFolder <- paste0(zalphaSelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(zalphaNeutralFolder))
##Get zalpha max matrix
zalphaMaxMatrixP<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.table(paste0(zalphaNeutralFolder,"zalpha_",sim,".txt"),col.names = c("SNPpos","zalpha"))
  zalphaMaxMatrixP[sim,1]<-0
  zalphaMaxMatrixP[sim,2]<-max(simDF$zalpha[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm = TRUE)
}
#Selected
for (sim in 1:nSims){
  simDF <- read.table(paste0(zalphaSelectedFolder,"zalpha_",sim,".txt"),col.names = c("SNPpos","zalpha"))
  zalphaMaxMatrixP[nSims+sim,1]<-1
  zalphaMaxMatrixP[nSims+sim,2]<-max(simDF$zalpha[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm=TRUE)
}

write(t(zalphaMaxMatrixP),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaMaxMatrixP.txt",ncolumns = 2)

require(pROC)
drawROCCurve("zalphaP","roc curve of zalpha",zalphaMaxMatrixP)