## gets the middle window of h12 data

midsection <- c(150001,350000) #For analysis
h12NeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12Neutral"
h12SelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12SelectedP"


if (substr(h12NeutralFolder,nchar(h12NeutralFolder),nchar(h12NeutralFolder))!="/"){
  h12NeutralFolder <- paste0(h12NeutralFolder,"/")
}
if (substr(h12SelectedFolder,nchar(h12SelectedFolder),nchar(h12SelectedFolder))!="/"){
  h12SelectedFolder <- paste0(h12SelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(h12NeutralFolder))
##Get h12 max matrix
h12MaxMatrixP<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.table(paste0(h12NeutralFolder,"h12_",sim,".txt"),col.names = c("SNPpos","h12"))
  h12MaxMatrixP[sim,1]<-0
  h12MaxMatrixP[sim,2]<-max(simDF$h12[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm = TRUE)
}
#Selected
for (sim in 1:nSims){
  simDF <- read.table(paste0(h12SelectedFolder,"h12_",sim,".txt"),col.names = c("SNPpos","h12"))
  h12MaxMatrixP[nSims+sim,1]<-1
  h12MaxMatrixP[nSims+sim,2]<-max(simDF$h12[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm=TRUE)
}

write(t(h12MaxMatrixP),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12MaxMatrixP.txt",ncolumns = 2)

require(pROC)
drawROCCurve("h12P","roc curve of h12P",h12MaxMatrixP)
