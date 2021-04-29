## gets the middle window of tSel data

midsection <- c(150001,350000) #For analysis
tSelNeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutputT"
tSelSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputPT"


if (substr(tSelNeutralFolder,nchar(tSelNeutralFolder),nchar(tSelNeutralFolder))!="/"){
  tSelNeutralFolder <- paste0(tSelNeutralFolder,"/")
}
if (substr(tSelSelectedFolder,nchar(tSelSelectedFolder),nchar(tSelSelectedFolder))!="/"){
  tSelSelectedFolder <- paste0(tSelSelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(tSelNeutralFolder))
##Get tSel max matrix
tSelMaxMatrixPT<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelNeutralFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixPT[sim,1]<-0
  tSelMaxMatrixPT[sim,2]<-max(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm = TRUE)
}
#Selected
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelSelectedFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixPT[nSims+sim,1]<-1
  tSelMaxMatrixPT[nSims+sim,2]<-max(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm=TRUE)
}

write(t(tSelMaxMatrixPT),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/tSelMaxMatrixPT.txt",ncolumns = 2)

require(pROC)
drawROCCurve("tSelPT","roc curve of tSel",tSelMaxMatrixPT)



## For AVERAGE rather than MAX
## gets the middle window of tSel data

midsection <- c(150001,350000) #For analysis
tSelNeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutput"
tSelSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutput"


if (substr(tSelNeutralFolder,nchar(tSelNeutralFolder),nchar(tSelNeutralFolder))!="/"){
  tSelNeutralFolder <- paste0(tSelNeutralFolder,"/")
}
if (substr(tSelSelectedFolder,nchar(tSelSelectedFolder),nchar(tSelSelectedFolder))!="/"){
  tSelSelectedFolder <- paste0(tSelSelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(tSelNeutralFolder))
##Get tSel max matrix
tSelMaxMatrixav<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelNeutralFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixav[sim,1]<-0
  tSelMaxMatrixav[sim,2]<-mean(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm = TRUE)
}
#Selected
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelSelectedFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixav[nSims+sim,1]<-1
  tSelMaxMatrixav[nSims+sim,2]<-mean(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm=TRUE)
}
tselrocav<- roc(tSelMaxMatrixav[,1],tSelMaxMatrixav[,2])
plot(tselrocav,main="tsel roc 90% with av")
auc(tselrocav,partial.auc=c(0.95,1))

#write(t(tSelMaxMatrix2),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/tSelMaxMatrixP.txt",ncolumns = 2)

#require(pROC)
#drawROCCurve("tSel2","roc curve of tSel",tSelMaxMatrix2)



## For MEDIAN rather than MAX
## gets the middle window of tSel data

midsection <- c(150001,350000) #For analysis
tSelNeutralFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutputT"
tSelSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputPT"


if (substr(tSelNeutralFolder,nchar(tSelNeutralFolder),nchar(tSelNeutralFolder))!="/"){
  tSelNeutralFolder <- paste0(tSelNeutralFolder,"/")
}
if (substr(tSelSelectedFolder,nchar(tSelSelectedFolder),nchar(tSelSelectedFolder))!="/"){
  tSelSelectedFolder <- paste0(tSelSelectedFolder,"/")
}


##count files in folder
nSims <- length(list.files(tSelNeutralFolder))
##Get tSel max matrix
tSelMaxMatrixTmed<-matrix(NA,ncol = 2,nrow=2*nSims)

#Neutral
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelNeutralFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixTmed[sim,1]<-0
  tSelMaxMatrixTmed[sim,2]<-median(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm = TRUE)
}
#Selected
for (sim in 1:nSims){
  simDF <- read.table(paste0(tSelSelectedFolder,"tselOutput_",sim,".txt"),col.names = c("SNPpos","tSel"))
  tSelMaxMatrixTmed[nSims+sim,1]<-1
  tSelMaxMatrixTmed[nSims+sim,2]<-median(simDF$tSel[simDF$SNPpos>=midsection[1] & simDF$SNPpos<=midsection[2]],na.rm=TRUE)
}
tselrocTmed<- roc(tSelMaxMatrixTmed[,1],tSelMaxMatrixTmed[,2])
plot(tselrocTmed,main="tselT roc")
auc(tselrocTmed,partial.auc=c(0.95,1))/0.05

#write(t(tSelMaxMatrix2),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/tSelMaxMatrixP.txt",ncolumns = 2)

#require(pROC)
#drawROCCurve("tSel2","roc curve of tSel",tSelMaxMatrix2)