## standardises nSL based on allele frequency bins


standardisensl <- function(inputFolder,outputFolder,prefix){
  if (substr(inputFolder,nchar(inputFolder),nchar(inputFolder))!="/"){
    inputFolder <- paste0(inputFolder,"/")
  }
  if (substr(outputFolder,nchar(outputFolder),nchar(outputFolder))!="/"){
    outputFolder <- paste0(outputFolder,"/")
  }
  
  ##count files in folder
  nSims <- length(list.files(inputFolder))
  #Neutral
  for (sim in 1:nSims){
    simDF <- read.delim(paste0(inputFolder,prefix,"_",sim,".txt"))
    binStats<-matrix(NA,nrow=100,ncol=3)
    for (i in 1:100){
      binStats[i,1]<- sum(simDF$DAF ==i)
      binStats[i,2]<- mean(simDF$SL[simDF$DAF ==i])
      binStats[i,3]<- sd(simDF$SL[simDF$DAF ==i])
    }
    nSL <- NULL
    for (i in 1:nrow(simDF)){
      nSL <- c(nSL,(simDF$SL[i]-binStats[simDF$DAF[i],2])/binStats[simDF$DAF[i],3])
    }
    simDF$nSL <- nSL
    #plot(simDF$SNPpos,simDF$nSL,main=paste(sim,"nSL"))
    write.table(simDF[,c(2,10)],paste0(outputFolder,prefix,"_",sim,".txt"),row.names = FALSE,quote=FALSE)
  }
}
nslNeutralFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslNeutral/"
nslNeutralStandFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslNeutralStand/"
standardisensl(nslNeutralFolder,nslNeutralStandFolder,"nslNeutral")

#nslSelectedFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelected/"
#nslSelectedStandFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedStand/"
#standardisensl(nslSelectedFolder,nslSelectedStandFolder,"nslSelect")

nslSelectedPFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedP/"
nslSelectedPStandFolder <-  "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedPStand/"
standardisensl(nslSelectedPFolder,nslSelectedPStandFolder,"nslSelect")
