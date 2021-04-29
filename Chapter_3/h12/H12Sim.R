## Reads in the data from the simulation and runs H12 over it
##  Returns a List with an element for each simulation
##  Element is a matrix where
##    column 1 is the bp location
##    column 2 is the H12 value

## Input:
##  Simulation Output folder
##  Window Size in Snps

h12Sim <- function (simFolder, windowSize, stepSize, outFile){
  
  ##count files in folder
  nSims <- length(list.files(simFolder))/2
  ##add "/" to end of folder name if not already
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  
  h12SimList <- vector("list",nSims)
  
  for (sim in 1:nSims){
    print.default(paste("sim",sim))
    simBPs <- unlist(read.table(paste0(simFolder,"simBPs_",sim,".txt")),use.names = FALSE)  
    simMatrix <-data.matrix(read.table(paste0(simFolder,"simMatrix_",sim,".txt")))
    colnames(simMatrix)<-NULL
    
    nSNPs <- length(simBPs)
    nChrm <- nrow(simMatrix)
    
    ##initialise H12 matrix
    h12Mat <- matrix(NA,nrow=nSNPs,ncol=2)  
    h12Mat[,1] <- simBPs
    
    if(nSNPs>=windowSize){
      for (snp in seq(windowSize/2,nSNPs-windowSize/2,stepSize)){
        haplotypeMat <- matrix(simMatrix[1,(snp-(windowSize/2)+1):(snp+(windowSize/2))],nrow=1)
        haplotypeCount <- 1
        for(chrm in 2:nChrm){
          testVec<-simMatrix[chrm,(snp-(windowSize/2)+1):(snp+(windowSize/2))]
          foundMatch <- FALSE
          for (check in 1:length(haplotypeCount)){
            if (isTRUE(all.equal(testVec,haplotypeMat[check,]))){
              foundMatch <- TRUE
              haplotypeCount[check] <- haplotypeCount[check]+1
              break
            }
          }
          if (foundMatch==FALSE){
            haplotypeMat <- rbind(haplotypeMat,testVec)
            haplotypeCount <- c(haplotypeCount,1)
          }
        }
        ##now calculate proportions
        nhap <- length(haplotypeCount)
        H1<-0
        mostCom <- 0
        secondMost <- 0
        for (i in 1:nhap){
          H1 <- H1 + (haplotypeCount[i]/nChrm)^2
          if(haplotypeCount[i]>secondMost){
            if(haplotypeCount[i]>mostCom){
              secondMost <- mostCom
              mostCom<-haplotypeCount[i]
            } else {
              secondMost <- haplotypeCount[i]
            }
          }
        }
        H12 <- H1 + 2*(mostCom/nChrm)*(secondMost/nChrm)
        h12Mat[snp,2] <- H12
      }
    }
    ##save h12 mat to list
    h12SimList[[sim]] <- h12Mat
    write.table(h12Mat,file=paste0(outFile,"h12_",sim,".txt"),col.name=FALSE,row.names=FALSE)
  }
  return(h12SimList)
}


##gets the max value for the mid section of a simulation run
geth12Max<-function(simMatrix,midsection){
  h12max <- max(simMatrix[simMatrix[,1]>=midsection[1] & simMatrix[,1]<=midsection[2],2],na.rm = TRUE)
  return(h12max)
}

##gets the matrix of maximum h12 values 
##column 1 indictates if value is from the neutral (0) or selected (1) simulation
##column 2 is the maximum h12 value
geth12MaxMatrix<-function(h12Neutral,h12Select,midsection){
  h12Matrix<-matrix(NA,ncol=2,nrow=length(h12Neutral)+length(h12Select))
  #neutral
  for (i in 1:(length(h12Neutral))){
    h12Matrix[i,1] <- 0
    h12Matrix[i,2] <-geth12Max(h12Neutral[[i]],midsection)
  }
  #selected
  for (i in 1:(length(h12Select))){
    h12Matrix[length(h12Neutral)+i,1] <- 1
    h12Matrix[length(h12Neutral)+i,2] <-geth12Max(h12Select[[i]],midsection)
  }  
  return(h12Matrix)
}





