
##  zalpha function
##  requires:
##    folder containing simulation files (text string giving location)
##    window size (number of kb both sides of the locus)
##    stepsize (in snps)
##    outFolder - where the output shoudl be saved
##  Returns a List with an element for each simulation
##  Element is a matrix where
##    column 1 is the bp location for a snp
##    column 2 is the zalpha value



zalphaSim <- function(simFolder, windowSize, stepSize, outFolder){
  
  ## count files in folder
  nSims <- length(list.files(simFolder))/2
  
  ## add "/" to end of folder name if not already
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  ## add "/" to end of folder name if not already
  if (substr(outFolder,nchar(outFolder),nchar(outFolder))!="/"){
    outFolder <- paste0(outFolder,"/")
  }
  
  ## initilise output list
  zalphaSimList <- vector("list",nSims)
  
  ## Loop over each simulation in turn
  for (sim in 1:nSims){
    ## Read in bp locations
    simBPs <- unlist(read.table(paste0(simFolder,"simBPs_",sim,".txt")),use.names = FALSE)  
    ## Read in matrix of snp status:
    ## rows are chromosomes
    ## columns are each snp
    ## 0 = ancenstral; 1 = derived
    simMatrix <-data.matrix(read.table(paste0(simFolder,"simMatrix_",sim,".txt")))
    colnames(simMatrix)<-NULL
    
    tempList <- removeRandomSNPs(simBPs,simMatrix,500)
    simBPs <- tempList[[1]]
    simMatrix<-tempList[[2]]
    rm(tempList)
    
    
    ## get the number of snps
    nSNPs <- length(simBPs)
    
    ## initialise zalpha matrix
    ## column one is the bp location
    ## column two is the zalpha value
    zAlphaMat <- matrix(NA,nrow=nSNPs,ncol=2)  
    zAlphaMat[,1] <- simBPs
    
    ## initialise r2 matrix for each snp vs each snp
    r2vec <- simBPs[simBPs >= simBPs[stepSize]-windowSize*1000/2 & simBPs <= simBPs[stepSize]+windowSize*1000/2]
    r2matrix <- matrix(-1,nrow=length(r2vec),ncol=length(r2vec))
    
    ## loop over all SNPs (ignoring the end bits)
    for (locus in seq(stepSize,(nSNPs-stepSize+1),stepSize)){
      
      #record bp location of locus
      bpLocus <- simBPs[locus]
      
      ##new r2matrix
      newr2vec <-simBPs[simBPs >= bpLocus-windowSize*1000/2 & simBPs <= bpLocus+windowSize*1000/2]
      while (r2vec[1] != newr2vec[1]){
        r2vec <- r2vec[c(-1)]
        r2matrix <- r2matrix[c(-1),c(-1)]
      }
      while (length(r2vec) !=length(newr2vec)){
        r2vec <- c(r2vec,newr2vec[length(r2vec)+1])
        r2matrix <- cbind(rbind(r2matrix,rep(-1,length(r2vec)-1)),rep(-1,length(r2vec)))
      }
      
      ## add all r2 in set L
      
      pairLocus <- locus - 1  # snp immediately to the left of current focus
      
      #initialise count
      sumr2L <- 0
      countr2L <- 0
      countL <- 0
      
      # Loop over all pairs to the left of the current locus, within the window
      while (pairLocus > 0 && simBPs[pairLocus] >= bpLocus-windowSize*1000/2){
        countL <- countL + 1
        r2locus <- pairLocus+1
        while (r2locus < locus){
          
          # read r2 value from the r2 matrix
          curr2matrix <- r2matrix[r2vec==simBPs[pairLocus],r2vec==simBPs[r2locus]][1]
          
          # if it has not yet been calculated, then calculate it
          if (is.na(curr2matrix) == FALSE && curr2matrix==-1){
            thisr2 <- r2value(create2By2MatrixSim(simMatrix,r2locus,pairLocus))
            r2matrix[r2vec==simBPs[pairLocus],r2vec==simBPs[r2locus]]<-thisr2
          } else {
            thisr2 <- curr2matrix
          }
          
          # if the r2 value is valid then add it to the sum of r2 values for L
          if (is.na(thisr2)==FALSE && thisr2 >=0){
            sumr2L <- sumr2L + thisr2
            countr2L <- countr2L+1
          }
          r2locus <- r2locus+1
        }
        pairLocus<- pairLocus-1
      }
      
      #add all r2 in set R 
      
      pairLocus <- locus + 1
      sumr2R <- 0
      countr2R <- 0
      countR <- 0
      while (pairLocus <= nSNPs && simBPs[pairLocus] <= bpLocus+windowSize*1000/2){
        countR <- countR + 1      
        r2locus <- pairLocus-1
        while (r2locus > locus){
          curr2matrix<-r2matrix[r2vec==simBPs[r2locus],r2vec==simBPs[pairLocus]][1]
          if (is.na(curr2matrix) == FALSE && curr2matrix==-1){
            thisr2 <- r2value(create2By2MatrixSim(simMatrix,r2locus,pairLocus))
            r2matrix[r2vec==simBPs[r2locus],r2vec==simBPs[pairLocus]]<-thisr2
          } else {
            thisr2 <- curr2matrix
          }        
          if (is.na(thisr2)==FALSE && thisr2 >=0){
            sumr2R <- sumr2R + thisr2
            countr2R <- countr2R+1
          }
          r2locus <- r2locus-1
        }
        pairLocus<- pairLocus+1
      }    
      if (countR >= 4 && countL >= 4 && countR*countL>=25){
        zAlphaMat[locus,2] <- (sumr2L/countr2L+sumr2R/countr2R)/2
      }
    }
    ## save zalpha matrix to output List
    write.table(zAlphaMat,file=paste0(outFolder,"zalpha_",sim,".txt"),col.name=FALSE,row.names=FALSE)
    zalphaSimList[[sim]] <- zAlphaMat
  }
  return(zalphaSimList)
}

## Creates a 2x2 matrix for a given pair of SNPs in a given simulated matrix
create2By2MatrixSim <- function(simMatrix,i,j){
  matrixf <- matrix(0,nrow = 2,ncol=2) #00 01 10 11
  for(x in 1:nrow(simMatrix)){
    matrixf[simMatrix[x,i]+1,simMatrix[x,j]+1]<-matrixf[simMatrix[x,i]+1,simMatrix[x,j]+1]+1
  }
  return(matrixf)
}

## gets the max value for the mid section of a simulation run
getZalphaMax<-function(simMatrix,midsection){
  zalphamax <- max(simMatrix[simMatrix[,1]>=midsection[1] & simMatrix[,1]<=midsection[2],2],na.rm = TRUE)
  return(zalphamax)
}

## gets the matrix of maximum zalpha values 
## column 1 indictates if value is from the neutral (0) or selected (1) simulation
## column 2 is the maximum zalpha value
getZalphaMaxMatrix<-function(zalphaNeutral,zalphaSelect,midsection){
  zalphaMatrix<-matrix(NA,ncol=2,nrow=length(zalphaNeutral)+length(zalphaSelect))
  #neutral
  for (i in 1:(length(zalphaNeutral))){
    zalphaMatrix[i,1] <- 0
    zalphaMatrix[i,2] <-getZalphaMax(zalphaNeutral[[i]],midsection)
  }
  #selected
  for (i in 1:(length(zalphaSelect))){
    zalphaMatrix[length(zalphaNeutral)+i,1] <- 1
    zalphaMatrix[length(zalphaNeutral)+i,2] <-getZalphaMax(zalphaSelect[[i]],midsection)
  }  
  return(zalphaMatrix)
}

#Removes snps randomly until the avergae distance between them is >= given bp
removeRandomSNPs <- function(simBPs,simMatrix,bp=500){
  while (averageDistSNPs(simBPs) < bp){
    #remove random SNP
    removeSNP <- sample(1:length(simBPs),1)
    simBPs <- simBPs[-c(removeSNP)]
    simMatrix <- simMatrix[,-c(removeSNP)]
  }
  return(list(simBPs,simMatrix))
}
averageDistSNPs <- function(x){
  sumDist <- 0
  for (i in 1:(length(x)-1)){
    sumDist <- sumDist + x[i+1]-x[i]
  }
  return(sumDist/(length(x)-1))
}









