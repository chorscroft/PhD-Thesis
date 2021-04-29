
##  zalpha function
##  requires:
##    tped file  (text string giving location)
##    window size (number of kb both sides of the equation)



zalpha <- function(tPed, windowSize){
  
  ##read in tPed file
  tPedDF <- gettPed(tPed)
  tPedDF <- removeRandomSNPs(tPedDF)
  nSNPs <- nrow(tPedDF)
  
  ##initialise zaplha matrix
  zAlphaMat <- matrix(NA,nrow=nSNPs,ncol=2)  
  
  ##initialise r2 matrix
  r2matrix <- matrix(-1,nrow=nSNPs,ncol=nSNPs)
  
  ##loop over all SNPs
  for (locus in 2:(nSNPs-1)){
    
    print.default(locus)
    
    ##ASSUMES TPED FILE IN BP ORDER SMALL TO BIG
    #add all r2 in set L 
    pairLocus <- locus - 1
    sumr2L <- 0
    countr2L <- 0
    countL <- 0
    bpLocus <- tPedDF$bp[locus]
    zAlphaMat[locus,1]<-bpLocus
    while (pairLocus > 0 && tPedDF$bp[pairLocus] >= bpLocus-windowSize*1000/2){
      countL <- countL + 1
      r2locus <- pairLocus
      while (r2locus < locus){
        curr2matrix <- r2matrix[pairLocus,r2locus]
        if (is.na(curr2matrix) == FALSE && curr2matrix==-1){
          thisr2 <- r2value(hillAlgorithm(create3By3Matrix(tPedDF,r2locus,pairLocus)))
          r2matrix[pairLocus,r2locus]<-thisr2
        } else {
          thisr2 <- curr2matrix
        }
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
    while (pairLocus <= nSNPs && tPedDF$bp[pairLocus] <= bpLocus+windowSize*1000/2){
      countR <- countR + 1
      r2locus <- pairLocus
      while (r2locus > locus){
        curr2matrix<-r2matrix[r2locus,pairLocus]
        if (is.na(curr2matrix) == FALSE && curr2matrix==-1){
          thisr2 <- r2value(hillAlgorithm(create3By3Matrix(tPedDF,r2locus,pairLocus)))
          r2matrix[r2locus,pairLocus]<-thisr2
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
  return(zAlphaMat)
}

