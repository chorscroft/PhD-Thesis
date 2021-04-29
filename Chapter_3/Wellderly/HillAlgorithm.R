## Given two loci with 2 alleles, work out the chromosomal arrangement and D
## Implementation of the algorithm in Hill (1974)
## Estimation of linkage disequilibrium in randomly mating populations

#Input a 3 by 3 matrix of observed values
#Rows are AA, Aa, aa
#Columns are BB, Bb, bb
#N <- matrix(c(57,140,101,39,224,226,3,54,156),nrow=3,byrow =TRUE)

hillAlgorithm <- function(N){
  if(is.matrix(N)==FALSE){
    #print.default("Hill Algorithm function not given a matrix as input")
    return(NA)
  } else if(nrow(N) !=3 || ncol(N) != 3){
    #print.default("Hill Algorithm function was not given a 3x3 matrix as input")
    return(NA)
  } else {  
    sumN <- sum(N)
    sumRows <- apply(N,1,sum)
    sumColumns <-apply(N,2,sum)
    
    if (sumRows[1]==sumN || sumRows[2]==sumN ||sumRows[3]==sumN
        || sumColumns[1]==sumN || sumColumns[2]==sumN ||sumColumns[3]==sumN){
      #print.default("Hill Algorithm failed as one SNP has no diversity")
      return(NA)
    }
    
    #Calculate estimate of p and q, the allele frequencies of A and B respectively
    p <- (sumRows[1]+0.5*sumRows[2])/sumN
    q <- (sumColumns[1]+0.5*sumColumns[2])/sumN
    
    #Derive X values (values needed later by calculation)
    X <- matrix(c(2*N[1,1]+N[1,2]+N[2,1],2*N[1,3]+N[1,2]+N[2,3],
                  2*N[3,1]+N[3,2]+N[2,1],2*N[3,3]+N[3,2]+N[2,3]),nrow=2,byrow=TRUE)
    
    #Initialize f11 (frequency of AB's in the population)
    f11 <- (1/(4*sumN))*(X[1,1]-X[1,2]-X[2,1]+X[2,2])+0.5-(1-p)*(1-q)
    
    #Loop until the difference between the old f11 and the new one is less than 10e-8
    itcount <- 0 #count the iterations
    maxitcount <- 15000 #number of iterations before timeout
    targetDiff <- 0.00000001
    diff <- 1
    while (diff > targetDiff && itcount < maxitcount){
      itcount <- itcount + 1
      oldf11 <- f11
      f11 <- (X[1,1]+N[2,2]*oldf11*(1-p-q+oldf11)/(oldf11*(1-p-q+oldf11)+(p-oldf11)*(q-oldf11)))/(2*sumN)
      diff <- abs(oldf11-f11)
    }
    if(itcount == maxitcount){
      #print.default("Timed Out")
    } else {
      #Return the matrix
      D <- f11 - p*q
      matrixf <- matrix(c(f11*sumN,(p*(1-q)-D)*sumN,
                          ((1-p)*q-D)*sumN,((1-p)*(1-q)+D)*sumN),nrow=2,byrow=TRUE)
      return(matrixf)
    }
  }
}
