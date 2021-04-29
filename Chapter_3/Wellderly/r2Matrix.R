## This function takes a 2x2 matrix and finds the r squared value
## r^2 = (D^2)/(pA1*pA2*pB1*PB2)
## D = pA1B1 - pA1*pB1

r2value <- function(matrixf){
  if (is.matrix(matrixf)==FALSE){
    #print.default("r2 value function was not given a matrix as input")
    return(NA)
  } else if(nrow(matrixf) !=2 || ncol(matrixf) != 2){
    #print.default("r2 value function was not given a 2x2 matrix as input")
    return(NA)
  } else {
    sumN <- sum(matrixf)
    freqRows <- apply(matrixf,1,sum)/sumN
    freqColumns <-apply(matrixf,2,sum)/sumN 
    freqMatrix <- matrixf/sumN
    return(((freqMatrix[1,1]-freqRows[1]*freqColumns[1])^2)/(freqRows[1]*freqRows[2]*freqColumns[1]*freqColumns[2]))
  }
}
#print.default(r2value(matrixf))

#matrixf <- matrix(c(10,0,0,10))
