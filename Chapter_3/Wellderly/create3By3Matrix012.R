#Creates the 3x3 matrix from a 012 tped file data frame given the lines
#x is a tped file data frame
#i is the first SNP
#j is the second SNP

create3By3Matrix <- function(x,i,j) {
  if (is.data.frame(x)==FALSE){
    #print.default("create3By3Matrix function was not given a data frame as input")
    return(NA)
  } else if (is.numeric(i)==FALSE || is.numeric(j)==FALSE){
    #print.default("create3By3Matrix function was not given two numbers for SNP rows as input")
    return(NA)
  } else {
    outMatrix <- matrix(0,nrow = 3, ncol = 3)
    for (xcol in seq(5,ncol(x),by=2)) {
      if (x[i,xcol]==0 || x[i,xcol+1]==0 || x[j,xcol]==0 || x[j,xcol+1]==0){
        #skip
      } else {
        genoA = 0 #1 = AA, 2 = Aa, 3 = aa
        genoB = 0 #1 = BB, 2 = Bb, 3 = bb
        
        #Define Genotype of A
        if (x[i,xcol] == 1 && x[i,xcol+1] == 1) {
          genoA = 1
        } else if (x[i,xcol] == 2 && x[i,xcol+1] == 2) {
          genoA = 3
        } else {
          genoA = 2
        }
        
        #Define Genotype of B
        if (x[j,xcol] == 1 && x[j,xcol+1] == 1) {
          genoB = 1
        } else if (x[j,xcol] == 2 && x[j,xcol+1] == 2) {
          genoB = 3
        } else {
          genoB = 2
        }
        
        #Update Matrix
        outMatrix[genoA,genoB] <- outMatrix[genoA,genoB] +1
      }
    }
    return(outMatrix)
  }
}
