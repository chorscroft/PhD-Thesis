##removes random SNPs until average spacing between 
##polymorphic sites is 2500bp
##returns a data frame


removeRandomSNPs <- function(x,bp=2500){
  if (is.data.frame(x)==FALSE){
    #print.default("removeRandomSNPs function was not given a data frame as input")
    return(NA)
  }
  if (is.numeric(bp)==FALSE){
    #print.default("removeRandomSNPs function was not given a number for bp as input")
    return(NA)    
  }
  while (averageDistSNPs(x) < bp){
    #remove random SNP
    removeSNP <- sample(1:nrow(x),1)
    x <- x[-c(removeSNP),]
  }
  return(x)
}
averageDistSNPs <- function(x){
  sumDist <- 0
  for (i in 1:(nrow(x)-1)){
    sumDist <- sumDist + x$bp[i+1]-x$bp[i]
  }
  return(sumDist/(nrow(x)-1))
}
