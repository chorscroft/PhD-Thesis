
require(fields)
simFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/RecombinationAcross/Sims/"
midsection <- c(200001,300000)

avrr2mat<-matrix(NA,nrow=100,ncol=100)
countr2mat<-matrix(0,nrow=100,ncol=100)
sumr2mat<-matrix(0,nrow=100,ncol=100)
removeMirror <- function(x){
  for (i in 1:nrow(x)){
    x[i,i:nrow(x)] <- NA
  }
  return(x)
}
round1000 <- function(x){
  return(floor(x/1000))
}

rotate  <- function(x, clockwise=T) {
  if (clockwise) { t( apply(x, 2, rev))
  } else {apply( t(x),2, rev)} 
}


for (sim in 1:100){
  ## Read in bp locations
  simBPs <- unlist(read.table(paste0(simFolder,"simBPs_",sim,".txt")),use.names = FALSE)  
  ## Read in matrix of snp status:
  ## rows are chromosomes
  ## columns are each snp
  ## 0 = ancenstral; 1 = derived
  simMatrix <-data.matrix(read.table(paste0(simFolder,"simMatrix_",sim,".txt")))
  colnames(simMatrix)<-NULL
  
  nSNPs <- length(simBPs)
  r2matrix <-cor(simMatrix)^2
  
  
  smallr2matrix <-r2matrix[simBPs>=midsection[1] & simBPs<=midsection[2],simBPs>=midsection[1] & simBPs<=midsection[2]]
  smallSimBPs <-simBPs[simBPs>=midsection[1] & simBPs<=midsection[2]]
  
  smallr2matrixmir <- removeMirror(smallr2matrix)
  for(i in 1:(nrow(smallr2matrixmir)-1)){
    ii <- round1000(smallSimBPs[i]-midsection[1])+1
    for(j in (i+1):nrow(smallr2matrixmir)){
      jj <- round1000(smallSimBPs[j]-midsection[1])+1
      countr2mat[jj,ii] <- countr2mat[jj,ii]+1
      sumr2mat[jj,ii]<-sumr2mat[jj,ii]+smallr2matrixmir[j,i]
    }
  }
}
for (i in 1:100){
  for(j in 1:100){
    if (countr2mat[i,j] > 0){
      avrr2mat[i,j] <- sumr2mat[i,j]/countr2mat[i,j]
    }
  }
}
#write(t(avrr2mat),file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/RecombinationAcross/avrr2mat.txt",ncolumns = 201)

bmp("//filestore.soton.ac.uk/users/ch19g17/mydocuments/RecombinationAcross/r2across.bmp",units = "cm",width = 17,height = 14,res=300)
par(mar=c(5,4.5,4,7))
image(rotate(avrr2mat), axes=F,col=tim.colors(),xlab="Kb",ylab="Kb")
mtext(text=seq(300,200,-10), side=2, line=0.3, at=seq(0,1,0.1), las=1, cex=0.8)
mtext(text=seq(200,300,10), side=1, line=0.3, at=seq(0,1,0.1), las=2, cex=0.8)
#axis(1,at=c(0,0.5,1),labels=c(2,2.5,3))
#axis(2,at=c(0,0.5,1),labels=c(2,2.5,3))
image.plot(avrr2mat, legend.only=T)
title(main=expression(paste(italic('r'^2)," values 200 generations after fixation")))
dev.off()



