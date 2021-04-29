outputfolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wellderly/"
outputfile <- read.table(paste0(outputfolder,"zalphaoutput4425.txt"))


shrinkfile <- outputfile[outputfile[,1]>=135100000 & outputfile[,1]<=137900000,]
                   
plot(shrinkfile,main=expression(paste("Z",''[alpha]," applied to Wellderly data")),ylab=expression(paste("Z",''[alpha])),xlab="bp on Chr2",type="l",xaxs="i",yaxs="i",xlim=c(135000000,138000000),ylim=c(0.05,0.35),las=1)

abline(v=136608646,col="red")
legend("topright","-13910C>T",col="red",lty=1)
