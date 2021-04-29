
require(qqman)

chrMap<-list()

for (i in 1:38){
  chrMap[[i]]<-read.table(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/LDhat/Maps/cM_map_chr",i,".txt"),stringsAsFactors = FALSE,header=TRUE)
  
}


SNP<-NULL
CHR<-NULL
BP<-NULL
cMMb<-NULL
for (i in 1:38){
  SNP<-c(SNP,paste0(i,"_",chrMap[[i]]$loci))
  CHR<-c(CHR,rep(i,nrow(chrMap[[i]])))
  BP<-c(BP,chrMap[[i]]$loci)
  cMMb<-c(cMMb,chrMap[[i]]$cMMb)
}

manhatdf<-data.frame(SNP=SNP,CHR=CHR,BP=BP,cMMb=cMMb)
atloc<-cumsum(as.numeric(tapply(manhatdf$BP,manhatdf$CHR,max)))-tapply(manhatdf$BP,manhatdf$CHR,quantile,0.5)

manhattan(manhatdf, p = "cMMb", logp = FALSE, ylab = "cM per Mb", genomewideline = FALSE,
          suggestiveline = FALSE, xaxt="n",type="l")#,highlight=highlight)
axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)


require(intervals)
# ----------------------------------------------------------------------------------------
getRate <- function(map, intI) {
  # get recombination rate in cM/Mb for a given map and interval
  rr1 <- approx( map$loci, map$cumulative_cM, xout=intI[,1] )
  rr2 <- approx( map$loci, map$cumulative_cM, xout=intI[,2] )
  rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # cM/Mb
  return(rate)
}
# ----------------------------------------------------------------------------------------
binGenome_per_chr <- function(bsize, df) {
  # zero-based indexing, half open intervals
  bpr <- c(1,df$loci[nrow(df)])
  bpr[1] <- bpr[1] - bpr[1]%%bsize
  bpr[2] <- bpr[2] + (bsize - bpr[2]%%bsize)
  x <- seq(bpr[1],bpr[2],by=bsize)
  seqbinI <- Intervals( cbind(x[-length(x)],x[-1]), closed=c(TRUE,FALSE), type="Z" )

  return(seqbinI)
}
#----------------------

chrMapSmooth<-list()
for (i in 1:38){
  seqint<-binGenome_per_chr(1000000,chrMap[[i]])
  intI<-seqint@.Data
  rate<-getRate(chrMap[[i]],intI)
  chrMapSmooth[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
  chrMapSmooth[[i]]<-chrMapSmooth[[i]][is.na(rate)==FALSE,]
}

## Smoothed Manhatten plot

SNP<-NULL
CHR<-NULL
BP<-NULL
cMMb<-NULL
for (i in 1:38){
  SNP<-c(SNP,paste0(i,"_",chrMapSmooth[[i]]$int_centre))
  CHR<-c(CHR,rep(i,nrow(chrMapSmooth[[i]])))
  BP<-c(BP,chrMapSmooth[[i]]$int_centre)
  cMMb<-c(cMMb,chrMapSmooth[[i]]$rate)
}

manhatdfsmooth<-data.frame(SNP=SNP,CHR=CHR,BP=BP,cMMb=cMMb)
atlocsmooth<-cumsum(as.numeric(tapply(manhatdfsmooth$BP,manhatdfsmooth$CHR,max)))-tapply(manhatdfsmooth$BP,manhatdfsmooth$CHR,quantile,0.5)

manhattan(manhatdfsmooth, p = "cMMb", logp = FALSE, ylab = "cM per Mb", genomewideline = FALSE,
          suggestiveline = FALSE, xaxt="n",type="l")#,highlight=highlight)
axis(1,at=atlocsmooth[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
axis(1,at=atlocsmooth[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)



## combined
jpeg("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/LDhat/RecombinationMap.jpeg",res=300,height=15,width=25,units='cm',quality=100)
par(mfrow=c(2,1),mar=c(1,6,4,2))
manhattan(manhatdf, p = "cMMb", logp = FALSE, ylab = "cM per Mb", genomewideline = FALSE,
          suggestiveline = FALSE,xlab="", xaxt="n",type="l",ylim=c(0,40))#,highlight=highlight)
#axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
#axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
par(mar=c(5,6,0,2))
manhattan(manhatdfsmooth, p = "cMMb", logp = FALSE, ylab = "cM per Mb\n1 Mb intervals", genomewideline = FALSE,
          suggestiveline = FALSE, xaxt="n",type="l")#,highlight=highlight)
axis(1,at=atlocsmooth[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
axis(1,at=atlocsmooth[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
dev.off()

