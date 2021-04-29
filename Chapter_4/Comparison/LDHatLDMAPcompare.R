
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22/LDMAP comparison")

bmp("LDhatLDMAPBherer.bmp",res=300,height=14,width=15,units='cm')
par(mfrow=c(1,1))
par(mar=c(5,4,4,6)+0.1)

### LDhat
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22")
startPos <- 20000.428
endPos <- 51219.641
ldhat <- read.table("res.txt",header=TRUE,colClasses = "numeric")
ldhat<-ldhat[-1,]
ldhat$kbpos<-ldhat$Loci+startPos-0.001


##  CM
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22")
source("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/convertLDhatToCMperMB2.R")
trimmedStart<- 20000428
trimmedEnd<- 51218377
map_length_cM<-60.67105
ldhatcM<-convertLDhatTocMperMb2(ldhat$Loci,ldhat$Mean_rho,startPos,endPos,trimmedStart,trimmedEnd,bp=0.001,kb=TRUE,map_length_cM)
plot(ldhatcM$loci,ldhatcM$cumulative_cM,type="l",main="Comparison of Maps",xlim=c(20000000,52000000),ylim=c(0,61),xlab="",ylab="",axes=FALSE,lty=3)
axis(2,ylim=c(0,70))
mtext("cM",side=2,line=2.5)
box()


### LDMAP
par(new=TRUE)
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22/LDMAP comparison")
ldmap <- read.csv("concat_CEU_csv_list.txt.ter.csv",header=TRUE,colClasses = c("numeric","character","numeric","numeric","character"))
plot(ldmap$KBmap,ldmap$LDUmap,type="l",col="red",axes=FALSE,xlab="",ylab="",ylim=c(0,1022),xlim=c(20000.000,52000.000),lty=1)
mtext("LDU",side=4,line=4)
axis(4,ylim=c(0,1100))

axis(1,xlim=c(20000.000,52000.000))
mtext("chromosome 22 Kb",side=1,line=2.5)




###bherer
par(new=TRUE)
sexavgChr22 <- read.delim("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wavelet/Bherer/Refined_EUR_genetic_map_b37/Refined_EUR_genetic_map_b37/sexavg_chr22.txt")
## trim top and recalculate cM
sexavgChr22<-sexavgChr22[910:nrow(sexavgChr22),]
sexavgChr22$cM[1]<-0
for (i in 2:nrow(sexavgChr22)){
  sexavgChr22$cM[i]<-sexavgChr22$rate[i-1]*(sexavgChr22$pos[i]-sexavgChr22$pos[i-1])/1000000+sexavgChr22$cM[i-1]
}
plot(sexavgChr22$pos/1000,sexavgChr22$cM,type="l",col="blue",axes=FALSE,xlab="",ylab="",ylim=c(0,61),xlim=c(20000.000,52000.000),lty=2)

legend("topleft",legend=c("LDhat","LDMAP","Bherer"),col=c("black","red","blue"),lty=c(3,1,2))
dev.off()


#######################

# Add loci and lDUperMb column to ldmap

ldmap$loci<-ldmap$KBmap*1000
ldmap$ldupermb <- 0
for (i in 1:(nrow(ldmap)-1)){
  ldmap$ldupermb[i] <-1000000*(ldmap$LDUmap[i+1]-ldmap$LDUmap[i])/(ldmap$loci[i+1]-ldmap$loci[i])
}

plot(ldmap$KBmap,ldmap$ldupermb)

#############################
#required functions 
# ----------------------------------------------------------------------------------------
getRate <- function(mappos, maprate, intI) { #loci, cumulative rate, interval file
  # get recombination rate in rate/Mb for a given map and interval
  rr1 <- approx( mappos, maprate, xout=intI[,1] )
  rr2 <- approx( mappos, maprate, xout=intI[,2] )
  rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # rate/Mb
  return(rate)
}


# ----------------------------------------------------------------------------------------
binGenome_per_chr <- function(chr, bsize, rm.gaps=TRUE,genomeData,gap=NULL) {
  # zero-based indexing, half open intervals
  bpr <- c(1,genomeData$nbases[ genomeData$chrom==chr ] )
  bpr[1] <- bpr[1] - bpr[1]%%bsize
  bpr[2] <- bpr[2] + (bsize - bpr[2]%%bsize)
  x <- seq(bpr[1],bpr[2],by=bsize)
  seqbinI <- Intervals( cbind(x[-length(x)],x[-1]), closed=c(TRUE,FALSE), type="Z" )
  if(rm.gaps) 
  {
    gapI <- with( gap[gap$chrom==chr,], {
      Intervals(cbind(chromStart, chromEnd), closed=c(TRUE,FALSE), type="Z" )})
    ol <- unlist(interval_overlap(from=gapI, to=seqbinI))
    if(length(ol)>0) seqbinI <- seqbinI[-ol, ]
  }
  return(seqbinI)
}
# ----------------------------------------------------------------------------------------

#############################
#
# put all three datsets in one intervalled data frame

require(intervals)
genomeData<-data.frame(chrom="chr22",nbases=51218000)
seqint<-binGenome_per_chr("chr22",1000,rm.gaps = FALSE,genomeData)
intI<-seqint@.Data
mergedData<-data.frame(
  loci=intI[,1],
  ldhat=getRate(ldhatcM$loci,ldhatcM$cumulative_cM,intI),
  ldmap=getRate(ldmap$loci,ldmap$LDUmap,intI),
  bherer=getRate(sexavgChr22$pos,sexavgChr22$cM,intI)
)
mergedData<-na.omit(mergedData)
plot(mergedData$loci,mergedData$ldhat,type="l")
lines(mergedData$loci,mergedData$ldmap,col="blue")
lines(mergedData$loci,mergedData$bherer,col="red")

cor(mergedData$ldhat,mergedData$ldmap)
cor(mergedData$ldhat,mergedData$bherer)
cor(mergedData$ldmap,mergedData$bherer)

summary(lm(mergedData$bherer~mergedData$ldhat))
anova(lm(mergedData$bherer~mergedData$ldhat))

summary(lm(mergedData$bherer~mergedData$ldmap))
anova(lm(mergedData$bherer~mergedData$ldmap))

summary(lm(mergedData$ldmap~mergedData$ldhat))
anova(lm(mergedData$ldmap~mergedData$ldhat))

bob<-hist(mergedData$bherer)
