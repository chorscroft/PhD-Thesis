## Gets the ldhat output and transforms it into 
## input for wavelet analysis

getWaveletFile<-function(res,startPos,endPos,trimmedStart,trimmedEnd,map_length_cM){
  res<-res[-1,]
  res$pos<-res$Loci+startPos-0.001
  rescM<-convertLDhatTocMperMb2(res$Loci,res$Mean_rho,startPos,endPos,trimmedStart,trimmedEnd,bp=0.001,kb=TRUE,map_length_cM)
  colnames(rescM)<-c("pos","cM","rate")
  
  genomeData<-data.frame(chrom="chr22",nbases=rescM$pos[nrow(rescM)])
  
  gap<-data.frame(chrom=character(),chromStart=double(),chromEnd=double(),stringsAsFactors = FALSE)
  gaps <- 0
  for (i in 1:(nrow(rescM)-1)){
    if (rescM$pos[i+1]-rescM$pos[i] >= 50000){
      gaps <- gaps + 1
      gap[gaps,]<-list("chr22",as.integer(rescM$pos[i]),as.integer(rescM$pos[i+1]))
    }
  }
  
  seqint<-binGenome_per_chr("chr22",1000,rm.gaps = TRUE,genomeData,gap)
  intI<-seqint@.Data
  
  rate<-getRate(rescM,intI)
  lograte<-log10(rate)
  binRates<-cbind(intI,rate,lograte)
  return(data.frame(binRates))
}


# ----------------------------------------------------------------------------------------
getRate <- function(map, intI) {
  # get recombination rate in cM/Mb for a given map and interval
  rr1 <- approx( map$pos, map$cM, xout=intI[,1] )
  rr2 <- approx( map$pos, map$cM, xout=intI[,2] )
  rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # cM/Mb
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
getPowerSpectrum<-function(dataDWT){
  powerSpectrum<-NULL
  for(i in 1:dataDWT$dictionary$n.levels){
    powerSpectrum[i]<-sum(dataDWT$data[[i]]^2,na.rm = TRUE)
    if (powerSpectrum[i]==0){
      powerSpectrum[i]<-NA
    }
  }
  powerSpectrum<-powerSpectrum/sum(powerSpectrum,na.rm=TRUE)
  return(powerSpectrum)
}
# ----------------------------------------------------------------------------------------
getDifference<-function(mergedFile,initial1,initial2){
  
  diff<-cbind(mergedFile[,c(1,2)],logdiff=mergedFile[[paste0(initial1,"_lograte")]]-mergedFile[[paste0(initial2,"_lograte")]])
  
  #trim NAs from top and bottom of file
  
  while (is.na(diff$logdiff[1])){
    diff<-diff[-1,]
  }
  while (is.na(diff$logdiff[nrow(diff)])){
    diff<-diff[-nrow(diff),]
  }
  
  row.names(diff)<-c(1:nrow(diff))
  
  #Pad with NAs to get a power of 2
  pwr<-1
  while(pwr<nrow(diff)){
    pwr<-pwr*2
  }
  if (nrow(diff)<pwr){
    for (i in (nrow(diff)+1):pwr){
      diff[i,]<-c(diff[i-1,1]+1000,diff[i-1,2]+1000,NA)
    }
  }
  return(diff)
}
# ----------------------------------------------------------------------------------------
getVariance<-function(dataDWT){
  calcvar<-0
  count<-0
  for (i in 1:dataDWT$dictionary$n.levels){
    calcvar<-calcvar+sum(dataDWT$data[[i]]^2,na.rm=TRUE)
    count<-count+sum(!is.na(dataDWT$data[[i]]))
  }
  return(calcvar/count)
}
# ----------------------------------------------------------------------------------------

getWaveletFormat<-function(data){
  #trim NAs from top and bottom of file
  
  while (is.na(data[1,3])){
    data<-data[-1,]
  }
  while (is.na(data[nrow(data),3])){
    data<-data[-nrow(data),]
  }
  
  row.names(data)<-c(1:nrow(data))
  
  #Pad with NAs to get a power of 2
  pwr<-1
  while(pwr<nrow(data)){
    pwr<-pwr*2
  }
  if (nrow(data)<pwr){
    for (i in (nrow(data)+1):pwr){
      data[i,]<-c(data[i-1,1]+1000,data[i-1,2]+1000,NA)
    }
  }
  return(data)
}
# ----------------------------------------------------------------------------------------
getWaveletCorrelation<-function(x,y){ ##x and y are DWT objects
  est<-NULL
  pval<-NULL
  for (i in 1:13){  ##need more than 2 observations
    rankCor<-cor.test(x$data[[i]],y$data[[i]],method="kendall")
    est[i]<-rankCor$estimate
    pval[i]<-rankCor$p.value
  }
  return(list(est=est,pval=pval))
}
# ----------------------------------------------------------------------------------------

plotLegendCWT<-function(x){
  x$power <- x$power.corr
  yvals <- log2(x$period)
  zvals <- log2(abs(x$power/x$sigma2))
  zlim <- range(c(-1, 1) * max(zvals))
  zvals[zvals < zlim[1]] <- zlim[1]
  locs <- pretty(range(zlim), n = 5)
  leg.lab<-2^locs
  leg.lab[leg.lab<1] <- paste0("1/",(1/leg.lab[leg.lab<1]))
  fill.cols <- c("#00007F", "blue", "#007FFF", "cyan", 
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  col.pal <- colorRampPalette(fill.cols)
  fill.colors <- col.pal(64)
  image.plot(x$t, yvals, t(zvals), 
             zlim = zlim, 
             ylim = rev(range(yvals)), 
             col = fill.colors, 
             horizontal = FALSE, 
             legend.only = TRUE, 
             axis.args = list(at = locs,labels = leg.lab),
             xpd = NA)
}

# ----------------------------------------------------------------------------------------
getCumulativeCM<-function(res,startPos,endPos,trimmedStart,trimmedEnd,map_length_cM){
  res<-res[-1,]
  res$pos<-res$Loci+startPos-0.001
  rescM<-convertLDhatTocMperMb2(res$Loci,res$Mean_rho,startPos,endPos,trimmedStart,trimmedEnd,bp=0.001,kb=TRUE,map_length_cM)
  colnames(rescM)<-c("pos","cM","rate")
  return(rescM)
  
}
