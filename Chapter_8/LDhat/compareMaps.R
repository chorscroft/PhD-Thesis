
getPlot<-function(map,start,end){
  mapPlot<-data.frame(pos=map$pos,cM=map$cM)
  # look for first bp in map that is greater than or equal to the start value
  for (x in 1:nrow(mapPlot)){
    if (mapPlot$pos[x]>=start){
      break
    }
  }
  
  if (mapPlot$pos[x]==start){ #if it is equal then trim off the previous loci
    mapPlot<-mapPlot[c(x:nrow(mapPlot)),]
  } else { # if not then adjust linerally between this snp and the previous
    cMstart<-mapPlot$cM[x-1]+(mapPlot$cM[x]-mapPlot$cM[x-1])*(start-mapPlot$pos[x-1])/(mapPlot$pos[x]-mapPlot$pos[x-1])
    mapPlot<-rbind(c(start,cMstart),mapPlot[c(x:nrow(mapPlot)),])
  }
  
  # look for first bp in map that is less than or equal to the end value
  for (x in nrow(mapPlot):1){
    if (mapPlot$pos[x]<=end){
      break
    }
  }
  
  if (map$pos[x]==end){ #if it is equal then trim off the end of the map
    mapPlot<-mapPlot[c(1:x),]
  } else { # if not then adjust linerally between this snp and the previous
    cMend<-mapPlot$cM[x]+(mapPlot$cM[x+1]-mapPlot$cM[x])*(end-mapPlot$pos[x])/(mapPlot$pos[x+1]-mapPlot$pos[x])
    mapPlot<-rbind(mapPlot[c(1:x),],c(end,cMend))
  }
  
  #adjust so it starts at 0
  mapPlot$cM<-mapPlot$cM-mapPlot$cM[1]
  return(mapPlot)
}

#######################################################################

## get counts for SNP density comparison
snpDensity<-data.frame(chr=c(1:38),bp=0,ldhat_snps=0,auton_snps=0,campbell_snps=0)

ldhatMapPlot<-list()
autonMapPlot<-list()
campbellMapPlot<-list()

## Loop for each chromosome
for (i in 1:38){
  
  ## get LDhat Map
  ldhatMap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",i,"/cM_map_chr",i,".txt"),header=TRUE,stringsAsFactors=FALSE)
  colnames(ldhatMap)<-c("pos","cM","rate")
  
  ## get auton Map
  autonMap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/AutonMaps/mark4_cleaned_chr",i,".cf3.1.sorted.txt"),header=FALSE,skip=1,stringsAsFactors=FALSE)
  colnames(autonMap)<-c("chr","pos","rate","cM")
  ## get Campbell Map
  campbellMap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/CampbellMaps/chr",i,"_average_canFam3.1.txt"),header=TRUE,stringsAsFactors=FALSE)
  
  ## get start and end from maps, so that each map fits
  start<-max(ldhatMap$pos[1],autonMap$pos[1],campbellMap$pos[1])
  end<-min(ldhatMap$pos[nrow(ldhatMap)],autonMap$pos[nrow(autonMap)],campbellMap$pos[nrow(campbellMap)])
  
  
  # create plots of ldhat vs auton vs campbell
  options(scipen=999)
  ldhatMapPlot[[i]]<-getPlot(ldhatMap,start,end)
  autonMapPlot[[i]]<-getPlot(autonMap,start,end)
  campbellMapPlot[[i]]<-getPlot(campbellMap,start,end)
  
  jpeg(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/compareMaps/cumulativecM/compare_ccM_chr",i,".jpeg"),res=300,height=14,width=14,units='cm',quality=100)
  plot(ldhatMapPlot[[i]]$pos/1000000,ldhatMapPlot[[i]]$cM,type="l",xlab="Position (Mb)",ylab="Cumulative cM",main=paste("Chromosome",i))
  lines(autonMapPlot[[i]]$pos/1000000,autonMapPlot[[i]]$cM,col="blue",lty=2)
  lines(campbellMapPlot[[i]]$pos/1000000,campbellMapPlot[[i]]$cM,col="red",lty=3)
  legend("topleft",legend=c("LDhat","Auton","Campbell"),lty=c(1,2,3),col=c("black","blue","red"))
  dev.off()
  
  
  # get SNP density information
  snpDensity$bp[i]<-end-start+1
  snpDensity$ldhat_snps[i]<-nrow(ldhatMapPlot[[i]])
  snpDensity$auton_snps[i]<-nrow(autonMapPlot[[i]])
  snpDensity$campbell_snps[i]<-nrow(campbellMapPlot[[i]])
}

#snp denity:
#ldhat:
1000000*sum(snpDensity$ldhat_snps)/sum(snpDensity$bp)
#auton:
1000000*sum(snpDensity$auton_snps)/sum(snpDensity$bp)
#campbell:
1000000*sum(snpDensity$campbell_snps)/sum(snpDensity$bp)

## interpolate into bins
require(intervals)
# ----------------------------------------------------------------------------------------
getRate <- function(map, intI) {
  # get recombination rate in cM/Mb for a given map and interval
  rr1 <- approx( map$pos, map$cM, xout=intI[,1] )
  rr2 <- approx( map$pos, map$cM, xout=intI[,2] )
  rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # cM/Mb
  return(rate)
}
# ----------------------------------------------------------------------------------------
binGenome_per_chr <- function(bsize, df) {
  # zero-based indexing, half open intervals
  bpr <- c(1,df$pos[nrow(df)])
  bpr[1] <- bpr[1] - bpr[1]%%bsize
  bpr[2] <- bpr[2] + (bsize - bpr[2]%%bsize)
  x <- seq(bpr[1],bpr[2],by=bsize)
  seqbinI <- Intervals( cbind(x[-length(x)],x[-1]), closed=c(TRUE,FALSE), type="Z" )
  
  return(seqbinI)
}
#----------------------


## WAVELET ANALYSIS
require(wmtsa)

ldhat_interp_w<-list()
auton_interp_w<-list()
campbell_interp_w<-list()
for (i in 1:38){
  seqint<-binGenome_per_chr(500000,ldhatMapPlot[[i]])
  intI<-seqint@.Data
  rate<-log10(getRate(ldhatMapPlot[[i]],intI))
  ldhat_interp_w[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
  ldhat_interp_w[[i]]<-ldhat_interp_w[[i]][is.na(rate)==FALSE,]
  
  seqint<-binGenome_per_chr(500000,autonMapPlot[[i]])
  intI<-seqint@.Data
  rate<-log10(getRate(autonMapPlot[[i]],intI))
  auton_interp_w[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
  auton_interp_w[[i]]<-auton_interp_w[[i]][is.na(rate)==FALSE,]
  
  seqint<-binGenome_per_chr(500000,campbellMapPlot[[i]])
  intI<-seqint@.Data
  rate<-log10(getRate(campbellMapPlot[[i]],intI))
  campbell_interp_w[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
  campbell_interp_w[[i]]<-campbell_interp_w[[i]][is.na(rate)==FALSE,]
}

getPowerSpectrum<-function(dataDWT){
  powerSpectrum<-NULL
  for(i in 1:dataDWT$dictionary$n.levels){
    powerSpectrum[i]<-sum(dataDWT$data[[i]]^2,na.rm=TRUE)
  }
  powerSpectrum<-powerSpectrum/sum(powerSpectrum,na.rm=TRUE)
  return(powerSpectrum)
}


ldhat_modwt<-list()
auton_modwt<-list()
campbell_modwt<-list()
for (i in 1:38){
  ldhat_modwt[[i]]<-wavMODWT(ldhat_interp_w[[i]]$rate,wavelet="haar",n.levels=5,keep.series=TRUE)
  auton_modwt[[i]]<-wavMODWT(auton_interp_w[[i]]$rate,wavelet="haar",n.levels=5,keep.series=TRUE)
  campbell_modwt[[i]]<-wavMODWT(campbell_interp_w[[i]]$rate,wavelet="haar",n.levels=5,keep.series=TRUE)
}

## Power Spectrum
ldhat_ps<-sapply(ldhat_modwt,getPowerSpectrum)
auton_ps<-sapply(auton_modwt,getPowerSpectrum)
campbell_ps<-sapply(campbell_modwt,getPowerSpectrum)

jpeg("/temp/hgig/EXOME_DATA/Clare/Dogs/compareMaps/PowerSpectrums.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(apply(ldhat_ps,1,mean,na.rm=TRUE),type="l",xlab="Scale Mb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(ldhat_ps,auton_ps,campbell_ps,na.rm=TRUE)))
axis(1,1:5,labels=2^(1:5)*0.5,las=2)
lines(apply(auton_ps,1,mean,na.rm=TRUE),col="blue",lty=2)
lines(apply(campbell_ps,1,mean,na.rm=TRUE),col="red",lty=3)
legend("topright",legend=c("LDhat","Auton","Campbell"),col=c("black","blue","red"),lty=c(1,2,3))
dev.off()

## Correlations
## Wavelet correlation
#function
getWaveletCorrelation<-function(x,y){ ##x and y are lists
  est<-NULL
  pval<-NULL
  for (i in 1:5){  ##need more than 2 observations
    rankCor<-cor.test(x[[i]],y[[i]],method="kendall")
    est[i]<-rankCor$estimate
    pval[i]<-rankCor$p.value
  }
  return(list(est=est,pval=pval))
}
formatData<-function(modwt){
  scaleData<-list()
  ## get each of the decompositions
  for (s in 1:modwt[[1]]$dictionary$n.levels){
    scaleData[[s]]<-modwt[[1]]$data[[s]]
    for (i in 2:38){
      scaleData[[s]]<-c(scaleData[[s]],modwt[[i]]$data[[s]]) 
    }
  }
  return(scaleData)
}
#get correlations
ldhat_auton_cor<-getWaveletCorrelation(formatData(ldhat_modwt),formatData(auton_modwt))
ldhat_campbell_cor<-getWaveletCorrelation(formatData(ldhat_modwt),formatData(campbell_modwt))
auton_campbell_cor<-getWaveletCorrelation(formatData(auton_modwt),formatData(campbell_modwt))

jpeg("/temp/hgig/EXOME_DATA/Clare/Dogs/compareMaps/Correlations.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(ldhat_auton_cor$est[1:5],type="l",ylim=c(0,1),col="black",xaxt="n",xlab="Scale Mb",ylab="Kendall's tau")
axis(1,at=c(1:5),labels=c(2^(1:5)*0.5),las=2)
points(which(ldhat_auton_cor$pval[1:5]<0.01),ldhat_auton_cor$est[which(ldhat_auton_cor$pval[1:5]<0.01)],col="black")
lines(ldhat_campbell_cor$est[1:5],col="blue",lty=2)
points(which(ldhat_campbell_cor$pval[1:5]<0.01),ldhat_campbell_cor$est[which(ldhat_campbell_cor$pval[1:5]<0.01)],col="blue")
lines(auton_campbell_cor$est[1:5],col="red",lty=3)
points(which(auton_campbell_cor$pval[1:5]<0.01),auton_campbell_cor$est[which(auton_campbell_cor$pval[1:5]<0.01)],col="red")
legend("topleft",legend=c("LDhat and Auton","LDhat and Campbell","Auton and Campbell"),col=c("black","blue","red"),lty=c(1,2,3))
legend("bottomright",legend="p-value < 0.01",pch=1)
dev.off()

## Correlations between original data
## (Remember wavelet correlations look at *change* whereas this is straight)


ldhat_orig<-ldhat_modwt[[1]]$series@data
auton_orig<-auton_modwt[[1]]$series@data
campbell_orig<-campbell_modwt[[1]]$series@data
for (i in 2:38){
  ldhat_orig<-c(ldhat_orig,ldhat_modwt[[i]]$series@data) 
  auton_orig<-c(auton_orig,auton_modwt[[i]]$series@data) 
  campbell_orig<-c(campbell_orig,campbell_modwt[[i]]$series@data) 
}
laCor<-cor.test(ldhat_orig,auton_orig,method="kendall")
lcCor<-cor.test(ldhat_orig,campbell_orig,method="kendall")
acCor<-cor.test(auton_orig,campbell_orig,method="kendall")
#ldhat auton
laCor
#ldhat campbell
lcCor
#auton campbell
acCor

sd(ldhat_orig)
sd(auton_orig,na.rm=TRUE)
sd(campbell_orig)

##### non wavelet code
# scaleCor<-data.frame(scale=c(500000,1000000,2000000,4000000,8000000),ldhat_auton=NA,ldhat_campbell=NA,auton_campbell=NA)
# for (s in 1:nrow(scaleCor)){
#   ldhat_interp<-list()
#   auton_interp<-list()
#   campbell_interp<-list()
#   for (i in 1:38){
#     seqint<-binGenome_per_chr(scaleCor$scale[s],ldhatMapPlot[[i]])
#     intI<-seqint@.Data
#     rate<-getRate(ldhatMapPlot[[i]],intI)
#     ldhat_interp[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
#     ldhat_interp[[i]]<-ldhat_interp[[i]][is.na(rate)==FALSE,]
#     
#     seqint<-binGenome_per_chr(scaleCor$scale[s],autonMapPlot[[i]])
#     intI<-seqint@.Data
#     rate<-getRate(autonMapPlot[[i]],intI)
#     auton_interp[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
#     auton_interp[[i]]<-auton_interp[[i]][is.na(rate)==FALSE,]
#     
#     seqint<-binGenome_per_chr(scaleCor$scale[s],campbellMapPlot[[i]])
#     intI<-seqint@.Data
#     rate<-getRate(campbellMapPlot[[i]],intI)
#     campbell_interp[[i]]<-data.frame(int1=intI[,1],int2=intI[,2],int_centre=(intI[,1]+intI[,2])/2,rate=rate)
#     campbell_interp[[i]]<-campbell_interp[[i]][is.na(rate)==FALSE,]
#   }
#   
#   
#   ldhat_concat<-NULL
#   auton_concat<-NULL
#   campbell_concat<-NULL
#   for (i in 1:38){
#     ldhat_concat<-c(ldhat_concat,ldhat_interp[[i]]$rate)
#     auton_concat<-c(auton_concat,auton_interp[[i]]$rate)
#     campbell_concat<-c(campbell_concat,campbell_interp[[i]]$rate)
#   }
#   scaleCor[s,2]<-cor(ldhat_concat,auton_concat)^2
#   scaleCor[s,3]<-cor(ldhat_concat,campbell_concat)^2
#   scaleCor[s,4]<-cor(auton_concat,campbell_concat)^2
# }
# scaleCor$scale
# plot(1:5,scaleCor$ldhat_auton,type="l",ylim=c(0,1))
# lines(1:5,scaleCor$ldhat_campbell,col="blue")
# lines(1:5,scaleCor$auton_campbell,col="red")

