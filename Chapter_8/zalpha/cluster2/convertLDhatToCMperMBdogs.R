

convertLDhatTocMperMb<-function(loci,rate,start,end,map_length_cM){
  
  loci<-loci*1000+start-1   #convert loci from map back to bp
  loci<-c(loci,end)         #add last value back on
  
  #work out cumulative rho
  cumulative_rho<-0
  for (i in 2:length(loci)){
    cumulative_rho[i]<-rate[i-1]*(loci[i]-loci[i-1])/1000+cumulative_rho[i-1]   #rates are in rho/Kb so adjust by 1000 to get rho/Mb
  }
  #total cumulative rho cumulative_rho[length(cumulative_rho)] should be the same as in the first line of res.txt
  conversion_factor<-map_length_cM/cumulative_rho[length(cumulative_rho)]
  cumulative_cM<-cumulative_rho*conversion_factor
  
  cMMb<-rep(0,length(loci))
  for (i in 1:(length(loci)-1)){
    cMMb[i]<-(cumulative_cM[i+1]-cumulative_cM[i])/((loci[i+1]-loci[i])/1000000)
  }
  cmReturn<-data.frame(loci,cumulative_cM,cMMb)
  
  return(cmReturn)
  
  
 
}
getMapLength<-function(start,end,autonMap){
  
  # look for first bp in autonMap that is greater than or equal to the start value
  for (x in 1:nrow(autonMap)){
    if (autonMap$V2[x]>=start){
      break
    }
  }
  
  if (autonMap$V2[x]==start){ #if it is equal then read off the cM value
    cMstart<-autonMap$V4[x]
  } else { # if not then adjust linerally between this snp and the previous
    cMstart<-autonMap$V4[x-1]+(autonMap$V4[x]-autonMap$V4[x-1])*(start-autonMap$V2[x-1])/(autonMap$V2[x]-autonMap$V2[x-1])
  }
  
  # look for first bp in autonMap that is less than or equal to the end value
  for (x in nrow(autonMap):1){
    if (autonMap$V2[x]<=end){
      break
    }
  }
  
  if (autonMap$V2[x]==end){ #if it is equal then read off the cM value
    cMend<-autonMap$V4[x]
  } else { # if not then adjust linerally between this snp and the previous
    cMend<-autonMap$V4[x]+(autonMap$V4[x+1]-autonMap$V4[x])*(end-autonMap$V2[x])/(autonMap$V2[x+1]-autonMap$V2[x])
  }
  mapLength<-cMend-cMstart
  return(mapLength)
}
getAutonPlot<-function(autonMap,start,end){
  autonMapPlot<-autonMap[,c(2,4)]
  # look for first bp in autonMap that is greater than or equal to the start value
  for (x in 1:nrow(autonMapPlot)){
    if (autonMapPlot$V2[x]>=start){
      break
    }
  }
  
  if (autonMapPlot$V2[x]==start){ #if it is equal then trim off the previous loci
    autonMapPlot<-autonMapPlot[c(x:nrow(autonMapPlot)),]
  } else { # if not then adjust linerally between this snp and the previous
    cMstart<-autonMapPlot$V4[x-1]+(autonMapPlot$V4[x]-autonMapPlot$V4[x-1])*(start-autonMapPlot$V2[x-1])/(autonMapPlot$V2[x]-autonMapPlot$V2[x-1])
    autonMapPlot<-rbind(c(start,cMstart),autonMapPlot[c(x:nrow(autonMapPlot)),])
  }
  
  # look for first bp in autonMap that is less than or equal to the end value
  for (x in nrow(autonMapPlot):1){
    if (autonMapPlot$V2[x]<=end){
      break
    }
  }
  
  if (autonMap$V2[x]==end){ #if it is equal then trim off the end of the map
    autonMapPlot<-autonMapPlot[c(1:x),]
  } else { # if not then adjust linerally between this snp and the previous
    cMend<-autonMapPlot$V4[x]+(autonMapPlot$V4[x+1]-autonMapPlot$V4[x])*(end-autonMapPlot$V2[x])/(autonMapPlot$V2[x+1]-autonMapPlot$V2[x])
    autonMapPlot<-rbind(autonMapPlot[c(1:x),],c(end,cMend))
  }
  
  #adjust so it starts at 0
  autonMapPlot$V4<-autonMapPlot$V4-autonMapPlot$V4[1]
  return(autonMapPlot)
}

#######################################################################

##get map lengths for each chromosome
rhoLengths<-rep(0,38)
cMLengths<-rep(0,38)

## get counts for SNP density comparison
snpDensity<-data.frame(chr=c(1:38),bp=0,auton_snps=0,ldhat_snps=0)

## Loop for each chromosome
for (i in 1:38){
  ## read in current map
  ldhatMap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/chr",i,"/res.txt"),header=TRUE)
  rhoLength<-ldhatMap$Mean_rho[1]
  ldhatMap<-ldhatMap[-1,]
  
  ## get start and end from vcf
  vcf<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/chr",i,"/cleanVCF.vcf"),header=FALSE,skip=6)
  start<-vcf$V2[1]
  end<-vcf$V2[nrow(vcf)]
  rm(vcf)
  
  ## get map length
  autonMap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/AutonMaps/mark4_cleaned_chr",i,".cf3.1.sorted.txt"),header=FALSE,skip=1)
  mapLength<-getMapLength(start,end,autonMap)
  
  ## convert map to cM
  newMap<-convertLDhatTocMperMb(loci=ldhatMap$Loci,rate=ldhatMap$Mean_rho,start,end,map_length_cM=mapLength)
  
  # write map to file
  write.table(newMap,paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/chr",i,"/cM_map_chr",i,".txt"),quote=FALSE,row.names=FALSE)
  
#  #store map lengths
#  rhoLengths[i]<-rhoLength
#  cMLengths[i]<-mapLength
#  
#  # create plots of ldhat vs Auton
#  options(scipen=999)
#  autonMapPlot<-getAutonPlot(autonMap,start,end)
#  jpeg(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/comparegraphs/map_compare_chr",i,".jpeg"),res=300,height=14,width=14,units='cm',quality=100)
#  plot(newMap$loci/1000000,newMap$cumulative_cM,type="l",xlab="Position (Mb)",ylab="Cumulative cM",main=paste("Chromosome",i))
#  lines(autonMapPlot$V2/1000000,autonMapPlot$V4,col="blue",lty=2)
#  legend("topleft",legend=c("LDhat","Auton"),lty=c(1,2),col=c("black","blue"))
#  dev.off()
#
#  ## get SNP density information
#  snpDensity$bp[i]<-end-start+1
#  snpDensity$auton_snps[i]<-nrow(autonMap)
#  snpDensity$ldhat_snps[i]<-nrow(newMap)

}

##snp denity:
##auton:
#1000000*sum(snpDensity$auton_snps)/sum(snpDensity$bp)
##ldhat:
#1000000*sum(snpDensity$ldhat_snps)/sum(snpDensity$bp)
#
### get scatterplot comparing the two map lengths
#jpeg("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/comparegraphs/compare_map_lengths.jpeg",res=300,height=14,width=14,units='cm',quality=100)
#plot(rhoLengths,cMLengths,xlab="LDhat map length (rho)",ylab="Auton map length (cM)")
#abline(lm(cMLengths~rhoLengths))
#corr<-paste(" =",round(cor(rhoLengths,cMLengths)^2,4))
#legend("bottomright",as.expression(substitute("r"^2*corrr,list(corrr=corr))))
#dev.off()
#cor.test(rhoLengths,cMLengths)$p.value




