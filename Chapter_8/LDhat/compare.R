
## create data frame for comparison

compareDF<-data.frame(chr=c(1:38),start=NA ,end=NA,SNPs=NA,res_count=NA,locs_n=NA,locs_topend=NA,locs_count=NA,locs_end=NA,start_a=NA,end_a=NA)

for (i in 1:38){

  compareDF$chr[i]<-i
  ## get vcf
  vcf<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",i,"/cleanVCF.vcf"),header=FALSE,skip=6)
  compareDF$start[i]<-vcf$V2[1]
  compareDF$end[i]<-vcf$V2[nrow(vcf)]
  compareDF$SNPs[i]<-nrow(vcf)
  rm(vcf)
  ## get res
  res<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",i,"/res.txt"),header=TRUE)
  compareDF$res_count[i]<-nrow(res)
  rm(res)
  ##get locs
  locs<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",i,"/locs.txt"),header=FALSE,nrows=1)
  compareDF$locs_n[i]=locs$V1
  compareDF$locs_topend[i]=locs$V2
  rm(locs)
  locs2<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",i,"/locs.txt"),header=FALSE,skip=1)
  compareDF$locs_count[i]=nrow(locs2)
  compareDF$locs_end[i]=locs2$V1[nrow(locs2)]
  rm(locs2)
  ## get Auton
  auton<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/AutonMaps/mark4_cleaned_chr",i,".cf3.1.sorted.txt"),header=FALSE,skip=1)
  compareDF$start_a[i]<-auton$V2[1]
  compareDF$end_a[i]<-auton$V2[nrow(auton)]
  rm(auton)
}

for (i in 1:38){
  compareDF$within[i]<-compareDF$start_a[i]<=compareDF$start[i] & compareDF$end_a[i] >= compareDF$end[i]
}

write.table(compareDF,file="comparison.txt",row.names=FALSE,quote=FALSE)
write.table(compareDF[,c(10,11)],"AutonBoundaries.txt",row.names=FALSE,col.names=FALSE)
