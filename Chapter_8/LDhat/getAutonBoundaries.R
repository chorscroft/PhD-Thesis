
## create data frame for comparison

boundDF<-data.frame(chr=c(1:38),start=NA,end=NA)

for (i in 1:38){
  ## get boudaries of map for each autosome
  auton<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/AutonMaps/mark4_cleaned_chr",i,".cf3.1.sorted.txt"),header=FALSE,skip=1)
  boundDF$start[i]<-auton$V2[1]
  boundDF$end[i]<-auton$V2[nrow(auton)]
  rm(auton)
}

write.table(boundDF,"/temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/DataForZalpha/AutonBoundaries.txt",row.names=FALSE,col.names=FALSE)
