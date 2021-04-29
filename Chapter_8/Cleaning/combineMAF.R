
for (i in 1:38){
  maf<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/MAF/chr",i,"/plink.frq"),header=TRUE)
  map<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/MAF/chr",i,"/plink.map"))
  colnames(map)<-c("CHR","SNP","CM","BP")
  merge<-merge(maf,map[,c(2,4)],by="SNP")
  write.table(merge,"/temp/hgig/EXOME_DATA/Clare/Dogs/MAF/MAFs.txt",append=TRUE,quote=FALSE,row.names=FALSE,col.names=ifelse(i==1,TRUE,FALSE))
}