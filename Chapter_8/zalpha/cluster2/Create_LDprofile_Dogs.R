
require(zalpha)

dist<-NULL
x<-NULL
for (i in 1:38){
  cM_map<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/chr",i,"/cM_map_chr",i,".txt"),header=TRUE,stringsAsFactors=FALSE)
  tped<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/chr",i,".tped"),header=FALSE,stringsAsFactors=FALSE)
  dist[[i]]<-cM_map$cumulative_cM
  x[[i]]<-as.matrix(tped[,-c(1:4)])
  # Set missing values to NA
  x[[i]][x[[i]] == 0]<-NA
}

LDprofile<-create_LDprofile(dist=dist,x=x,bin_size=0.0001,max_dist=2,beta_params=TRUE)
options(scipen=999)
write.table(LDprofile,file="/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/LDprofile/LDprofile.txt",quote=FALSE,row.names=FALSE )
