
chrm<-commandArgs(T)[1]

## get the zalpha package
require(zalpha)

## read in data
data<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/DataForZalpha/chr",chrm,".tped"),header=FALSE,stringsAsFactors=FALSE)
data[,5:ncol(data)][data[,5:ncol(data)]==0]<-NA

## get genetic distances
LDmap<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr",chrm,"/cM_map_chr",chrm,".txt"),header=TRUE,stringsAsFactors=FALSE)

## read in LD profile
LDprofile<-read.table("/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDprofile.txt",header=TRUE,stringsAsFactors=FALSE)

all_zalpha<-Zalpha_all(pos=data$V4,ws=300000,x=as.matrix(data[,5:ncol(data)]),dist=LDmap$cumulative_cM,LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,LDprofile_Beta_b=LDprofile$Beta_b)
write.table(all_zalpha,paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/zalpha_all/chr",chrm,"/zalpha_chr_",chrm,".txt"),quote=FALSE,row.names=FALSE)
