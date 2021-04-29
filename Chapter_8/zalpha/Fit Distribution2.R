## Finds Candidate regions
require(qqman)

## get the formulae for each of the stats
statsFormula<-c(
  expression(paste("Z",''[alpha])),
  expression(paste("Z",''[beta])),
  expression(paste(bgroup("(",atop("|L|",2),")"),"+",bgroup("(",atop("|R|",2),")"))),
  "|L||R|",
  expression(paste("Z",''[alpha],"+","Z",''[beta])), #"Zalpha_plus_Zbeta",
  expression(paste("Z",''[alpha],"-","Z",''[beta])), #"Zalpha_minus_Zbeta",
  expression(paste("Z",''[alpha],"Z",''[beta])), #"ZalphaZbeta",
  expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])))), #"Zalpha_over_Zbeta",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore)), #"Zalpha_Zscore",
  expression(paste("Z"[alpha]^BetaCDF)), #"Zalpha_BetaCDF",
  expression(paste("Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"))), #"Zalpha_minus_Zalpha_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected",
  expression(paste("Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zbeta_rsq_over_expected",
  expression(paste("Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zbeta_log_rsq_over_expected",
  expression(paste("Z"[beta]^ZScore)), #"Zbeta_Zscore",
  expression(paste("Z"[beta]^BetaCDF)), #"Zbeta_BetaCDF",
  expression(paste("Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"))), #"Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"+","Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"+","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore,"+","Z"[beta]^ZScore)), #"Zalpha_Zscore_plus_Zbeta_Zscore",
  expression(paste("Z"[alpha]^BetaCDF,"+","Z"[beta]^BetaCDF)), #"Zalpha_BetaCDF_plus_Zbeta_BetaCDF",
  expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")+","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")), #"Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"+",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"-","Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"-","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore,"-","Z"[beta]^ZScore)), #"Zalpha_Zscore_minus_Zbeta_Zscore",
  expression(paste("Z"[alpha]^BetaCDF,"-","Z"[beta]^BetaCDF)), #"Zalpha_BetaCDF_minus_Zbeta_BetaCDF",
  expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")-","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")), #"Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"-",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF)), #"Zalpha_BetaCDFZbeta_BetaCDF",
  expression(paste(frac("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF))), #"Zalpha_BetaCDF_over_Zbeta_BetaCDF",
  expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])),"-",frac("Z"[alpha]^paste("E[",r^2,"]"),"Z"[beta]^paste("E[",r^2,"]")))) #"Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected"
)
statsName<-c(
  "Zalpha",
  "Zbeta",
  "L_plus_R",
  "LR",
  "Zalpha_plus_Zbeta",
  "Zalpha_minus_Zbeta",
  "ZalphaZbeta",
  "Zalpha_over_Zbeta",
  "Zalpha_rsq_over_expected",
  "Zalpha_log_rsq_over_expected",
  "Zalpha_Zscore",
  "Zalpha_BetaCDF",
  "Zalpha_minus_Zalpha_expected",
  "Zalpha_over_Zalpha_expected",
  "Zbeta_rsq_over_expected",
  "Zbeta_log_rsq_over_expected",
  "Zbeta_Zscore",
  "Zbeta_BetaCDF",
  "Zbeta_minus_Zbeta_expected",
  "Zbeta_over_Zbeta_expected",
  "Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected",
  "Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected",
  "Zalpha_Zscore_plus_Zbeta_Zscore",
  "Zalpha_BetaCDF_plus_Zbeta_BetaCDF",
  "Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected",
  "Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected",
  "Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected",
  "Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected",
  "Zalpha_Zscore_minus_Zbeta_Zscore",
  "Zalpha_BetaCDF_minus_Zbeta_BetaCDF",
  "Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected",
  "Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected",
  "Zalpha_BetaCDFZbeta_BetaCDF",
  "Zalpha_BetaCDF_over_Zbeta_BetaCDF",
  "Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected"
)

## Get data
chr_stats<-list()
for (i in 1:38){
  chr_stats[[i]]<-read.table(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/zalpha_all_statistics/all_stats_chr_",i,".txt"),header=TRUE,stringsAsFactors = FALSE)
}

SNP<-NULL
CHR<-NULL
BP<-NULL
addLeading0s<-function(x,len){
  while (nchar(x)<len){
    x<-paste0("0",x) 
  }
  return(x)
}

for (i in 1:38){
  formatSNPs<-sapply(chr_stats[[i]]$loci,addLeading0s,8)
  SNP<-c(SNP,paste0(addLeading0s(i,2),"_",formatSNPs))
  CHR<-c(CHR,rep(i,nrow(chr_stats[[i]])))
  BP<-c(BP,chr_stats[[i]]$loci)
}

## Fit distribution to Zalpha
s<-1
STAT<-NULL
for (i in 1:38){
  STAT<-c(STAT,chr_stats[[i]][,1+s])
}
manhatdf_Z<-data.frame(SNP=SNP,CHR=CHR,BP=BP,Z=STAT)
manhatdf_Z<-manhatdf_Z[is.na(manhatdf_Z$Z)==FALSE,]

## Calculate empirical distribution
manhatdf_Z$Z_EmpDist<-rank(manhatdf_Z$Z)/nrow(manhatdf_Z)

## Fit distribution to Zalpha_rsq_over_expected
s<-9
STAT<-NULL
for (i in 1:38){
  STAT<-c(STAT,chr_stats[[i]][,1+s])
}
manhatdf_Zroe<-data.frame(SNP=SNP,CHR=CHR,BP=BP,Zroe=STAT)
manhatdf_Zroe<-manhatdf_Zroe[is.na(manhatdf_Zroe$Zroe)==FALSE,]

## Calculate empirical distribution
manhatdf_Zroe$Zroe_EmpDist<-rank(manhatdf_Zroe$Zroe)/nrow(manhatdf_Zroe)

## Fit distribution to Zalpha_BetaCDF
s<-12
STAT<-NULL
for (i in 1:38){
  STAT<-c(STAT,chr_stats[[i]][,1+s])
}
manhatdf_Zb<-data.frame(SNP=SNP,CHR=CHR,BP=BP,Zb=STAT)
manhatdf_Zb<-manhatdf_Zb[is.na(manhatdf_Zb$Zb)==FALSE,]

## Calculate empirical distribution
manhatdf_Zb$Zb_EmpDist<-rank(manhatdf_Zb$Zb)/nrow(manhatdf_Zb)

## Fit distribution to Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected
s<-27
STAT<-NULL
for (i in 1:38){
  STAT<-c(STAT,chr_stats[[i]][,1+s])
}
manhatdf_ZroeMZbroe<-data.frame(SNP=SNP,CHR=CHR,BP=BP,ZroeMZbroe=STAT)
manhatdf_ZroeMZbroe<-manhatdf_ZroeMZbroe[is.na(manhatdf_ZroeMZbroe$ZroeMZbroe)==FALSE,]

## Calculate empirical distribution
manhatdf_ZroeMZbroe$ZroeMZbroe_EmpDist<-rank(manhatdf_ZroeMZbroe$ZroeMZbroe)/nrow(manhatdf_ZroeMZbroe)

## Fit distribution to Zalpha_BetaCDF_over_Zbeta_BetaCDF
s<-34
STAT<-NULL
for (i in 1:38){
  STAT<-c(STAT,chr_stats[[i]][,1+s])
}
manhatdf_ZbOZbb<-data.frame(SNP=SNP,CHR=CHR,BP=BP,ZbOZbb=STAT)
manhatdf_ZbOZbb<-manhatdf_ZbOZbb[is.na(manhatdf_ZbOZbb$ZbOZbb)==FALSE,]

## Calculate empirical distribution
manhatdf_ZbOZbb$ZbOZbb_EmpDist<-rank(manhatdf_ZbOZbb$ZbOZbb)/nrow(manhatdf_ZbOZbb)




## manhattans
# par(mfrow=c(1,1))
# manhattan(manhatdf_Z,p="Z",logp = FALSE,suggestiveline = quantile(manhatdf_Z$Z,0.999),main=statsFormula[1],ylab="")
# manhattan(manhatdf_Zroe,p="Zroe",logp = FALSE,suggestiveline = quantile(manhatdf_Zroe$Zroe,0.999),main=statsFormula[9],ylab="")
# manhattan(manhatdf_Zb,p="Zb",logp = FALSE,suggestiveline = quantile(manhatdf_Zb$Zb,0.999),main=statsFormula[12],ylab="")

## merge datasets
merged<-merge(manhatdf_Z,manhatdf_Zroe,all=TRUE)
merged<-merge(merged,manhatdf_Zb,all=TRUE)

highlightSNPs<-merged$SNP[(merged$Z_EmpDist>0.999 & is.na(merged$Z_EmpDist)==FALSE) | 
                          (merged$Zroe_EmpDist>0.999 & is.na(merged$Zroe_EmpDist)==FALSE) | 
                          (merged$Zb_EmpDist>0.999 & is.na(merged$Zb_EmpDist)==FALSE)]


## plots with highlights
# jpeg("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/tophighlight_new.jpeg",res=300,height=14,width=16,units='cm',quality=100)
#   layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))
#   atloc<-cumsum(as.numeric(tapply(merged$BP,merged$CHR,max)))-tapply(merged$BP,merged$CHR,quantile,0.5)
#   par(mar=c(0,5,2,2)+0.1,oma=c(6,0,0,0))
#   manhattan(merged,p="Z",logp = FALSE,suggestiveline = quantile(merged$Z,0.999,na.rm=TRUE),highlight = highlightSNPs,xaxt="n",xlab="",ylab=statsFormula[1])
#   legend("topright",legend=c("top 0.1% of SNPs","top 0.1% for at least one statistic"),col = c("blue","green3"),lty=c(1,NA),pch=c(NA,16),inset = c(0,-0.2),xpd=TRUE)
#   manhattan(merged[is.na(merged$Zroe)==FALSE,],p="Zroe",logp=FALSE,suggestiveline = quantile(merged$Zroe,0.999,na.rm=TRUE),highlight = highlightSNPs,xaxt="n",xlab="",ylab=statsFormula[9])
#   manhattan(merged[is.na(merged$Zb)==FALSE,],p="Zb",logp=FALSE,suggestiveline = quantile(merged$Zb,0.999,na.rm=TRUE),highlight = highlightSNPs,xaxt="n",xlab="",ylab=statsFormula[12])
#   axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
#   axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
#   mtext("Chromosome",1,3,cex=0.75)
# dev.off()

## Filter by the candidate SNPs (top 0.1% for at least one stat)
candSNPs<-merged[(merged$Z_EmpDist>0.999 & is.na(merged$Z_EmpDist)==FALSE) | 
                       (merged$Zroe_EmpDist>0.999 & is.na(merged$Zroe_EmpDist)==FALSE) | 
                       (merged$Zb_EmpDist>0.999 & is.na(merged$Zb_EmpDist)==FALSE),]

## read in MAF data and merge on
MAF<-read.table("H:/Dogs/Candidate Regions/MAFs.txt",header = TRUE)
MAF<-MAF[,c(2,7,3,4,5)]
candSNPs<-merge(candSNPs,MAF)
candSNPs<-candSNPs[order(candSNPs$SNP),c(3,1,2,4:12)]

write.table(candSNPs[,c(-1)],"H:/Dogs/Candidate Regions/CandSNPs.txt",row.names = FALSE,col.names = TRUE,quote = FALSE)

#merge on final stats
merged<-merge(merged,manhatdf_ZroeMZbroe,all=TRUE)
merged<-merge(merged,manhatdf_ZbOZbb,all=TRUE)

## Write table of the five stats
write.table(merged,"H:/Dogs/Candidate Regions/merged.txt",row.names = FALSE,col.names = TRUE,quote = FALSE)


## Candidate Regions

## get recombination map
chrMap<-list()

for (i in 1:38){
  chrMap[[i]]<-read.table(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/LDhat/Maps/cM_map_chr",i,".txt"),stringsAsFactors = FALSE,header=TRUE)
}
plot(chrMap[[3]]$loci[chrMap[[3]]$loci<20000000],chrMap[[3]]$cMMb[chrMap[[3]]$loci<20000000],type="l")


## First create Zalpha manhattan
jpeg("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/Zalpha.jpeg",res=300,height=10,width=20,units='cm',quality=100)
atloc<-cumsum(as.numeric(tapply(manhatdf_Z$BP,manhatdf_Z$CHR,max)))-tapply(manhatdf_Z$BP,manhatdf_Z$CHR,quantile,0.5)
manhattan(manhatdf_Z,p="Z",logp = FALSE,suggestiveline = quantile(manhatdf_Z$Z,0.999),ylab=statsFormula[1],xaxt="n",xlab="",ylim=c(0,0.8))
axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
mtext("Chromosome",1,3,cex=0.75)
legend("topright",legend = "top 0.1% of SNPs",col="blue",lty=1,bty="n")
dev.off()

## zoom in on chromosome 3 bp 0-2000000
chr3plot<-function(manhat_BP,manhat_CHR,manhat_STAT,ylim,statno,last){
  if (last==FALSE){
    plot(manhat_BP[manhat_CHR==3 & manhat_BP<20000000],manhat_STAT[manhat_CHR==3 & manhat_BP<20000000],xaxt="n",xlab="",ylab=statsFormula[statno],ylim=ylim,type="n",las=1)
  } else {
    plot(manhat_BP[manhat_CHR==3 & manhat_BP<20000000],manhat_STAT[manhat_CHR==3 & manhat_BP<20000000],ylab=statsFormula[statno],ylim=ylim,type="n",las=1)
    mtext("Chromosome 3 bp", side = 1, line = 3, outer = TRUE)
  }
  adj=75000
  rect(xleft=2602488-adj, ybottom=-1, xright=2602488+adj, ytop=5,density=NULL,col="lightgreen",border=NA)
  rect(xleft=13415724-adj, ybottom=-1, xright=13415724+adj, ytop=5,density=NULL,col="lightblue",border=NA)
  rect(xleft=17474659-adj, ybottom=-1, xright=17554116+adj, ytop=5,density=NULL,col="pink",border=NA)
  points(manhat_BP[manhat_CHR==3 & manhat_BP<20000000],manhat_STAT[manhat_CHR==3 & manhat_BP<20000000],ylim=ylim,pch=16,cex=0.5)
  abline(h=quantile(manhat_STAT,0.999),col="blue")
  text(x=19500000,y=quantile(manhat_STAT,0.999),adj=c(0,-0.6),labels = "top 0.1%",col = "blue")
  mtext("A",side=3,line=0,at=2602488,col="lightgreen",cex=0.75)
  mtext("B",side=3,line=0,at=13415724,col="lightblue",cex=0.75)
  mtext("C",side=3,line=0,at=17514388,col="pink",cex=0.75)
  box()
}

jpeg("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/chr320MB.jpeg",res=300,height=10,width=20,units='cm',quality=100)
options(scipen=999)
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))
par(mar=c(0,5,1,2)+0.1,oma=c(6,0,0,0))
chr3plot(manhatdf_Z$BP,manhatdf_Z$CHR,manhatdf_Z$Z,c(0,0.6),1,FALSE)
chr3plot(manhatdf_Zroe$BP,manhatdf_Zroe$CHR,manhatdf_Zroe$Zroe,c(0,4),9,FALSE)
chr3plot(manhatdf_Zb$BP,manhatdf_Zb$CHR,manhatdf_Zb$Zb,c(0.2,1),12,TRUE)
dev.off()

## create Zalpha_roe and Z_b manhattans
jpeg("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/Zroe_Zb_manhat.jpeg",res=300,height=18,width=25,units='cm',quality=100)
layout(matrix(c(1,2), 2, 1, byrow = TRUE))
atloc<-cumsum(as.numeric(tapply(manhatdf_Z$BP,manhatdf_Z$CHR,max)))-tapply(manhatdf_Z$BP,manhatdf_Z$CHR,quantile,0.5)
par(mar=c(0,5,2,2)+0.1,oma=c(6,0,0,0))
manhattan(merged[is.na(merged$Zroe)==FALSE,],p="Zroe",logp=FALSE,suggestiveline = quantile(merged$Zroe,0.999,na.rm=TRUE),xaxt="n",xlab="",ylab=statsFormula[9])
legend("topright",legend = "top 0.1% of SNPs",col="blue",lty=1,bty="n")
manhattan(merged[is.na(merged$Zb)==FALSE,],p="Zb",logp=FALSE,suggestiveline = quantile(merged$Zb,0.999,na.rm=TRUE),xaxt="n",xlab="",ylab=statsFormula[12])
axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
mtext("Chromosome",1,3,cex=0.75)
dev.off()





### Get only SNPs for Zroe and Zb
finalCandSNPs <- candSNPs[candSNPs$Zroe_EmpDist>0.999 | candSNPs$Zb_EmpDist>0.999,]
write.table(finalCandSNPs[,c(-1)],"H:/Dogs/Candidate Regions/FinalCandSNPs.txt",row.names = FALSE,col.names = TRUE,quote = FALSE)





## Annotate SNPs using ChIP
require(ChIPpeakAnno)
require(biomaRt)
ensembl<-useMart("ensembl",dataset = "clfamiliaris_gene_ensembl")

testdf<-data.frame(chr=paste0("chr",finalCandSNPs$CHR),start=finalCandSNPs$BP,end=finalCandSNPs$BP,strand=rep("*",nrow(finalCandSNPs)))
DogGrange<-makeGRangesFromDataFrame(testdf)

anno<-annotatePeakInBatch(DogGrange,mart=ensembl,select = "all",output="both")
annogene<-addGeneIDs(anno,mart=ensembl,IDs2Add=c("external_gene_name","go_id","name_1006","namespace_1003"))

write.csv(annogene,"\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/ChIPoutput.csv",quote = FALSE)
