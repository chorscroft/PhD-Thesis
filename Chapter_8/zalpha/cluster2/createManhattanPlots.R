## load the qqman package

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

chr_stats<-list()
for (i in 1:38){
  chr_stats[[i]]<-read.table(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/cluster2/zalpha_all_statistics/all_stats_chr_",i,".txt"),header=TRUE,stringsAsFactors = FALSE)
}

SNP<-NULL
CHR<-NULL
BP<-NULL
for (i in 1:38){
  SNP<-c(SNP,paste0(i,"_",chr_stats[[i]]$loci))
  CHR<-c(CHR,rep(i,nrow(chr_stats[[i]])))
  BP<-c(BP,chr_stats[[i]]$loci)
}

 
## loop around stats
for (s in 1:35){  #1:35
  
  STAT<-NULL
  for (i in 1:38){
    STAT<-c(STAT,chr_stats[[i]][,1+s])
  }

  manhatdf<-data.frame(SNP=SNP,CHR=CHR,BP=BP,STAT=STAT)
  manhatdf<-manhatdf[is.na(manhatdf$STAT)==FALSE,]
  
  if (s ==3 | s == 4){
    #highlight<-manhatdf$SNP[rank(manhatdf$STAT)<=100]
    percent0_1<-quantile(manhatdf$STAT,0.001)
  } else {
    #highlight<-manhatdf$SNP[rank(manhatdf$STAT)>=(nrow(manhatdf)-99)]
    percent0_1<-quantile(manhatdf$STAT,0.999)
  }
  atloc<-cumsum(as.numeric(tapply(manhatdf$BP,manhatdf$CHR,max)))-tapply(manhatdf$BP,manhatdf$CHR,quantile,0.5)
  
  jpeg(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/cluster2/manhattan/",statsName[s],".jpg"),res=300,height=14,width=25,units='cm',quality=100)
  par(mar=c(5,7,4,2)+0.1) #Change margins so formula fits on the left
  if (min(manhatdf$STAT,na.rm=TRUE)<0){
    if (min(manhatdf$STAT,na.rm=TRUE)>-1 & max(manhatdf$STAT,na.rm=TRUE)<1)  {
      ylim<-c(floor(min(manhatdf$STAT,na.rm=TRUE)*10)/10,ceiling(max(manhatdf$STAT,na.rm=TRUE)*10)/10)
    } else {    
      ylim<-c(floor(min(manhatdf$STAT,na.rm=TRUE)),ceiling(max(manhatdf$STAT,na.rm=TRUE)))
    }
  } else {
    ylim<-c(0,ceiling(max(manhatdf$STAT,na.rm=TRUE)))
  }
  manhattan(manhatdf, p = "STAT", logp = FALSE, ylab = statsFormula[s], genomewideline = percent0_1, 
            suggestiveline = FALSE, main = statsFormula[s],xaxt="n",ylim=ylim)#,highlight=highlight)
  axis(1,at=atloc[seq(2,38,2)],labels=c(paste0(" ",seq(2,8,2)),seq(10,38,2)),adj=1,hadj=0.6,las=3)
  axis(1,at=atloc[seq(1,37,2)],labels=c(paste0(seq(1,9,2)," "),seq(11,37,2)),adj=1,hadj=1.8,las=3,tcl=-1)
  dev.off()
}





