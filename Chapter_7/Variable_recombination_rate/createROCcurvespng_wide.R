## Get required package
require(pROC)

## Set working directory

setwd("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/VariableRecombinationRates")

## Get data

neutralResults<-read.csv("neutralResults_wide.csv")
selectedResults<-read.csv("selectedResults_wide.csv")
selectedMidResults<-read.csv("selectedMidResults_wide.csv")

## Get area under the curve and pAUC tables

aucDF<-data.frame(stat=colnames(neutralResults[,-1]),AUC=NA,pAUC=NA)
for (i in 2:36){
  rocStat<-roc(controls=neutralResults[,i],cases=selectedResults[,i])
  aucDF$AUC[i-1]<-auc(rocStat)
  aucDF$AUClowCI[i-1]<-ci.auc(rocStat)[1]
  aucDF$AUChighCI[i-1]<-ci.auc(rocStat)[3]
  aucDF$pAUC[i-1]<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  aucDF$pAUClowCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[1]/0.05
  aucDF$pAUChighCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[3]/0.05
}

aucMidNeutDF<-data.frame(stat=colnames(neutralResults[,-1]),AUC=NA,pAUC=NA)
for (i in 2:36){
  rocStat<-roc(controls=neutralResults[,i],cases=selectedMidResults[,i])
  aucMidNeutDF$AUC[i-1]<-auc(rocStat)
  aucMidNeutDF$AUClowCI[i-1]<-ci.auc(rocStat)[1]
  aucMidNeutDF$AUChighCI[i-1]<-ci.auc(rocStat)[3]
  aucMidNeutDF$pAUC[i-1]<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  aucMidNeutDF$pAUClowCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[1]/0.05
  aucMidNeutDF$pAUChighCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[3]/0.05
}
aucMidEndDF<-data.frame(stat=colnames(neutralResults[,-1]),AUC=NA,pAUC=NA)
for (i in 2:36){
  rocStat<-roc(controls=selectedResults[,i],cases=selectedMidResults[,i])
  aucMidEndDF$AUC[i-1]<-auc(rocStat)
  aucMidEndDF$AUClowCI[i-1]<-ci.auc(rocStat)[1]
  aucMidEndDF$AUChighCI[i-1]<-ci.auc(rocStat)[3]
  aucMidEndDF$pAUC[i-1]<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  aucMidEndDF$pAUClowCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[1]/0.05
  aucMidEndDF$pAUChighCI[i-1]<- ci.auc(rocStat,partial.auc=c(0.95,1),reuse.auc=FALSE)[3]/0.05
}

## Write out AUC tables to .csv
setwd("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/VariableRecombinationRates/wide")
write.csv(aucDF,"auc_Selected_Neutral_png_wide.csv")
write.csv(aucMidNeutDF,"auc_MidSelected_Neutral_png_wide.csv")
write.csv(aucMidEndDF,"auc_MidSelected_EndSelected_png_wide.csv")

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

## Create roc curves for each statistic

setwd("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/VariableRecombinationRates/wide/rocNeutralSelected")
for (i in 1:35){
  png(paste0(colnames(neutralResults)[i+1],".png"),res=300,width=14,height=14,units="cm")
  rocStat<-roc(controls=neutralResults[,i+1],cases=selectedResults[,i+1])
  par(pty="s")
  plot(0,xlim=c(1,0),ylim=c(0,1),type="n",xlab="Specificity",ylab="Sensitivity",main=statsFormula[i])
  AUC<-auc(rocStat)
  pAUC<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  legend("bottomright",legend=c(paste0("AUC = ",round(AUC*100,1),"%"),paste0("pAUC = ",round(pAUC*100,1),"%")))
  plot(ci.sp(rocStat,boot.n=2000,sensitivities = seq(0, 1, .01),progress="none"),type="s",lty=2,no.roc=TRUE)
  abline(v=0.95,col="gray")
  abline(a=1,b=-1,col="gray")
  plot(rocStat,add=TRUE)
  dev.off()
}
setwd("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/VariableRecombinationRates/wide/rocNeutralMidSelected")
for (i in 1:35){
  png(paste0(colnames(neutralResults)[i+1],".png"),res=300,width=14,height=14,units="cm")
  rocStat<-roc(controls=neutralResults[,i+1],cases=selectedMidResults[,i+1])
  par(pty="s")
  plot(0,xlim=c(1,0),ylim=c(0,1),type="n",xlab="Specificity",ylab="Sensitivity",main=statsFormula[i])
  AUC<-auc(rocStat)
  pAUC<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  legend("bottomright",legend=c(paste0("AUC = ",round(AUC*100,1),"%"),paste0("pAUC = ",round(pAUC*100,1),"%")))
  plot(ci.sp(rocStat,boot.n=2000,sensitivities = seq(0, 1, .01),progress="none"),type="s",lty=2,no.roc=TRUE)
  abline(v=0.95,col="gray")
  abline(a=1,b=-1,col="gray")
  plot(rocStat,add=TRUE)
  dev.off()
}
setwd("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/VariableRecombinationRates/wide/rocMidSelectedEndSelected")
for (i in 1:35){
  png(paste0(colnames(selectedResults)[i+1],".png"),res=300,width=14,height=14,units="cm")
  rocStat<-roc(controls=selectedMidResults[,i+1],cases=selectedResults[,i+1])
  par(pty="s")
  plot(0,xlim=c(1,0),ylim=c(0,1),type="n",xlab="Specificity",ylab="Sensitivity",main=statsFormula[i])
  AUC<-auc(rocStat)
  pAUC<-auc(rocStat,partial.auc=c(0.95,1))/0.05
  legend("bottomright",legend=c(paste0("AUC = ",round(AUC*100,1),"%"),paste0("pAUC = ",round(pAUC*100,1),"%")))
  plot(ci.sp(rocStat,boot.n=2000,sensitivities = seq(0, 1, .01),progress="none"),type="s",lty=2,no.roc=TRUE)
  abline(v=0.95,col="gray")
  abline(a=1,b=-1,col="gray")
  plot(rocStat,add=TRUE)
  dev.off()
}