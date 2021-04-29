
## Define functions for setting up the results data frame and populating it
setUpDataframe<-function(position){
  data.frame(loci=position,
             Zalpha=NA,
             Zbeta=NA,
             L_plus_R=NA,
             LR=NA,
             Zalpha_plus_Zbeta=NA,
             Zalpha_minus_Zbeta=NA,
             ZalphaZbeta=NA,
             Zalpha_over_Zbeta=NA,
             Zalpha_rsq_over_expected=NA,
             Zalpha_log_rsq_over_expected=NA,
             Zalpha_Zscore=NA,
             Zalpha_BetaCDF=NA,
             Zalpha_minus_Zalpha_expected=NA,
             Zalpha_over_Zalpha_expected=NA,
             Zbeta_rsq_over_expected=NA,
             Zbeta_log_rsq_over_expected=NA,
             Zbeta_Zscore=NA,
             Zbeta_BetaCDF=NA,
             Zbeta_minus_Zbeta_expected=NA,
             Zbeta_over_Zbeta_expected=NA,
             Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected=NA,
             Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected=NA,
             Zalpha_Zscore_plus_Zbeta_Zscore=NA,
             Zalpha_BetaCDF_plus_Zbeta_BetaCDF=NA,
             Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected=NA,
             Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected=NA,
             Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected=NA,
             Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected=NA,
             Zalpha_Zscore_minus_Zbeta_Zscore=NA,
             Zalpha_BetaCDF_minus_Zbeta_BetaCDF=NA,
             Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected=NA,
             Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected=NA,
             Zalpha_BetaCDFZbeta_BetaCDF=NA,
             Zalpha_BetaCDF_over_Zbeta_BetaCDF=NA,
             Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected=NA
  )
}
calculateResults<-function(Zalpha_all){
  resultsDataframe<-setUpDataframe(Zalpha_all$position)
  resultsDataframe$Zalpha<-Zalpha_all$Zalpha
  resultsDataframe$Zbeta<-Zalpha_all$Zbeta
  resultsDataframe$L_plus_R<-Zalpha_all$L_plus_R
  resultsDataframe$LR<-Zalpha_all$LR
  resultsDataframe$Zalpha_plus_Zbeta<-Zalpha_all$Zalpha+Zalpha_all$Zbeta
  resultsDataframe$Zalpha_minus_Zbeta<-Zalpha_all$Zalpha-Zalpha_all$Zbeta
  resultsDataframe$ZalphaZbeta<-Zalpha_all$Zalpha*Zalpha_all$Zbeta
  resultsDataframe$Zalpha_over_Zbeta<-Zalpha_all$Zalpha/Zalpha_all$Zbeta
  resultsDataframe$Zalpha_rsq_over_expected<-Zalpha_all$Zalpha_rsq_over_expected
  resultsDataframe$Zalpha_log_rsq_over_expected<-Zalpha_all$Zalpha_log_rsq_over_expected
  resultsDataframe$Zalpha_Zscore<-Zalpha_all$Zalpha_Zscore
  resultsDataframe$Zalpha_BetaCDF<-Zalpha_all$Zalpha_BetaCDF
  resultsDataframe$Zalpha_minus_Zalpha_expected<-Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected
  resultsDataframe$Zalpha_over_Zalpha_expected<-Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected
  resultsDataframe$Zbeta_rsq_over_expected<-Zalpha_all$Zbeta_rsq_over_expected
  resultsDataframe$Zbeta_log_rsq_over_expected<-Zalpha_all$Zbeta_log_rsq_over_expected
  resultsDataframe$Zbeta_Zscore<-Zalpha_all$Zbeta_Zscore
  resultsDataframe$Zbeta_BetaCDF<-Zalpha_all$Zbeta_BetaCDF
  resultsDataframe$Zbeta_minus_Zbeta_expected<-Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected
  resultsDataframe$Zbeta_over_Zbeta_expected<-Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
  resultsDataframe$Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected<-Zalpha_all$Zalpha_rsq_over_expected+Zalpha_all$Zbeta_rsq_over_expected
  resultsDataframe$Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected<-Zalpha_all$Zalpha_log_rsq_over_expected+Zalpha_all$Zbeta_log_rsq_over_expected
  resultsDataframe$Zalpha_Zscore_plus_Zbeta_Zscore<-Zalpha_all$Zalpha_Zscore+Zalpha_all$Zbeta_Zscore
  resultsDataframe$Zalpha_BetaCDF_plus_Zbeta_BetaCDF<-Zalpha_all$Zalpha_BetaCDF+Zalpha_all$Zbeta_BetaCDF
  resultsDataframe$Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected<-(Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)+(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected)
  resultsDataframe$Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected<-Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected+Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
  resultsDataframe$Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected<-Zalpha_all$Zalpha_rsq_over_expected-Zalpha_all$Zbeta_rsq_over_expected
  resultsDataframe$Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected<-Zalpha_all$Zalpha_log_rsq_over_expected-Zalpha_all$Zbeta_log_rsq_over_expected
  resultsDataframe$Zalpha_Zscore_minus_Zbeta_Zscore<-Zalpha_all$Zalpha_Zscore-Zalpha_all$Zbeta_Zscore
  resultsDataframe$Zalpha_BetaCDF_minus_Zbeta_BetaCDF<-Zalpha_all$Zalpha_BetaCDF-Zalpha_all$Zbeta_BetaCDF
  resultsDataframe$Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected<-(Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)-(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected)
  resultsDataframe$Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected<-Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected-Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
  resultsDataframe$Zalpha_BetaCDFZbeta_BetaCDF<-Zalpha_all$Zalpha_BetaCDF*Zalpha_all$Zbeta_BetaCDF
  resultsDataframe$Zalpha_BetaCDF_over_Zbeta_BetaCDF<-Zalpha_all$Zalpha_BetaCDF/Zalpha_all$Zbeta_BetaCDF
  resultsDataframe$Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected<-Zalpha_all$Zalpha/Zalpha_all$Zbeta-Zalpha_all$Zalpha_expected/Zalpha_all$Zbeta_expected
  return(resultsDataframe)
}

## Loop over chromosomes 
for (chr in 1:38){

  Zalpha_all<-read.table(paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr",chr,"/zalpha_chr_",chr,".txt"),header=TRUE, stringsAsFactors=FALSE)
  dfResults<-calculateResults(Zalpha_all)
  write.table(dfResults,paste0("/temp/hgig/EXOME_DATA/Clare/Dogs/cluster2/zalpha_all/chr",chr,"/all_stats_chr_",chr,".txt"),row.names=FALSE,quote=FALSE)  
  
}

