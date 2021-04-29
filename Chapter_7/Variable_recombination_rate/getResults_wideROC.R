
## get zalpha package
require(zalpha)

## read in LDprofile
LDprofile<-read.csv("/home/ch19g17/Sim_for_variable_recombination_zalpha/LDprofile/LDprofile.csv")

## Read in variable recombination rate
recombination<-read.table("/home/ch19g17/Sim_for_variable_recombination_zalpha/recombinationRates.txt",col.names=c("end","rate"))
cumRecomb<-data.frame(pos=c(0,recombination$end),rate=c(0,recombination$rate),cM=0)
for (i in 2:nrow(cumRecomb)){
  cumRecomb$cM[i]<-(cumRecomb$pos[i]-cumRecomb$pos[i-1])*cumRecomb$rate[i]/0.01+cumRecomb$cM[i-1]
}

## function for getting the centimorgan distances given the cumulative recombination rate and bp positions
getCM<-function(positions,cumRecomb){
  cM<-NULL
  j<-1
  for(i in 1:length(positions)){
    # find upper bound j (also assumes positions are sorted)
    while(positions[i]>=cumRecomb$pos[j]){
      j<-j+1
    }
    # calculate cM
    calcCM<-cumRecomb$cM[j-1]+(positions[i]-cumRecomb$pos[j-1])*cumRecomb$rate[j]/0.01 
    # append to vector
    cM<-c(cM,calcCM)
  }
  return(cM)
}

## Define functions for setting up the results data frame and populating it
setUpDataframe<-function(){
  data.frame(sim=1:100,
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
calculateResults<-function(resultsDataframe,Zalpha_all,sim){
  
  resultsDataframe$Zalpha[sim]<-max(Zalpha_all$Zalpha, na.rm=TRUE)
  resultsDataframe$Zbeta[sim]<-max(Zalpha_all$Zbeta, na.rm=TRUE)
  resultsDataframe$L_plus_R[sim]<-min(Zalpha_all$L_plus_R, na.rm=TRUE)
  resultsDataframe$LR[sim]<-min(Zalpha_all$LR, na.rm=TRUE)
  resultsDataframe$Zalpha_plus_Zbeta[sim]<-max(Zalpha_all$Zalpha+Zalpha_all$Zbeta, na.rm=TRUE)
  resultsDataframe$Zalpha_minus_Zbeta[sim]<-max(Zalpha_all$Zalpha-Zalpha_all$Zbeta, na.rm=TRUE)
  resultsDataframe$ZalphaZbeta[sim]<-max(Zalpha_all$Zalpha*Zalpha_all$Zbeta, na.rm=TRUE)
  resultsDataframe$Zalpha_over_Zbeta[sim]<-max(Zalpha_all$Zalpha/Zalpha_all$Zbeta, na.rm=TRUE)
  resultsDataframe$Zalpha_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_log_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_log_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_Zscore[sim]<-max(Zalpha_all$Zalpha_Zscore, na.rm=TRUE)
  resultsDataframe$Zalpha_BetaCDF[sim]<-max(Zalpha_all$Zalpha_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zalpha_minus_Zalpha_expected[sim]<-max(Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_over_Zalpha_expected[sim]<-max(Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected, na.rm=TRUE)
  resultsDataframe$Zbeta_rsq_over_expected[sim]<-max(Zalpha_all$Zbeta_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zbeta_log_rsq_over_expected[sim]<-max(Zalpha_all$Zbeta_log_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zbeta_Zscore[sim]<-max(Zalpha_all$Zbeta_Zscore, na.rm=TRUE)
  resultsDataframe$Zbeta_BetaCDF[sim]<-max(Zalpha_all$Zbeta_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zbeta_minus_Zbeta_expected[sim]<-max(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected, na.rm=TRUE)
  resultsDataframe$Zbeta_over_Zbeta_expected[sim]<-max(Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_rsq_over_expected+Zalpha_all$Zbeta_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_log_rsq_over_expected+Zalpha_all$Zbeta_log_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_Zscore_plus_Zbeta_Zscore[sim]<-max(Zalpha_all$Zalpha_Zscore+Zalpha_all$Zbeta_Zscore, na.rm=TRUE)
  resultsDataframe$Zalpha_BetaCDF_plus_Zbeta_BetaCDF[sim]<-max(Zalpha_all$Zalpha_BetaCDF+Zalpha_all$Zbeta_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected[sim]<-max((Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)+(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected), na.rm=TRUE)
  resultsDataframe$Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected[sim]<-max(Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected+Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_rsq_over_expected-Zalpha_all$Zbeta_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected[sim]<-max(Zalpha_all$Zalpha_log_rsq_over_expected-Zalpha_all$Zbeta_log_rsq_over_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_Zscore_minus_Zbeta_Zscore[sim]<-max(Zalpha_all$Zalpha_Zscore-Zalpha_all$Zbeta_Zscore, na.rm=TRUE)
  resultsDataframe$Zalpha_BetaCDF_minus_Zbeta_BetaCDF[sim]<-max(Zalpha_all$Zalpha_BetaCDF-Zalpha_all$Zbeta_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected[sim]<-max((Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)-(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected), na.rm=TRUE)
  resultsDataframe$Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected[sim]<-max(Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected-Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected, na.rm=TRUE)
  resultsDataframe$Zalpha_BetaCDFZbeta_BetaCDF[sim]<-max(Zalpha_all$Zalpha_BetaCDF*Zalpha_all$Zbeta_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zalpha_BetaCDF_over_Zbeta_BetaCDF[sim]<-max(Zalpha_all$Zalpha_BetaCDF/Zalpha_all$Zbeta_BetaCDF, na.rm=TRUE)
  resultsDataframe$Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected[sim]<-max(Zalpha_all$Zalpha/Zalpha_all$Zbeta-Zalpha_all$Zalpha_expected/Zalpha_all$Zbeta_expected, na.rm=TRUE)
  return(resultsDataframe)
}

## Set up dataframes to store the max or min of the statistics for the centre region
selectedResults<-setUpDataframe()
selectedMidResults<-setUpDataframe()
neutralResults<-setUpDataframe()

## Loop over the number of simulations (100)

for (sim in 1:100){

  ## Read in the end of the selected dataset
  
  setwd(paste0("/home/ch19g17/Sim_for_variable_recombination_zalpha/selected/sim_",sim))
  
  ## read in simulation output
  
  positions<-as.vector(read.table("outputEndPositions.txt",header=FALSE)[,1])
  nosnps<-length(positions)
  
  snps<-read.fwf("outputEnd.txt",rep(1,nosnps),header=FALSE)
  snps<-t(as.matrix(snps))
  nogenomes<-dim(snps)[2]
  
  ## remove snps with MAF < 0.05
  keep<-rep(TRUE,nosnps)
  for (i in 1:nosnps){
    if (sum(snps[i,])/nogenomes < 0.05){
      keep[i]<-FALSE
    }
  }
  nosnps<-sum(keep)
  snps<-snps[keep,]
  positions<-positions[keep]
  
  ## get cM distance
  cM<-getCM(positions,cumRecomb)
  
  endZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
                            LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
                            LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
                            LDprofile_Beta_b=LDprofile$Beta_b,X=c(200000,800000))
  selectedResults<-calculateResults(selectedResults,endZalpha_all,sim)
  
  ## Read in the Mid of the selected dataset
  
  ## read in simulation output
  positions<-as.vector(read.table("outputMidPositions.txt",header=FALSE)[,1])
  nosnps<-length(positions)
  
  snps<-read.fwf("outputMid.txt",rep(1,nosnps),header=FALSE,skip=3)
  snps<-t(as.matrix(snps))
  nogenomes<-dim(snps)[2]
  
  ## remove snps with MAF < 0.05
  keep<-rep(TRUE,nosnps)
  for (i in 1:nosnps){
    if (sum(snps[i,])/nogenomes < 0.05){
      keep[i]<-FALSE
    }
  }
  nosnps<-sum(keep)
  snps<-snps[keep,]
  positions<-positions[keep]
  
  ## get cM distance
  cM<-getCM(positions,cumRecomb)
  
  midZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
                                       LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
                                       LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
                                       LDprofile_Beta_b=LDprofile$Beta_b,X=c(200000,800000))
  selectedMidResults<-calculateResults(selectedMidResults,midZalpha_all,sim)
  
  ## Read in the end of the neutral dataset
  
  setwd(paste0("/home/ch19g17/Sim_for_variable_recombination_zalpha/neutral/sim_",sim))
  
  ## read in simulation output
  
  positions<-as.vector(read.table("outputEndPositions.txt",header=FALSE)[,1])
  nosnps<-length(positions)
  
  snps<-read.fwf("outputEnd.txt",rep(1,nosnps),header=FALSE)
  snps<-t(as.matrix(snps))
  nogenomes<-dim(snps)[2]
  
  ## remove snps with MAF < 0.05
  keep<-rep(TRUE,nosnps)
  for (i in 1:nosnps){
    if (sum(snps[i,])/nogenomes < 0.05){
      keep[i]<-FALSE
    }
  }
  nosnps<-sum(keep)
  snps<-snps[keep,]
  positions<-positions[keep]
  
  ## get cM distance
  cM<-getCM(positions,cumRecomb)
  
  neutralZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
                            LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
                            LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
                            LDprofile_Beta_b=LDprofile$Beta_b,X=c(200000,800000))
  neutralResults<-calculateResults(neutralResults,neutralZalpha_all,sim)
}

## Write out results to .csv

write.csv(selectedResults,"/home/ch19g17/Sim_for_variable_recombination_zalpha/selectedResults_wide.csv",row.names=FALSE)
write.csv(selectedMidResults,"/home/ch19g17/Sim_for_variable_recombination_zalpha/selectedMidResults_wide.csv",row.names=FALSE)
write.csv(neutralResults,"/home/ch19g17/Sim_for_variable_recombination_zalpha/neutralResults_wide.csv",row.names=FALSE)

