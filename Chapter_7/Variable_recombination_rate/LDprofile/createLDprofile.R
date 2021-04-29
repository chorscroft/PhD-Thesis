## Read in the end of the dataset for creating the LDprofile

setwd("/home/ch19g17/Sim_for_zalpha_paper/LDprofile/")

## Load required package
require(fitdistrplus)

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

## set up diff and rsq concatenated vectors
concatDiffs<-NULL
concatRsq<-NULL

## loop around simulations

for (sim in 1:100){

  ## read in simulation output
  
  positions<-as.vector(read.table(paste0("/home/ch19g17/Sim_for_variable_recombination_zalpha/LDprofile/sim_",sim,"/outputEndPositions.txt"),header=FALSE)[,1])
  nosnps<-length(positions)
  
  snps<-read.fwf(paste0("/home/ch19g17/Sim_for_variable_recombination_zalpha/LDprofile/sim_",sim,"/outputEnd.txt"),rep(1,nosnps),header=FALSE)
  snps<-t(as.matrix(snps))
  nogenomes<-dim(snps)[2]
  
  ##remove snps with MAF < 0.05
  keep<-rep(TRUE,nosnps)
  for (i in 1:nosnps){
    if (sum(snps[i,])/nogenomes < 0.05){
      keep[i]<-FALSE
    }
  }
  nosnps<-sum(keep)
  snps<-snps[keep,]
  positions<-positions[keep]
  
  cM<-getCM(positions,cumRecomb)
  
  lower_triangle<-function(x){
    x[lower.tri(x)]
  }
  
  diffs<-lower_triangle(outer(cM,cM,"-"))
  rsq<-lower_triangle(cor(t(snps))^2)
  
  rsq<-rsq[diffs<2]
  diffs<-diffs[diffs<2]
  
  concatDiffs<-c(concatDiffs,diffs)
  concatRsq<-c(concatRsq,rsq)
}
  
### Create LD profile

## Inputs:
## Data frame of pairs of SNPs
## 1: the distance between SNPs
## 2: rsquared values
##
## Output:
## Data frame with columns:
## 1: bin in steps of 0.00001
## 2: Expected rsquared
## 3: standard deviation
## 4: Beta distribution a
## 5: Beta distribution b

assign_bins<-function(bin_size,number){
  ceilingTemp<-ceiling(number/bin_size)
  if(isTRUE(all.equal(ceilingTemp,number/bin_size))){
    return(ceilingTemp*bin_size)
  } else {
    return(floor(number/bin_size)*bin_size)
  }
}
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
getLDprofile<-function(diffs, rsq){
  bin<-assign_bins(0.0001,diffs)
  LDprofile<-data.frame(bin=seq(0,1.9999,0.0001),rsq=NA,sd=NA,Beta_a=NA,Beta_b=NA,n=NA)
  for (i in 1:nrow(LDprofile)){
    LDprofile$n[i]<-sum(bin==LDprofile$bin[i])
    if (LDprofile$n[i]>0){
      temprsq<-rsq[bin==LDprofile$bin[i]]
      LDprofile$rsq[i]<-mean(temprsq)
      LDprofile$sd[i]<-sd(temprsq)
      if (LDprofile$n[i]>1 & LDprofile$sd[i]>0){
        if (sum(temprsq==1 | temprsq==0)>0){
          # if there are any 0s or 1s adjust the data
          temprsq<-(temprsq*(length(temprsq)-1)+0.5)/length(temprsq)
        }
        betafit<-try(fitdist(temprsq,"beta"))
        if (class(betafit) != "try-error"){
          LDprofile$Beta_a[i]<-betafit$estimate[1]
          LDprofile$Beta_b[i]<-betafit$estimate[2]
        } else {
          startBetaParams<-estBetaParams(LDprofile$rsq[i], LDprofile$sd[i]^2)
          betafit<-try(fitdist(temprsq,"beta",start=list(shape1=startBetaParams$alpha, shape2=startBetaParams$beta)))
          if (class(betafit) != "try-error"){
            LDprofile$Beta_a[i]<-betafit$estimate[1]
            LDprofile$Beta_b[i]<-betafit$estimate[2]
          } else {
            LDprofile$Beta_a[i]<-NA
            LDprofile$Beta_b[i]<-NA
          }
        }
      }
    }
  }
  
  return(LDprofile)
}

LDprofile<-getLDprofile(diffs=concatDiffs,rsq=concatRsq)
write.csv(LDprofile,"/home/ch19g17/Sim_for_variable_recombination_zalpha/LDprofile/LDprofile.csv",row.names=FALSE)

