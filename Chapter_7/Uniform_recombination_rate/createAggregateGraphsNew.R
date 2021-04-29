# ## get zalpha package
# require(zalpha)
# 
# ## read in LDprofile
# LDprofile<-read.csv("/home/ch19g17/Sim_for_zalpha_paper/LDprofile/LDprofile.csv")
# 
# ## Loop over the number of simulations (100)
# 
# selectedResultsList<-list()
# selectedMidResultsList<-list()
# neutralResultsList<-list()
# 
# ## Define functions for setting up the results data frame and populating it
# calculateResults<-function(Zalpha_all){
#   resultsList<-list()
#   resultsList$position=Zalpha_all$position
#   resultsList$Zalpha=Zalpha_all$Zalpha
#   resultsList$Zbeta=Zalpha_all$Zbeta
#   resultsList$L_plus_R=Zalpha_all$L_plus_R
#   resultsList$LR=Zalpha_all$LR
#   resultsList$Zalpha_plus_Zbeta=Zalpha_all$Zalpha+Zalpha_all$Zbeta
#   resultsList$Zalpha_minus_Zbeta=Zalpha_all$Zalpha-Zalpha_all$Zbeta
#   resultsList$ZalphaZbeta=Zalpha_all$Zalpha*Zalpha_all$Zbeta
#   resultsList$Zalpha_over_Zbeta=Zalpha_all$Zalpha/Zalpha_all$Zbeta
#   resultsList$Zalpha_rsq_over_expected=Zalpha_all$Zalpha_rsq_over_expected
#   resultsList$Zalpha_log_rsq_over_expected=Zalpha_all$Zalpha_log_rsq_over_expected
#   resultsList$Zalpha_Zscore=Zalpha_all$Zalpha_Zscore
#   resultsList$Zalpha_BetaCDF=Zalpha_all$Zalpha_BetaCDF
#   resultsList$Zalpha_minus_Zalpha_expected=Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected
#   resultsList$Zalpha_over_Zalpha_expected=Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected
#   resultsList$Zbeta_rsq_over_expected=Zalpha_all$Zbeta_rsq_over_expected
#   resultsList$Zbeta_log_rsq_over_expected=Zalpha_all$Zbeta_log_rsq_over_expected
#   resultsList$Zbeta_Zscore=Zalpha_all$Zbeta_Zscore
#   resultsList$Zbeta_BetaCDF=Zalpha_all$Zbeta_BetaCDF
#   resultsList$Zbeta_minus_Zbeta_expected=Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected
#   resultsList$Zbeta_over_Zbeta_expected=Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
#   resultsList$Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected=Zalpha_all$Zalpha_rsq_over_expected+Zalpha_all$Zbeta_rsq_over_expected
#   resultsList$Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected=Zalpha_all$Zalpha_log_rsq_over_expected+Zalpha_all$Zbeta_log_rsq_over_expected
#   resultsList$Zalpha_Zscore_plus_Zbeta_Zscore=Zalpha_all$Zalpha_Zscore+Zalpha_all$Zbeta_Zscore
#   resultsList$Zalpha_BetaCDF_plus_Zbeta_BetaCDF=Zalpha_all$Zalpha_BetaCDF+Zalpha_all$Zbeta_BetaCDF
#   resultsList$Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected=(Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)+(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected)
#   resultsList$Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected=Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected+Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
#   resultsList$Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected=Zalpha_all$Zalpha_rsq_over_expected-Zalpha_all$Zbeta_rsq_over_expected
#   resultsList$Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected=Zalpha_all$Zalpha_log_rsq_over_expected-Zalpha_all$Zbeta_log_rsq_over_expected
#   resultsList$Zalpha_Zscore_minus_Zbeta_Zscore=Zalpha_all$Zalpha_Zscore-Zalpha_all$Zbeta_Zscore
#   resultsList$Zalpha_BetaCDF_minus_Zbeta_BetaCDF=Zalpha_all$Zalpha_BetaCDF-Zalpha_all$Zbeta_BetaCDF
#   resultsList$Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected=(Zalpha_all$Zalpha-Zalpha_all$Zalpha_expected)-(Zalpha_all$Zbeta-Zalpha_all$Zbeta_expected)
#   resultsList$Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected=Zalpha_all$Zalpha/Zalpha_all$Zalpha_expected-Zalpha_all$Zbeta/Zalpha_all$Zbeta_expected
#   resultsList$Zalpha_BetaCDFZbeta_BetaCDF=Zalpha_all$Zalpha_BetaCDF*Zalpha_all$Zbeta_BetaCDF
#   resultsList$Zalpha_BetaCDF_over_Zbeta_BetaCDF=Zalpha_all$Zalpha_BetaCDF/Zalpha_all$Zbeta_BetaCDF
#   resultsList$Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected=Zalpha_all$Zalpha/Zalpha_all$Zbeta-Zalpha_all$Zalpha_expected/Zalpha_all$Zbeta_expected
#   return(resultsList)
# }
# 
# for (sim in 1:100){
#   
#   ## Read in the end of the selected dataset
#   
#   setwd(paste0("/home/ch19g17/Sim_for_zalpha_paper/selected/sim_",sim))
#   
#   ## read in simulation output
#   
#   positions<-as.vector(read.table("outputEndPositions.txt",header=FALSE)[,1])
#   nosnps<-length(positions)
#   
#   snps<-read.fwf("outputEnd.txt",rep(1,nosnps),header=FALSE)
#   snps<-t(as.matrix(snps))
#   nogenomes<-dim(snps)[2]
#   
#   ## remove snps with MAF < 0.05
#   keep<-rep(TRUE,nosnps)
#   for (i in 1:nosnps){
#     if (sum(snps[i,])/nogenomes < 0.05){
#       keep[i]<-FALSE
#     }
#   }
#   nosnps<-sum(keep)
#   snps<-snps[keep,]
#   positions<-positions[keep]
#   
#   ## get cM distance
#   cM<-positions*1e-6
#   
#   endZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
#                             LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
#                             LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
#                             LDprofile_Beta_b=LDprofile$Beta_b)
#   
#   selectedResultsList[[paste0("sim_",sim)]]<-calculateResults(endZalpha_all)
#   write.csv(selectedResultsList[[paste0("sim_",sim)]],"simResults.csv",row.names=FALSE)
#   
#   ## Read in the Mid of the selected dataset
#   
#   ## read in simulation output
#   positions<-as.vector(read.table("outputMidPositions.txt",header=FALSE)[,1])
#   nosnps<-length(positions)
#   
#   snps<-read.fwf("outputMid.txt",rep(1,nosnps),header=FALSE,skip=3)
#   snps<-t(as.matrix(snps))
#   nogenomes<-dim(snps)[2]
#   
#   ## remove snps with MAF < 0.05
#   keep<-rep(TRUE,nosnps)
#   for (i in 1:nosnps){
#     if (sum(snps[i,])/nogenomes < 0.05){
#       keep[i]<-FALSE
#     }
#   }
#   nosnps<-sum(keep)
#   snps<-snps[keep,]
#   positions<-positions[keep]
#   
#   ## get cM distance
#   cM<-positions*1e-6
#   
#   midZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
#                             LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
#                             LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
#                             LDprofile_Beta_b=LDprofile$Beta_b)
#   selectedMidResultsList[[paste0("sim_",sim)]]<-calculateResults(midZalpha_all)  
#   write.csv(selectedMidResultsList[[paste0("sim_",sim)]],"simMidResults.csv",row.names=FALSE)
#   
#   ## Read in the end of the neutral dataset
#   
#   setwd(paste0("/home/ch19g17/Sim_for_zalpha_paper/neutral/sim_",sim))
#   
#   ## read in simulation output
#   
#   positions<-as.vector(read.table("outputEndPositions.txt",header=FALSE)[,1])
#   nosnps<-length(positions)
#   
#   snps<-read.fwf("outputEnd.txt",rep(1,nosnps),header=FALSE)
#   snps<-t(as.matrix(snps))
#   nogenomes<-dim(snps)[2]
#   
#   ## remove snps with MAF < 0.05
#   keep<-rep(TRUE,nosnps)
#   for (i in 1:nosnps){
#     if (sum(snps[i,])/nogenomes < 0.05){
#       keep[i]<-FALSE
#     }
#   }
#   nosnps<-sum(keep)
#   snps<-snps[keep,]
#   positions<-positions[keep]
#   
#   ## get cM distance
#   cM<-positions*1e-6
#   
#   neutralZalpha_all<-Zalpha_all(pos=positions,x=snps,ws=200000,dist=cM,
#                                 LDprofile_bins=LDprofile$bin,LDprofile_rsq=LDprofile$rsq,
#                                 LDprofile_sd=LDprofile$sd,LDprofile_Beta_a=LDprofile$Beta_a,
#                                 LDprofile_Beta_b=LDprofile$Beta_b)
#   neutralResultsList[[paste0("sim_",sim)]]<-calculateResults(neutralZalpha_all)
#   write.csv(neutralResultsList[[paste0("sim_",sim)]],"simResults.csv",row.names=FALSE)
# }

########################################################################################

## Read in.csvs with simResults
selectedResultsList<-list()
selectedMidResultsList<-list()
neutralResultsList<-list()
for (sim in 1:100){
  tempFile<-read.csv(paste0("/home/ch19g17/Sim_for_zalpha_paper/neutral/sim_",sim,"/simResults.csv"))
  neutralResultsList[[paste0("sim_",sim)]]<-list()
  for (i in 1:36){
    neutralResultsList[[paste0("sim_",sim)]][[colnames(tempFile)[i]]]<-tempFile[,i]
  }
  tempFile<-read.csv(paste0("/home/ch19g17/Sim_for_zalpha_paper/selected/sim_",sim,"/simResults.csv"))
  selectedResultsList[[paste0("sim_",sim)]]<-list()
  for (i in 1:36){
    selectedResultsList[[paste0("sim_",sim)]][[colnames(tempFile)[i]]]<-tempFile[,i]
  }
  tempFile<-read.csv(paste0("/home/ch19g17/Sim_for_zalpha_paper/selected/sim_",sim,"/simMidResults.csv"))
  selectedMidResultsList[[paste0("sim_",sim)]]<-list()
  for (i in 1:36){
    selectedMidResultsList[[paste0("sim_",sim)]][[colnames(tempFile)[i]]]<-tempFile[,i]
  }
}

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

# Set working directory for the graphs to go in
setwd("/home/ch19g17/Sim_for_zalpha_paper/Graphs/new")

#stop R ploting scientific notation
options(scipen=999)
#Adjust location of title so it doesn't overlap the legend
titleAdj<-rep(0.5,35)
titleAdj[c(22,24,25,28,30,31,33)]<-0.35
for (i in 1:35){
  binEndData<-data.frame(bin_pos=seq(0,990000,10000),mean=NA,sd=NA,n=NA,lowCI=NA,highCI=NA)
  binMidData<-data.frame(bin_pos=seq(0,990000,10000),mean=NA,sd=NA,n=NA,lowCI=NA,highCI=NA)
  binNeutData<-data.frame(bin_pos=seq(0,990000,10000),mean=NA,sd=NA,n=NA,lowCI=NA,highCI=NA)
  
  for(bin in 1:100){
    stat<-NULL
    for (sim in 1:100){
      stat <- c(stat,selectedResultsList[[paste0("sim_",sim)]][[i+1]][selectedResultsList[[paste0("sim_",sim)]]$position>=binEndData$bin_pos[bin] & selectedResultsList[[paste0("sim_",sim)]]$position<binEndData$bin_pos[bin]+10000])
    }
    stat<-stat[is.na(stat)==FALSE]
    binEndData$mean[bin]<-mean(stat)
    binEndData$sd[bin]<-sd(stat)
    binEndData$n[bin]<-length(stat)
    binEndData$lowCI[bin]<-binEndData$mean[bin]-qt(0.975,binEndData$n[bin]-1)*binEndData$sd[bin]/sqrt(binEndData$n[bin])
    binEndData$highCI[bin]<-binEndData$mean[bin]+qt(0.975,binEndData$n[bin]-1)*binEndData$sd[bin]/sqrt(binEndData$n[bin])
  }
  for(bin in 1:100){
    stat<-NULL
    for (sim in 1:100){
      stat <- c(stat,selectedMidResultsList[[paste0("sim_",sim)]][[i+1]][selectedMidResultsList[[paste0("sim_",sim)]]$position>=binMidData$bin_pos[bin] & selectedMidResultsList[[paste0("sim_",sim)]]$position<binMidData$bin_pos[bin]+10000])
    }
    stat<-stat[is.na(stat)==FALSE]
    binMidData$mean[bin]<-mean(stat)
    binMidData$sd[bin]<-sd(stat)
    binMidData$n[bin]<-length(stat)
    binMidData$lowCI[bin]<-binMidData$mean[bin]-qt(0.975,binMidData$n[bin]-1)*binMidData$sd[bin]/sqrt(binMidData$n[bin])
    binMidData$highCI[bin]<-binMidData$mean[bin]+qt(0.975,binMidData$n[bin]-1)*binMidData$sd[bin]/sqrt(binMidData$n[bin])
  }
  for(bin in 1:100){
    stat<-NULL
    for (sim in 1:100){
      stat <- c(stat,neutralResultsList[[paste0("sim_",sim)]][[i+1]][neutralResultsList[[paste0("sim_",sim)]]$position>=binNeutData$bin_pos[bin] & neutralResultsList[[paste0("sim_",sim)]]$position<binNeutData$bin_pos[bin]+10000])
    }
    stat<-stat[is.na(stat)==FALSE]
    binNeutData$mean[bin]<-mean(stat)
    binNeutData$sd[bin]<-sd(stat)
    binNeutData$n[bin]<-length(stat)
    binNeutData$lowCI[bin]<-binNeutData$mean[bin]-qt(0.975,binNeutData$n[bin]-1)*binNeutData$sd[bin]/sqrt(binNeutData$n[bin])
    binNeutData$highCI[bin]<-binNeutData$mean[bin]+qt(0.975,binNeutData$n[bin]-1)*binNeutData$sd[bin]/sqrt(binNeutData$n[bin])
  }
  if (names(selectedResultsList$sim_1)[i+1] == "L_plus_R" || names(selectedResultsList$sim_1)[i+1] == "LR"|| names(selectedResultsList$sim_1)[i+1] == "Zalpha_BetaCDF_over_Zbeta_BetaCDF"){
    ylims<-c(min(binNeutData$lowCI,binMidData$lowCI,binEndData$lowCI,na.rm=TRUE),max(binNeutData$highCI,binMidData$highCI,binEndData$highCI,na.rm=TRUE))
  } else {
    ylims<-c(min(0,binNeutData$lowCI,binMidData$lowCI,binEndData$lowCI,na.rm=TRUE),max(binNeutData$highCI,binMidData$highCI,binEndData$highCI,na.rm=TRUE))
  }
  
  pdf(paste0(names(selectedResultsList$sim_1)[i+1],".pdf"))
  par(mar=c(5,6,4,2)+0.1,xpd=TRUE,cex=1.25) #Increase left margin for formulae to fit
  plot(binEndData$bin_pos+5000,binEndData$mean,type="n",xaxp=c(0,1000000,4),ylim=ylims,xlab="Position (bp)",ylab=statsFormula[i])
  title(main=statsFormula[i],adj=titleAdj[i])
  polygonDF<-data.frame(x=c(binNeutData$bin_pos+5000,rev(binNeutData$bin_pos+5000)),y=c(binNeutData$lowCI,rev(binNeutData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="lightblue",border = NA,lty=1,density=50,angle=-30)
  polygonDF<-data.frame(x=c(binMidData$bin_pos+5000,rev(binMidData$bin_pos+5000)),y=c(binMidData$lowCI,rev(binMidData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="lightpink",border = NA,lty=1,density=50,angle=45)
  polygonDF<-data.frame(x=c(binEndData$bin_pos+5000,rev(binEndData$bin_pos+5000)),y=c(binEndData$lowCI,rev(binEndData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="gray",border = NA,lty=1,density=50,angle=-45)
  lines(binNeutData$bin_pos+5000,binNeutData$mean,col="darkblue")
  lines(binMidData$bin_pos+5000,binMidData$mean,col="red")
  lines(binEndData$bin_pos+5000,binEndData$mean,col="black")
  legend("topright",c("End","Mid-way","Neutral"),col=c("black","red","darkblue"),lty=1,inset=c(0,-0.213))
  dev.off()
  tiff(paste0(names(selectedResultsList$sim_1)[i+1],".tiff"),res=350,height=17,width=17,units="cm",compression="lzw")
  par(mar=c(5,6,4,2)+0.1,xpd=TRUE,cex=1.25) #Increase left margin for formulae to fit
  plot(binEndData$bin_pos+5000,binEndData$mean,type="n",xaxp=c(0,1000000,4),ylim=ylims,xlab="Position (bp)",ylab=statsFormula[i])
  title(main=statsFormula[i],adj=titleAdj[i])
  polygonDF<-data.frame(x=c(binNeutData$bin_pos+5000,rev(binNeutData$bin_pos+5000)),y=c(binNeutData$lowCI,rev(binNeutData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="lightblue",border = NA,lty=1,density=50,angle=-30)
  polygonDF<-data.frame(x=c(binMidData$bin_pos+5000,rev(binMidData$bin_pos+5000)),y=c(binMidData$lowCI,rev(binMidData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="lightpink",border = NA,lty=1,density=50,angle=45)
  polygonDF<-data.frame(x=c(binEndData$bin_pos+5000,rev(binEndData$bin_pos+5000)),y=c(binEndData$lowCI,rev(binEndData$highCI)))
  polygonDF<-polygonDF[is.na(polygonDF$y)==FALSE,]
  polygon(polygonDF$x,polygonDF$y,col="gray",border = NA,lty=1,density=50,angle=-45)
  lines(binNeutData$bin_pos+5000,binNeutData$mean,col="darkblue")
  lines(binMidData$bin_pos+5000,binMidData$mean,col="red")
  lines(binEndData$bin_pos+5000,binEndData$mean,col="black")
  legend("topright",c("End","Mid-way","Neutral"),col=c("black","red","darkblue"),lty=1,inset=c(0,-0.228))
  dev.off()
}
