

convertLDhatTocMperMb2<-function(loci,rate,start,end,trimmedStart,trimmedEnd,bp=0.001,kb=TRUE,map_length_cM){
  
  loci<-loci+start-bp
  loci<-c(loci,end)
  if (kb==TRUE){
    loci<-loci*1000
    rate<-rate/1000
  }
  
  #trim front and back
  rate<-rate[loci>=trimmedStart]
  loci<-loci[loci>=trimmedStart]
  rate<-rate[loci<=trimmedEnd]
  rate<-rate[-length(rate)]
  loci<-loci[loci<=trimmedEnd]
  
  #work out cumulative rho
  cumulative_rho<-0
  for (i in 2:length(loci)){
    cumulative_rho[i]<-rate[i-1]*(loci[i]-loci[i-1])/1000000+cumulative_rho[i-1]
  }
  conversion_factor<-map_length_cM/cumulative_rho[length(cumulative_rho)]
  cumulative_cM<-cumulative_rho*conversion_factor
  
  cMMb<-rep(0,length(loci))
  for (i in 1:(length(loci)-1)){
    cMMb[i]<-(cumulative_cM[i+1]-cumulative_cM[i])*1000000/(loci[i+1]-loci[i])
  }
  cmReturn<-data.frame(loci,cumulative_cM,cMMb)
  
  return(cmReturn)
  
  
 
}

#example:

# bob<-convertLDhatTocMperMb(resC$Loci[-1],resC$Mean_rho[-1],TRUE,10000)
# loci<-resC$Loci[-1]
# rate<-resC$Mean_rho[-1]
# start <- 20000.428
# end <- 51219.641
# trimmedStart<-20000428
# trimmedEnd<-51218377
# bp=0.001
# map_length_cM <- 60.67105
# 
# 


##testing files

# resCloci<-(resC$Loci[-1]+startPosC-0.001)*1000
# resWloci<-(resW$Loci[-1]+startPosW-0.001)*1000
# resBloci<-(resB$Loci[-1]+startPosB-0.001)*1000
# resZloci<-(resZ$Loci[-1]+startPosZ-0.001)*1000
# 
# write.csv(resCloci,"resC.csv",row.names=FALSE)
# write.csv(resWloci,"resW.csv",row.names=FALSE)
# write.csv(resBloci,"resB.csv",row.names=FALSE)
# write.csv(resZloci,"resZ.csv",row.names=FALSE)

