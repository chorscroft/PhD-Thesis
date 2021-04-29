###nSL flags:
# 0 = not nsl
# 1 = positive nsl only
# 2 = negative nsl only
# 3 = absoulte nsl
# 4 = nsl raw


getAveragePlot <- function(simFolder,midsection,prefix,nSLFlag){
  ### average graph of stat
  
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  nSims <- length(list.files(simFolder))
  sim <- 1
  mergeDF <- read.table(paste0(simFolder,prefix,"_",sim,".txt"),header=isTRUE(nSLFlag>0))
  if (nSLFlag ==1){
    mergeDF <- mergeDF[mergeDF[,2]>0,] #pos
  }
  if (nSLFlag ==2){
    mergeDF <- mergeDF[mergeDF[,2]<0,] #neg
  }
  if (nSLFlag ==3){
    mergeDF[mergeDF[,2]<0,2] <- -mergeDF[mergeDF[,2]<0,2] #abs
  }  
  colnames(mergeDF) = c("SNPpos",prefix)
  for (sim in 2:nSims){
    simDF <- read.table(paste0(simFolder,prefix,"_",sim,".txt"),header=isTRUE(nSLFlag>0))
    if (nSLFlag == 1){
      simDF <- simDF[simDF[,2]>0,] #pos
    }
    if (nSLFlag == 2){
      simDF <- simDF[simDF[,2]<0,] #neg
    }
    if (nSLFlag ==3){
      simDF[simDF[,2]<0,2] <- -simDF[simDF[,2]<0,2] #abs
    }
    colnames(simDF)<-c("SNPpos",prefix)
    mergeDF <- merge(mergeDF,simDF,by=1,all=TRUE)
  }
  averageStat <- cbind(mergeDF$SNPpos,apply(mergeDF[,c(-1)],1,sum,na.rm=TRUE),apply(mergeDF[,c(-1)],1,count),apply(mergeDF[,c(-1)],1,sumofsquare))
  averageStat <- cbind(sapply(averageStat[,1],round1000),averageStat)
  binStat <- seq(0,max(averageStat[,1]),1000)
  binStat <- cbind(binStat,rep(0,length(binStat)),rep(0,length(binStat)),rep(0,length(binStat)))
  for (i in 1:nrow(averageStat)){
    binStat[averageStat[i,1]/1000+1,2] <- binStat[averageStat[i,1]/1000+1,2]+averageStat[i,3]
    binStat[averageStat[i,1]/1000+1,3] <- binStat[averageStat[i,1]/1000+1,3]+averageStat[i,4]
    binStat[averageStat[i,1]/1000+1,4] <- binStat[averageStat[i,1]/1000+1,4]+averageStat[i,5]
  }
  binStat <- cbind(binStat,rep(NA,nrow(binStat)),rep(NA,nrow(binStat)),rep(NA,nrow(binStat)),rep(NA,nrow(binStat)))
  for (i in 1:nrow(binStat)){
    if (binStat[i,3]>0){
      binStat[i,5] <- binStat[i,2]/binStat[i,3]
      binStat[i,6] <- sqrt((binStat[i,4]-((binStat[i,2]^2)/binStat[i,3]))/(binStat[i,3]-1))
      binStat[i,7] <- binStat[i,5]+binStat[i,6]
      binStat[i,8] <- binStat[i,5]-binStat[i,6]
    }
  }
  binStat <- binStat[binStat[,1]>=midsection[1] & binStat[,1]<=midsection[2],]
  binStat[,1]<-binStat[,1]/100000
  return(binStat[,c(1,5,7,8)])
}
round1000 <- function(x){
  return(round(x/1000,0)*1000)
}
count <- function(x){
  return(sum(!is.na(x)))
}
sumofsquare <- function(x){
  return(sum(x^2,na.rm = TRUE))
}

tiff(file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/180412 Version after Peer Review/Version 3/figure3_v3tiff.tiff",res=300,height=17,width=17,units='cm',compression="lzw")

par(mar=c(4,4,4,4))
layout(matrix(c(1,2,3,4,5,5),3,2,byrow=TRUE),heights=c(3,3,1))

midsection <- c(150001,350000)
zalphaNeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaNeutral/"
h12NeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12Neutral/"
zalphaSelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaSelectedP/"
h12SelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12SelectedP/"
tSelNeutralTFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutputT/"
tSelSelectedPTFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputPT/"

#zalpha
zalphaNeutralPlot <-getAveragePlot(zalphaNeutralFolder,midsection,"zalpha",0)
zalphaSelectedPPlot <-getAveragePlot(zalphaSelectedPFolder,midsection,"zalpha",0)

#combined
plot(zalphaSelectedPPlot[,c(1,2)],xlab=expression(paste('10'^5," bp")),main=expression(paste("Z",''[alpha])),xlim=c(1.5,3.5),ylim=c(0.010,0.024),type="n",xaxs="i",yaxs="i",las=1)
arrows(zalphaSelectedPPlot[,1],zalphaSelectedPPlot[,3],zalphaSelectedPPlot[,1],zalphaSelectedPPlot[,4],length=0.05,angle=90,code=3,col="pink")
arrows(zalphaNeutralPlot[,1],zalphaNeutralPlot[,3],zalphaNeutralPlot[,1],zalphaNeutralPlot[,4],length=0.05,angle=90,code=3,col="gray")
points(zalphaSelectedPPlot[,c(1,2)],col="red",pch=20)
points(zalphaNeutralPlot[,c(1,2)],col="black",pch=20)
box()

#H12
H12NeutralPlot <-getAveragePlot(h12NeutralFolder,midsection,"h12",0)
H12SelectedPPlot <-getAveragePlot(h12SelectedPFolder,midsection,"h12",0)

#combined
plot(H12SelectedPlot[,c(1,2)],xlab=expression(paste('10'^5," bp")),main="H12",xlim=c(1.5,3.5),ylim=c(0,0.25),type="n",xaxs="i",yaxs="i",las=1,font.main=1)
arrows(H12SelectedPPlot[,1],H12SelectedPPlot[,3],H12SelectedPPlot[,1],H12SelectedPPlot[,4],length=0.05,angle=90,code=3,col="pink")
arrows(H12NeutralPlot[,1],H12NeutralPlot[,3],H12NeutralPlot[,1],H12NeutralPlot[,4],length=0.05,angle=90,code=3,col="gray")
points(H12SelectedPPlot[,c(1,2)],col="red",pch=20)
points(H12NeutralPlot[,c(1,2)],col="black",pch=20)
box()

###Standardise nSL
nslNeutralStandFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslNeutralStand/"
nslSelectedPStandFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedPStand/"

nslNeutralStandPlotabs <-getAveragePlot(nslNeutralStandFolder,midsection,"nslNeutral",3)
nslSelectedPStandPlotabs <-getAveragePlot(nslSelectedPStandFolder,midsection,"nslSelect",3)

#Plot abs
plot(nslSelectedPStandPlotabs[,c(1,2)],xlab=expression(paste('10'^5," bp")),main=expression(paste("|nS",''[L],"|")),xlim=c(1.5,3.5),ylim=c(0,3),type="n",xaxs="i",yaxs="i",las=1)
arrows(nslSelectedPStandPlotabs[,1],nslSelectedPStandPlotabs[,3],nslSelectedPStandPlotabs[,1],nslSelectedPStandPlotabs[,4],length=0.05,angle=90,code=3,col="pink")
arrows(nslNeutralStandPlotabs[,1],nslNeutralStandPlotabs[,3],nslNeutralStandPlotabs[,1],nslNeutralStandPlotabs[,4],length=0.05,angle=90,code=3,col="gray")
points(nslSelectedPStandPlotabs[,c(1,2)],col="red",pch=20)
points(nslNeutralStandPlotabs[,c(1,2)],col="black",pch=20)
box()

#tSel
tSelNeutralTPlot <-getAveragePlot(tSelNeutralTFolder,midsection,"tSelOutput",FALSE)
tSelSelectedPTPlot <-getAveragePlot(tSelSelectedPTFolder,midsection,"tSelOutput",FALSE)
#combined
plot(zalphaSelectedPPlot[,c(1,2)],xlab=expression(paste('10'^5," bp")),main="TSel",xlim=c(1.5,3.5),ylim=c(0,25),type="n",xaxs="i",yaxs="i",las=1,font.main=1)
arrows(tSelSelectedPTPlot[,1],tSelSelectedPTPlot[,3],tSelSelectedPTPlot[,1],tSelSelectedPTPlot[,4],length=0.05,angle=90,code=3,col="pink")
arrows(tSelNeutralTPlot[,1],tSelNeutralTPlot[,3],tSelNeutralTPlot[,1],tSelNeutralTPlot[,4],length=0.05,angle=90,code=3,col="gray")
points(tSelSelectedPTPlot[,c(1,2)],col="red",pch=20)
points(tSelNeutralTPlot[,c(1,2)],col="black",pch=20)
box()

par(mar=c(1,1,1,1))
plot(1,1, type="n", axes=FALSE, xlab="", ylab="")
legend("center",c("Neutral","Neutral SD","Selected","Selected SD"),col=c("black","gray","red","pink"),lty=1,cex=1,horiz=TRUE)

dev.off()




par(mfrow=c(1,1))





















nslNeutralStandPlot <-getAveragePlot(nslNeutralStandFolder,midsection,"nslNeutral",4)
nslSelectedPStandPlot <-getAveragePlot(nslSelectedPStandFolder,midsection,"nslSelect",4)
nslNeutralStandPlotpos <-getAveragePlot(nslNeutralStandFolder,midsection,"nslNeutral",1)
nslSelectedPStandPlotpos <-getAveragePlot(nslSelectedPStandFolder,midsection,"nslSelect",1)
nslNeutralStandPlotneg <-getAveragePlot(nslNeutralStandFolder,midsection,"nslNeutral",2)
nslSelectedPStandPlotneg <-getAveragePlot(nslSelectedPStandFolder,midsection,"nslSelect",2)


#Plot raw
plot(nslSelectedPStandPlot[,c(1,2)],xlab="bp",main=expression(paste("nS",''[L]," raw")),xlim=c(150000,350000),ylim=c(min(nslSelectedPStandPlot[,4])-0.5,max(nslSelectedPStandPlot[,3])+0.5),type="n",xaxs="i",yaxs="i",las=1)
arrows(nslSelectedPStandPlot[,1],nslSelectedPStandPlot[,3],nslSelectedPStandPlot[,1],nslSelectedPStandPlot[,4],length=0.05,angle=90,code=3,col="gray")
arrows(nslNeutralStandPlot[,1],nslNeutralStandPlot[,3],nslNeutralStandPlot[,1],nslNeutralStandPlot[,4],length=0.05,angle=90,code=3,col="pink")
points(nslSelectedPStandPlot[,c(1,2)],col="black",pch=20)
points(nslNeutralStandPlot[,c(1,2)],col="red",pch=20)
box()
legend("topright",c("Neutral","Selected","SD"),col=c("red","black","gray"),lty=1,cex=0.5)

#Plot pos/neg
plot(nslSelectedPStandPlotpos[,c(1,2)],xlab="bp",main=expression(paste("nS",''[L]," pos/neg")),xlim=c(150000,350000),ylim=c(min(nslSelectedPStandPlotneg[,4]-0.5),max(nslSelectedPStandPlotpos[,3])+0.5),type="n",xaxs="i",yaxs="i",las=1)
arrows(nslSelectedPStandPlotpos[,1],nslSelectedPStandPlotpos[,3],nslSelectedPStandPlotpos[,1],nslSelectedPStandPlotpos[,4],length=0.05,angle=90,code=3,col="gray")
arrows(nslNeutralStandPlotpos[,1],nslNeutralStandPlotpos[,3],nslNeutralStandPlotpos[,1],nslNeutralStandPlotpos[,4],length=0.05,angle=90,code=3,col="pink")
arrows(nslSelectedPStandPlotneg[,1],nslSelectedPStandPlotneg[,3],nslSelectedPStandPlotneg[,1],nslSelectedPStandPlotneg[,4],length=0.05,angle=90,code=3,col="gray")
arrows(nslNeutralStandPlotneg[,1],nslNeutralStandPlotneg[,3],nslNeutralStandPlotneg[,1],nslNeutralStandPlotneg[,4],length=0.05,angle=90,code=3,col="pink")
points(nslSelectedPStandPlotpos[,c(1,2)],col="black",pch=20)
points(nslNeutralStandPlotpos[,c(1,2)],col="red",pch=20)
points(nslSelectedPStandPlotneg[,c(1,2)],col="black",pch=20)
points(nslNeutralStandPlotneg[,c(1,2)],col="red",pch=20)
box()
legend("topright",c("Neutral","Selected","SD"),col=c("red","black","gray"),lty=1,cex=0.5)






plot(nslNeutralStandPlot,xlab="BP",ylab="nsl",main="nsl Neutral")
#plot(nslSelectedStandPlot,xlab="BP",ylab="nsl",main="nsl: Current 90% freq")
plot(nslSelectedPStandPlot,xlab="BP",ylab="nsl",main="nsl: Current 50% freq")
plot(nslNeutralStandPlotneg,xlab="BP",ylab="nsl",main="nsl Neutral")
#plot(nslSelectedStandPlotneg,xlab="BP",ylab="nsl",main="nsl: Current 90% freq")
plot(nslSelectedPStandPlotneg,xlab="BP",ylab="nsl",main="nsl: Current 50% freq")

plot(nslSelectedStandPlot,xlab="BP",ylab="nsl",main="nsl",type="l",ylim=c(-2,max(nslSelectedStandPlot[,2])))
lines(nslNeutralStandPlot,col="red")
lines(nslSelectedPStandPlot,col="green")
lines(nslSelectedStandPlotneg,col="black")
lines(nslNeutralStandPlotneg,col="red")
lines(nslSelectedPStandPlotneg,col="green")
legend("topright",c("Neutral","90% freq","50% freq"),col=c("red","black","green"),lty=1,cex=0.75)
