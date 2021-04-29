##Wavelet comparison of European and African data

source("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/convertLDhatToCMperMB2.R")
source("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wavelet/Comparison/getWaveletFile.R")
require(intervals)
require(wmtsa)
require(biwavelet)
require(fields)
#.libPaths("/home/ch19g17/R/R-packages")

### Parameters
### \\filestore.soton.ac.uk\users\ch19g17\mydocuments\LDHat\Map Length\COMBO.xlsx for workings
trimmedStart<- 20000428
trimmedEnd<- 51218377
map_length_cM<-60.67105


### Get CEU data ###

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22")
resC <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosC <- 20000.428
endPosC <- 51219.641
binC<-getWaveletFile(resC,startPosC,endPosC,trimmedStart,trimmedEnd,map_length_cM)

### Wellderly ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Wellderly")
resW <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosW <- 16051.295
endPosW <- 51218.377       ##16051.295+35167.083-0.001
binW<-getWaveletFile(resW,startPosW,endPosW,trimmedStart,trimmedEnd,map_length_cM)

### Baganda ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Baganda")
resB <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosB <- 16051.347
endPosB <- 51238.130       ##16051.347+35186.784-0.001
binB<-getWaveletFile(resB,startPosB,endPosB,trimmedStart,trimmedEnd,map_length_cM)

### Zulu ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Zulu")
resZ <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosZ <- 16051.347
endPosZ <- 51238.130       ##16051.347+35186.784-0.001
binZ<-getWaveletFile(resZ,startPosZ,endPosZ,trimmedStart,trimmedEnd,map_length_cM)

### Set working directory for output
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wavelet/Comparison New")



### Put all together ###
genomeData<-data.frame(chrom="chr22",nbases=max(binC[nrow(binC),2],binW[nrow(binW),2],binB[nrow(binB),2],binZ[nrow(binZ),2]))
mergedBin<-binGenome_per_chr("chr22",1000,rm.gaps = FALSE,genomeData)
mergedBin<-data.frame(mergedBin@.Data)
mergedBin<-merge(mergedBin,binC[,c(1,3,4)],by.x="X1",by.y="V1",all=TRUE)
colnames(mergedBin)[c(3,4)]<-c("C_rate","C_lograte")
mergedBin<-merge(mergedBin,binW[,c(1,3,4)],by.x="X1",by.y="V1",all=TRUE)
colnames(mergedBin)[c(5,6)]<-c("W_rate","W_lograte")
mergedBin<-merge(mergedBin,binB[,c(1,3,4)],by.x="X1",by.y="V1",all=TRUE)
colnames(mergedBin)[c(7,8)]<-c("B_rate","B_lograte")
mergedBin<-merge(mergedBin,binZ[,c(1,3,4)],by.x="X1",by.y="V1",all=TRUE)
colnames(mergedBin)[c(9,10)]<-c("Z_rate","Z_lograte")
mergedBin<-na.omit(mergedBin)

##just use section where all datasets have complete data
mergedBin <- mergedBin[mergedBin$X1>=25074000 & mergedBin$X1<50415000,]

region<-4454:20837

# Find C DWT
C_data<-getWaveletFormat(mergedBin[region,c(1,2,3)])
C_DWT<-wavDWT(C_data[,3],wavelet="haar")
#plot(C_DWT)

## Get power spectrum
C_PS<-getPowerSpectrum(C_DWT)
#plot(C_PS,type="l")


## Find W DWT
W_data<-getWaveletFormat(mergedBin[region,c(1,2,5)])
W_DWT<-wavDWT(W_data[,3],wavelet="haar")
#plot(W_DWT)

## Get power spectrum
W_PS<-getPowerSpectrum(W_DWT)
#plot(W_PS,type="l")


## Find B DWT
B_data<-getWaveletFormat(mergedBin[region,c(1,2,7)])
B_DWT<-wavDWT(B_data[,3],wavelet="haar")
#plot(B_DWT)

## Get power spectrum
B_PS<-getPowerSpectrum(B_DWT)
#plot(B_PS,type="l")


## Find Z DWT
Z_data<-getWaveletFormat(mergedBin[region,c(1,2,9)])
Z_DWT<-wavDWT(Z_data[,3],wavelet="haar")
#plot(Z_DWT)

## Get power spectrum
Z_PS<-getPowerSpectrum(Z_DWT)
#plot(Z_PS,type="l")


### Plot all power spectrums on same graph
jpeg("PowerSpectrum.jpg",res=300,height=14,width=14,units='cm',quality=100)
plot(C_PS,type="l",xlim=c(1,14),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(C_PS,W_PS,B_PS,Z_PS,na.rm=TRUE)),main="Proportion of Variance Comparison")
axis(1,1:14,labels=2^(1:14),las=2)
lines(W_PS,col="blue")
lines(B_PS,col="orange")
lines(Z_PS,col="red")
legend("topright",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
dev.off()

## log-10 transformed data

## Find C DWT
C_data_log<-getWaveletFormat(mergedBin[region,c(1,2,4)])
C_DWT_log<-wavDWT(C_data_log[,3],wavelet="haar")
#plot(C_DWT_log)

## Get power spectrum
C_PS_log<-getPowerSpectrum(C_DWT_log)
#plot(C_PS_log,type="l")


## Find W DWT
W_data_log<-getWaveletFormat(mergedBin[region,c(1,2,6)])
W_DWT_log<-wavDWT(W_data_log[,3],wavelet="haar")
#plot(W_DWT_log)

## Get power spectrum
W_PS_log<-getPowerSpectrum(W_DWT_log)
#plot(W_PS_log,type="l")


## Find B DWT
B_data_log<-getWaveletFormat(mergedBin[region,c(1,2,8)])
B_DWT_log<-wavDWT(B_data_log[,3],wavelet="haar")
#plot(B_DWT_log)

## Get power spectrum
B_PS_log<-getPowerSpectrum(B_DWT_log)
#plot(B_PS_log,type="l")


## Find Z DWT
Z_data_log<-getWaveletFormat(mergedBin[region,c(1,2,10)])
Z_DWT_log<-wavDWT(Z_data_log[,3],wavelet="haar")
#plot(Z_DWT_log)

## Get power spectrum
Z_PS_log<-getPowerSpectrum(Z_DWT_log)
#plot(Z_PS_log,type="l")


### Plot all power spectrums on same graph
bmp("PowerSpectrumLog.bmp",res=300,height=14,width=14,units='cm')#,compression="lzw")
plot(C_PS_log,type="l",xlim=c(1,14),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(C_PS_log,W_PS_log,B_PS_log,Z_PS_log,na.rm=TRUE)),main="Proportion of Variance Comparison (logged)")
axis(1,1:14,labels=2^(1:14),las=2)
lines(W_PS_log,col="blue")
lines(B_PS_log,col="orange")
lines(Z_PS_log,col="red")
legend("topright",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
dev.off()

## Get the rank correlations for the DWT
CW_cor_log<-getWaveletCorrelation(C_DWT_log,W_DWT_log)
CB_cor_log<-getWaveletCorrelation(C_DWT_log,B_DWT_log)
CZ_cor_log<-getWaveletCorrelation(C_DWT_log,Z_DWT_log)
WB_cor_log<-getWaveletCorrelation(W_DWT_log,B_DWT_log)
WZ_cor_log<-getWaveletCorrelation(W_DWT_log,Z_DWT_log)
BZ_cor_log<-getWaveletCorrelation(B_DWT_log,Z_DWT_log)

## plot all correlations on same graph
bmp("Correlations.bmp",res=300,height=17,width=17,units='cm')
plot(CW_cor_log$est[1:12],type="l",main="Correlations",ylim=c(0,1),col="black",xaxt="n",xlab="Scale Kb",ylab="Kendall's tau")
axis(1,at=c(1:12),labels=c(2^(1:12)),las=2)
points(which(CW_cor_log$pval[1:12]<0.01),CW_cor_log$est[which(CW_cor_log$pval[1:12]<0.01)],col="black")
lines(CB_cor_log$est[1:12],col="blue")
points(which(CB_cor_log$pval[1:12]<0.01),CB_cor_log$est[which(CB_cor_log$pval[1:12]<0.01)],col="blue")
lines(CZ_cor_log$est[1:12],col="orange")
points(which(CZ_cor_log$pval[1:12]<0.01),CZ_cor_log$est[which(CZ_cor_log$pval[1:12]<0.01)],col="orange")
lines(WB_cor_log$est[1:12],col="purple")
points(which(WB_cor_log$pval[1:12]<0.01),WB_cor_log$est[which(WB_cor_log$pval[1:12]<0.01)],col="purple")
lines(WZ_cor_log$est[1:12],col="red")
points(which(WZ_cor_log$pval[1:12]<0.01),WZ_cor_log$est[which(WZ_cor_log$pval[1:12]<0.01)],col="red")
lines(BZ_cor_log$est[1:12],col="orange")
points(which(BZ_cor_log$pval[1:12]<0.01),BZ_cor_log$est[which(BZ_cor_log$pval[1:12]<0.01)],col="orange")
legend("bottomright",legend=c("CEU & Wellderly","CEU & Baganda","CEU & Zulu","Wellderly & Baganda","Wellderly & Zulu","Baganda & Zulu"),col=c("black","blue","orange","purple","red","orange"),lty=1)
legend("bottomleft",legend="p<0.01",col="black",pch=1)
dev.off()




###########################################################
## MODWT

#Get MODWT
C_MODWT<-wavMODWT(C_data[,3],wavelet="haar")
W_MODWT<-wavMODWT(W_data[,3],wavelet="haar")
B_MODWT<-wavMODWT(B_data[,3],wavelet="haar")
Z_MODWT<-wavMODWT(Z_data[,3],wavelet="haar")

## Get power spectrum
C_MOPS<-getPowerSpectrum(C_MODWT)
W_MOPS<-getPowerSpectrum(W_MODWT)
B_MOPS<-getPowerSpectrum(B_MODWT)
Z_MOPS<-getPowerSpectrum(Z_MODWT)

#Plot power spectrum
bmp("powerspectrum.bmp",res=300,height=14,width=14,units='cm')
plot(C_MOPS,type="l",xlim=c(1,14),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(C_MOPS,W_MOPS,B_MOPS,Z_MOPS,na.rm=TRUE)),main="Proportion of Variance Comparison MODWT")
axis(1,1:14,labels=2^(1:14),las=2)
lines(W_MOPS,col="blue")
lines(B_MOPS,col="orange")
lines(Z_MOPS,col="red")
legend("topright",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
dev.off()
## log-10 transformed data MODWT

## Find MODWT for log10 data
C_MODWT_log<-wavMODWT(C_data_log[,3],wavelet="haar")
W_MODWT_log<-wavMODWT(W_data_log[,3],wavelet="haar")
B_MODWT_log<-wavMODWT(B_data_log[,3],wavelet="haar")
Z_MODWT_log<-wavMODWT(Z_data_log[,3],wavelet="haar")

## Get power spectrum
C_MOPS_log<-getPowerSpectrum(C_MODWT_log)
W_MOPS_log<-getPowerSpectrum(W_MODWT_log)
B_MOPS_log<-getPowerSpectrum(B_MODWT_log)
Z_MOPS_log<-getPowerSpectrum(Z_MODWT_log)

#Plot power spectrum
bmp("powerspectrumlog.bmp",res=300,height=14,width=14,units='cm')
plot(C_MOPS_log,type="l",xlim=c(1,14),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(C_MOPS_log,W_MOPS_log,B_MOPS_log,Z_MOPS_log,na.rm=TRUE)),main="Proportion of Variance Comparison (logged) MODWT")
axis(1,1:14,labels=2^(1:14),las=2)
lines(W_MOPS_log,col="blue")
lines(B_MOPS_log,col="orange")
lines(Z_MOPS_log,col="red")
legend("topright",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
dev.off()

## Get the rank correlations for the DWT
CW_MOcor_log<-getWaveletCorrelation(C_MODWT_log,W_MODWT_log)
CB_MOcor_log<-getWaveletCorrelation(C_MODWT_log,B_MODWT_log)
CZ_MOcor_log<-getWaveletCorrelation(C_MODWT_log,Z_MODWT_log)
WB_MOcor_log<-getWaveletCorrelation(W_MODWT_log,B_MODWT_log)
WZ_MOcor_log<-getWaveletCorrelation(W_MODWT_log,Z_MODWT_log)
BZ_MOcor_log<-getWaveletCorrelation(B_MODWT_log,Z_MODWT_log)

## plot all correlations on same graph
bmp("correlations.bmp",res=300,height=14,width=14,units='cm')
plot(CW_MOcor_log$est[1:12],type="l",main="Correlations MODWT",ylim=c(0,1),col="black",xaxt="n",xlab="Scale Kb",ylab="Kendall's tau")
axis(1,at=c(1:12),labels=c(2^(1:12)),las=2)
points(which(CW_MOcor_log$pval[1:12]<0.01),CW_MOcor_log$est[which(CW_MOcor_log$pval[1:12]<0.01)],col="black")
lines(CB_MOcor_log$est[1:12],col="blue")
points(which(CB_MOcor_log$pval[1:12]<0.01),CB_MOcor_log$est[which(CB_MOcor_log$pval[1:12]<0.01)],col="blue")
lines(CZ_MOcor_log$est[1:12],col="orange")
points(which(CZ_MOcor_log$pval[1:12]<0.01),CZ_MOcor_log$est[which(CZ_MOcor_log$pval[1:12]<0.01)],col="orange")
lines(WB_MOcor_log$est[1:12],col="purple")
points(which(WB_MOcor_log$pval[1:12]<0.01),WB_MOcor_log$est[which(WB_MOcor_log$pval[1:12]<0.01)],col="purple")
lines(WZ_MOcor_log$est[1:12],col="pink")
points(which(WZ_MOcor_log$pval[1:12]<0.01),WZ_MOcor_log$est[which(WZ_MOcor_log$pval[1:12]<0.01)],col="pink")
lines(BZ_MOcor_log$est[1:12],col="red")
points(which(BZ_MOcor_log$pval[1:12]<0.01),BZ_MOcor_log$est[which(BZ_MOcor_log$pval[1:12]<0.01)],col="red")
legend("bottomright",legend=c("CEU & Wellderly","CEU & Baganda","CEU & Zulu","Wellderly & Baganda","Wellderly & Zulu","Baganda & Zulu"),col=c("black","blue","orange","purple","pink","red"),lty=1)
legend("bottomleft",legend="p<0.01",col="black",pch=1)
dev.off()













# ##CWT
# 
# #CEU
# wtC<-wt(cbind(C_data$X1/1000,C_data$C_rate),mother="morlet",param=6)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
# plot(wtC, plot.cb = TRUE, plot.phase = FALSE,xlab="KB",main="CEU")
# plot(wtC, plot.cb = TRUE, plot.phase = FALSE,xlab="KB",main="CEU",xlim=c(37000,38000))
# 
# #Wellderly
# wtW<-wt(cbind(W_data$X1/1000,W_data$W_rate),mother="morlet",param=6)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
# plot(wtW, plot.cb = TRUE, plot.phase = FALSE,xlab="KB",main="Wellderly")
# 
# #Baganda
# wtB<-wt(cbind(B_data$X1/1000,B_data$B_rate),mother="morlet",param=6)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
# plot(wtB, plot.cb = TRUE, plot.phase = FALSE,xlab="KB",main="Baganda")
# 
# #Zulu
# wtZ<-wt(cbind(Z_data$X1/1000,Z_data$Z_rate),mother="morlet",param=6)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
# plot(wtZ, plot.cb = TRUE, plot.phase = FALSE,xlab="KB",main="Zulu")

##CWT LOG

#CEU
wtClog<-wt(cbind(C_data_log$X1/1000,C_data_log$C_lograte),mother="morlet",param=6)
bmp("CEUCWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtClog, plot.cb = FALSE, plot.phase = FALSE,main="CEU",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtClog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(C_data_log$X1/1000,C_data_log$C_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(C_data_log$X1[1]/1000,C_data_log$X1[length(C_data_log$X1)]/1000),xaxs="i",las=1)
dev.off()

#Wellderly
wtWlog<-wt(cbind(W_data_log$X1/1000,W_data_log$W_lograte),mother="morlet",param=6)
bmp("WellderlyCWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtWlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Wellderly",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtWlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(W_data_log$X1/1000,W_data_log$W_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(W_data_log$X1[1]/1000,W_data_log$X1[length(W_data_log$X1)]/1000),xaxs="i",las=1)
dev.off()

#Baganda
wtBlog<-wt(cbind(B_data_log$X1/1000,B_data_log$B_lograte),mother="morlet",param=6)
bmp("BagandaCWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtBlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Baganda",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtBlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(B_data_log$X1/1000,B_data_log$B_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(B_data_log$X1[1]/1000,B_data_log$X1[length(B_data_log$X1)]/1000),xaxs="i",las=1)
dev.off()

#Zulu
wtZlog<-wt(cbind(Z_data_log$X1/1000,Z_data_log$Z_lograte),mother="morlet",param=6)
bmp("ZuluCWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtZlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Zulu",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtZlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(Z_data_log$X1/1000,Z_data_log$Z_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(Z_data_log$X1[1]/1000,Z_data_log$X1[length(Z_data_log$X1)]/1000),xaxs="i",las=1)
dev.off()

## plot all on same graph
bmp("combineCWT.bmp",res=300,height=34,width=46,units='cm')
layout(matrix(c(1:8),4,2,byrow=TRUE))
par(oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtClog, plot.cb = FALSE, plot.phase = FALSE,main="CEU",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtClog)
plot(wtWlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Wellderly",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtWlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(C_data_log$X1/1000,C_data_log$C_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(C_data_log$X1[1]/1000,C_data_log$X1[length(C_data_log$X1)]/1000),xaxs="i",las=1)
plot(W_data_log$X1/1000,W_data_log$W_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(W_data_log$X1[1]/1000,W_data_log$X1[length(W_data_log$X1)]/1000),xaxs="i",las=1)
par(mar = c(0, 4, 4, 5) + 0.1)
plot(wtBlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Baganda",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtBlog)
plot(wtZlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Zulu",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(wtZlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(B_data_log$X1/1000,B_data_log$B_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(B_data_log$X1[1]/1000,B_data_log$X1[length(B_data_log$X1)]/1000),xaxs="i",las=1)
plot(Z_data_log$X1/1000,Z_data_log$Z_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xlim=c(Z_data_log$X1[1]/1000,Z_data_log$X1[length(Z_data_log$X1)]/1000),xaxs="i",las=1)
dev.off()




##window
zoomregion<-c(37000,38000)
bmp("zoomCWT.bmp",res=300,height=34,width=46,units='cm')
layout(matrix(c(1:8),4,2,byrow=TRUE))
par(oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(wtClog, plot.cb = FALSE, plot.phase = FALSE,main="CEU",xlab="",xaxt="n",lwd.sig=2,las=1,xlim=zoomregion)
plotLegendCWT(wtClog)
plot(wtWlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Wellderly",xaxt="n",lwd.sig=2,las=1,xlim=zoomregion)
plotLegendCWT(wtWlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(C_data_log$X1/1000,C_data_log$C_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xaxs="i",las=1,xlim=zoomregion)
plot(W_data_log$X1/1000,W_data_log$W_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xaxs="i",las=1,xlim=zoomregion)
par(mar = c(0, 4, 4, 5) + 0.1)
plot(wtBlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Baganda",xaxt="n",lwd.sig=2,las=1,xlim=zoomregion)
plotLegendCWT(wtBlog)
plot(wtZlog, plot.cb = FALSE, plot.phase = FALSE,xlab="",main="Zulu",xaxt="n",lwd.sig=2,las=1,xlim=zoomregion)
plotLegendCWT(wtZlog)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(B_data_log$X1/1000,B_data_log$B_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xaxs="i",las=1,xlim=zoomregion)
plot(Z_data_log$X1/1000,Z_data_log$Z_lograte,type="l",ylab="log(cM per Mb)",xlab="Kb",xaxs="i",las=1,xlim=zoomregion)
dev.off()



# ## wavelet coherance
# wavCohCW<-wtc(cbind(C_data$X1/1000,C_data$C_rate),cbind(W_data$X1/1000,W_data$W_rate),mother="morlet",param=6,nrands = 1000)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
# plot(wavCohCW, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
#      lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
#      plot.cb = TRUE, main = "Wavelet Coherence: C vs W")
# 
# wavCohCB<-wtc(cbind(C_data$X1/1000,C_data$C_rate),cbind(B_data$X1/1000,B_data$B_rate),mother="morlet",param=6,nrands = 1000)
# par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
# plot(wavCohCB, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
#      lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
#      plot.cb = TRUE, main = "Wavelet Coherence: C vs B")


## THIS SECTION RUNS ON IRIDIS
## log wavelet coherance
wavCohCWlog<-wtc(cbind(C_data_log$X1/1000,C_data_log$C_lograte),cbind(W_data_log$X1/1000,W_data_log$W_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohCW.bmp",res=300,height=17,width=23,units='cm')
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohCWlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Wellderly")
dev.off()

wavCohCBlog<-wtc(cbind(C_data_log$X1/1000,C_data_log$C_lograte),cbind(B_data_log$X1/1000,B_data_log$B_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohCB.bmp",res=300,height=17,width=23,units='cm',compression="lzw")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohCBlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Baganda")
dev.off()




wavCohCZlog<-wtc(cbind(C_data_log$X1/1000,C_data_log$C_lograte),cbind(Z_data_log$X1/1000,Z_data_log$Z_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohCZ.bmp",res=300,height=17,width=23,units='cm',compression="lzw")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohCZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Zulu")
dev.off()

wavCohWBlog<-wtc(cbind(W_data_log$X1/1000,W_data_log$W_lograte),cbind(B_data_log$X1/1000,B_data_log$B_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohWB.bmp",res=300,height=17,width=23,units='cm',compression="lzw")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohWBlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Wellderly & Baganda")
dev.off()

wavCohWZlog<-wtc(cbind(W_data_log$X1/1000,W_data_log$W_lograte),cbind(Z_data_log$X1/1000,Z_data_log$Z_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohWZ.bmp",res=300,height=17,width=23,units='cm',compression="lzw")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohWZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Wellderly & Zulu")
dev.off()

wavCohBZlog<-wtc(cbind(B_data_log$X1/1000,B_data_log$B_lograte),cbind(Z_data_log$X1/1000,Z_data_log$Z_lograte),mother="morlet",param=6,nrands = 1000)
bmp("wavCohBZ.bmp",res=300,height=17,width=23,units='cm',compression="lzw")
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohBZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Baganda & Zulu")
dev.off()

## combined graph 
bmp("wavCohcombined.bmp",res=300,height=51,width=46,units='cm')
layout(matrix(c(1:6),3,2,byrow = TRUE))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wavCohCWlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Wellderly")
plot(wavCohCBlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Baganda")
plot(wavCohCZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: CEU & Zulu")
plot(wavCohWBlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Wellderly & Baganda")
plot(wavCohWZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Wellderly & Zulu")
plot(wavCohBZlog, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Baganda & Zulu")
dev.off()







####
##correlation between bins in region

getKendallCorrelations<-function(mergedBin,col1,col2,region){
  returnData<-list()
  temp<-cor.test(mergedBin[region,col1],mergedBin[region,col2],method="kendall")
  returnData$cor<-temp$est
  returnData$pvalue<-temp$p.value
  return(returnData)
}
# KcorCW<-getKendallCorrelations(mergedBin,3,5,region)
# KcorCB<-getKendallCorrelations(mergedBin,3,7,region)
# KcorWB<-getKendallCorrelations(mergedBin,5,7,region)
# KcorCZ<-getKendallCorrelations(mergedBin,3,9,region)
# KcorWZ<-getKendallCorrelations(mergedBin,5,9,region)
# KcorBZ<-getKendallCorrelations(mergedBin,7,9,region)

KcorCW<-getKendallCorrelations(mergedBin,4,6,region)
KcorCB<-getKendallCorrelations(mergedBin,4,8,region)
KcorWB<-getKendallCorrelations(mergedBin,6,8,region)
KcorCZ<-getKendallCorrelations(mergedBin,4,10,region)
KcorWZ<-getKendallCorrelations(mergedBin,6,10,region)
KcorBZ<-getKendallCorrelations(mergedBin,8,10,region)

##plot correlations
require(corrplot)
Mr<-matrix(0,nrow=4,ncol=4,dimnames=list(c("CEU","Wellderly","Baganda","Zulu"),c("CEU","Wellderly","Baganda","Zulu")))
Mr[1,2]<-KcorCW$cor
Mr[1,3]<-KcorCB$cor
Mr[1,4]<-KcorCZ$cor
Mr[2,3]<-KcorWB$cor
Mr[2,4]<-KcorWZ$cor
Mr[3,4]<-KcorBZ$cor
Mr[2,1]<-KcorCW$cor
Mr[3,1]<-KcorCB$cor
Mr[4,1]<-KcorCZ$cor
Mr[3,2]<-KcorWB$cor
Mr[4,2]<-KcorWZ$cor
Mr[4,3]<-KcorBZ$cor
Mp<-matrix(0,nrow=4,ncol=4,dimnames=list(c("CEU","Wellderly","Baganda","Zulu"),c("CEU","Wellderly","Baganda","Zulu")))
Mp[1,2]<-KcorCW$pvalue
Mp[1,3]<-KcorCB$pvalue
Mp[1,4]<-KcorCZ$pvalue
Mp[2,3]<-KcorWB$pvalue
Mp[2,4]<-KcorWZ$pvalue
Mp[3,4]<-KcorBZ$pvalue
Mp[2,1]<-KcorCW$pvalue
Mp[3,1]<-KcorCB$pvalue
Mp[4,1]<-KcorCZ$pvalue
Mp[3,2]<-KcorWB$pvalue
Mp[4,2]<-KcorWZ$pvalue
Mp[4,3]<-KcorBZ$pvalue
par(mfrow=c(1,1))
corrplot(Mr,method="color",type="upper",addCoef.col = "black",diag=FALSE,tl.col="black",p.mat=Mp,sig.level=0.01,cl.lim=c(0,1))



### Investigate low coherence area around 40Mb
invest<-c(40000,45000)
plot(C_data_log$X1/1000,C_data_log$C_lograte,xlim=invest,type="l")
lines(W_data_log$X1/1000,W_data_log$W_lograte,xlim=invest,col="red")

rollav<-data.frame(x=C_data_log$X1/1000,c=C_data_log$C_lograte,w=W_data_log$W_lograte)
rollav<-rollav[rollav$x>=invest[1] & rollav$x<=invest[2],]
rollav$croll <- NA
rollav$wroll<-NA
for (i in 40500:44500){
  rollav$croll[rollav$x==i]<-mean(rollav$c[rollav$x>=(i-500) & rollav$x<=(i+500)])
  rollav$wroll[rollav$x==i]<-mean(rollav$w[rollav$x>=(i-500) & rollav$x<=(i+500)])
}

plot(rollav$x,rollav$croll,type="l",ylim=c(-2,0))
lines(rollav$x,rollav$wroll,col="red")

rollav$cchange<-NA
rollav$wchange<-NA
for (i in 1:5000){
  rollav$cchange[i]<-rollav$c[i+1]-rollav$c[i]
  rollav$wchange[i]<-rollav$w[i+1]-rollav$w[i]
  
}
plot(rollav$x,rollav$cchange,type="l")
lines(rollav$x,rollav$wchange,col="red")
plot(rollav$cchange,rollav$wchange)


plot(1:16,C_DWT_log$data$d10,type="l")
lines(1:16,W_DWT_log$data$d10,col="red")
plot(C_DWT_log$data$d10,W_DWT_log$data$d10)

write.csv(wavCohCWlog$rsq,"wavcoh.csv")

#----------------------------------------------------------
## Get cumulative cM 

bmp("cM.bmp",res=300,height=17,width=17,units='cm')

ccmC<-getCumulativeCM(resC,startPosC,endPosC,trimmedStart,trimmedEnd,map_length_cM)
plot(ccmC$pos/1000000,ccmC$cM,type="n",xlab="Position (Mb) on Chromosome 22", ylab="cM", xaxs="i",yaxs="i",xlim=c(20,51.218377),main="Cumulative centimorgan maps")
rect(29.527,0,45.911,65,col="azure2",lty=0)
box()
lines(ccmC$pos/1000000,ccmC$cM,col="black")
ccmW<-getCumulativeCM(resW,startPosW,endPosW,trimmedStart,trimmedEnd,map_length_cM)
lines(ccmW$pos/1000000,ccmW$cM,col="blue")
ccmB<-getCumulativeCM(resB,startPosB,endPosB,trimmedStart,trimmedEnd,map_length_cM)
lines(ccmB$pos/1000000,ccmB$cM,col="orange")
ccmZ<-getCumulativeCM(resZ,startPosZ,endPosZ,trimmedStart,trimmedEnd,map_length_cM)
lines(ccmZ$pos/1000000,ccmZ$cM,col="red")
legend("topleft",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)

dev.off()


##### Cold spot around 32.5 Mb in European data and not African

plot(ccmC$pos[ccmC$pos>=32300000 & ccmC$pos<=32550000]/1000000,ccmC$cM[ccmC$pos>=32300000 & ccmC$pos<=32550000],type="n",xlab="Position (Mb) on Chromosome 22", ylab="cM", xaxs="i",yaxs="i",main="Cumulative centiMorgan maps",ylim=c(20,24))
points(ccmC$pos[ccmC$pos>=32300000 & ccmC$pos<=32550000]/1000000,ccmC$cM[ccmC$pos>=32300000 & ccmC$pos<=32550000],col="black")
points(ccmZ$pos[ccmZ$pos>=32300000 & ccmZ$pos<=32550000]/1000000,ccmZ$cM[ccmZ$pos>=32300000 & ccmZ$pos<=32550000],col="red")

plot(C_data$X1[C_data$X1>=32300000 & C_data$X1<=32550000],C_data$C_rate[C_data$X1>=32300000 & C_data$X1<=32550000])
plot(Z_data$X1[Z_data$X1>=32300000 & Z_data$X1<=32550000],Z_data$Z_rate[Z_data$X1>=32300000 & Z_data$X1<=32550000])

## Plot for European cold spot
bmp("Europeancoldspot.bmp",res=300,height=17,width=17,units='cm')
plot(C_data_log$X1[C_data_log$X1/1000>=32300 & C_data_log$X1/1000<=32550]/1000,C_data_log$C_lograte[C_data_log$X1/1000>=32300 & C_data_log$X1/1000<=32550],type="l",col="black",ylab="log(cM per Mb)",xlab="Kb",main="European cold spot")
lines(W_data_log$X1[W_data_log$X1/1000>=32300 & W_data_log$X1/1000<=32550]/1000,W_data_log$W_lograte[W_data_log$X1/1000>=32300 & W_data_log$X1/1000<=32550],col="blue")
lines(B_data_log$X1[B_data_log$X1/1000>=32300 & B_data_log$X1/1000<=32550]/1000,B_data_log$B_lograte[B_data_log$X1/1000>=32300 & B_data_log$X1/1000<=32550],col="orange")
lines(Z_data_log$X1[Z_data_log$X1/1000>=32300 & Z_data_log$X1/1000<=32550]/1000,Z_data_log$Z_lograte[Z_data_log$X1/1000>=32300 & Z_data_log$X1/1000<=32550],col="red")
legend("bottomleft",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
dev.off()

# zooming in on regions
xval<-c(31500,32000)
plot(C_data_log$X1[C_data_log$X1/1000>=xval[1] & C_data_log$X1/1000<=xval[2]]/1000,C_data_log$C_lograte[C_data_log$X1/1000>=xval[1] & C_data_log$X1/1000<=xval[2]],type="l",col="black",ylab="log(cM per Mb)",xlab="Kb",main="European cold spot")
lines(W_data_log$X1[W_data_log$X1/1000>=xval[1] & W_data_log$X1/1000<=xval[2]]/1000,W_data_log$W_lograte[W_data_log$X1/1000>=xval[1] & W_data_log$X1/1000<=xval[2]],col="blue")
lines(B_data_log$X1[B_data_log$X1/1000>=xval[1] & B_data_log$X1/1000<=xval[2]]/1000,B_data_log$B_lograte[B_data_log$X1/1000>=xval[1] & B_data_log$X1/1000<=xval[2]],col="orange")
lines(Z_data_log$X1[Z_data_log$X1/1000>=xval[1] & Z_data_log$X1/1000<=xval[2]]/1000,Z_data_log$Z_lograte[Z_data_log$X1/1000>=xval[1] & Z_data_log$X1/1000<=xval[2]],col="red")
legend("bottomleft",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)





plot((resC$Loci[(resC$Loci+startPosC-0.001)>=32300 & (resC$Loci+startPosC-0.001)<=32550]+startPosC-0.001)/1000,resC$Mean_rho[(resC$Loci+startPosC-0.001)>=32300 & (resC$Loci+startPosC-0.001)<=32550])
plot((resZ$Loci[(resZ$Loci+startPosZ-0.001)>=32300 & (resZ$Loci+startPosZ-0.001)<=32550]+startPosZ-0.001)/1000,resZ$Mean_rho[(resZ$Loci+startPosZ-0.001)>=32300 & (resZ$Loci+startPosZ-0.001)<=32550])


##########################################
# Graphs for other regions just to check
region<-4454:20837
region<-1:16384
region<-8907:25290

C_data_log<-getWaveletFormat(mergedBin[region,c(1,2,4)])
W_data_log<-getWaveletFormat(mergedBin[region,c(1,2,6)])
B_data_log<-getWaveletFormat(mergedBin[region,c(1,2,8)])
Z_data_log<-getWaveletFormat(mergedBin[region,c(1,2,10)])
C_MODWT_log<-wavMODWT(C_data_log[,3],wavelet="haar")
W_MODWT_log<-wavMODWT(W_data_log[,3],wavelet="haar")
B_MODWT_log<-wavMODWT(B_data_log[,3],wavelet="haar")
Z_MODWT_log<-wavMODWT(Z_data_log[,3],wavelet="haar")

## Get power spectrum
C_MOPS_log<-getPowerSpectrum(C_MODWT_log)
W_MOPS_log<-getPowerSpectrum(W_MODWT_log)
B_MOPS_log<-getPowerSpectrum(B_MODWT_log)
Z_MOPS_log<-getPowerSpectrum(Z_MODWT_log)

#Plot power spectrum
#bmp("powerspectrumlog.bmp",res=300,height=14,width=14,units='cm')
plot(C_MOPS_log,type="l",xlim=c(1,14),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(C_MOPS_log,W_MOPS_log,B_MOPS_log,Z_MOPS_log,na.rm=TRUE)),main="Proportion of Variance Comparison (logged) MODWT")
axis(1,1:14,labels=2^(1:14),las=2)
lines(W_MOPS_log,col="blue")
lines(B_MOPS_log,col="orange")
lines(Z_MOPS_log,col="red")
#legend("topright",legend=c("CEU","Wellderly","Baganda","Zulu"),col=c("black","blue","orange","red"),lty=1)
#dev.off()

## Get the rank correlations for the DWT
CW_MOcor_log<-getWaveletCorrelation(C_MODWT_log,W_MODWT_log)
CB_MOcor_log<-getWaveletCorrelation(C_MODWT_log,B_MODWT_log)
CZ_MOcor_log<-getWaveletCorrelation(C_MODWT_log,Z_MODWT_log)
WB_MOcor_log<-getWaveletCorrelation(W_MODWT_log,B_MODWT_log)
WZ_MOcor_log<-getWaveletCorrelation(W_MODWT_log,Z_MODWT_log)
BZ_MOcor_log<-getWaveletCorrelation(B_MODWT_log,Z_MODWT_log)

## plot all correlations on same graph
#bmp("correlations.bmp",res=300,height=14,width=14,units='cm')
plot(CW_MOcor_log$est[1:12],type="l",main="Correlations MODWT",ylim=c(0,1),col="black",xaxt="n",xlab="Scale Kb",ylab="Kendall's tau")
axis(1,at=c(1:12),labels=c(2^(1:12)),las=2)
points(which(CW_MOcor_log$pval[1:12]<0.01),CW_MOcor_log$est[which(CW_MOcor_log$pval[1:12]<0.01)],col="black")
lines(CB_MOcor_log$est[1:12],col="blue")
points(which(CB_MOcor_log$pval[1:12]<0.01),CB_MOcor_log$est[which(CB_MOcor_log$pval[1:12]<0.01)],col="blue")
lines(CZ_MOcor_log$est[1:12],col="orange")
points(which(CZ_MOcor_log$pval[1:12]<0.01),CZ_MOcor_log$est[which(CZ_MOcor_log$pval[1:12]<0.01)],col="orange")
lines(WB_MOcor_log$est[1:12],col="purple")
points(which(WB_MOcor_log$pval[1:12]<0.01),WB_MOcor_log$est[which(WB_MOcor_log$pval[1:12]<0.01)],col="purple")
lines(WZ_MOcor_log$est[1:12],col="pink")
points(which(WZ_MOcor_log$pval[1:12]<0.01),WZ_MOcor_log$est[which(WZ_MOcor_log$pval[1:12]<0.01)],col="pink")
lines(BZ_MOcor_log$est[1:12],col="red")
points(which(BZ_MOcor_log$pval[1:12]<0.01),BZ_MOcor_log$est[which(BZ_MOcor_log$pval[1:12]<0.01)],col="red")
#legend("bottomright",legend=c("CEU & Wellderly","CEU & Baganda","CEU & Zulu","Wellderly & Baganda","Wellderly & Zulu","Baganda & Zulu"),col=c("black","blue","orange","purple","pink","red"),lty=1)
#legend("bottomleft",legend="p<0.01",col="black",pch=1)
#dev.off()


#####
plotWaveletMODecomposition(C_MODWT_log,getApproxCoeffMO(C_data_log$C_lograte),"MODWT_C.jpeg","MODWT of CEU")


getApproxCoeffMO <- function(x){
  approx<-list(s0=x)
  n<-length(x)
  level<-0
  for (i in 1:14){
    level<-level+1
    s<-rep(NA,n)
    for (t in 1:n){
      index<-t-2^(i-1)
      if (index<=0){
        index <- 16384+index
      }
      s[t]<-(approx[[level]][t]+approx[[level]][index])/(2)  #rescaled by multiplying by another 1/sqrt(2)
    }
    approx[[paste0("s",level)]]<-s
  }
  #check final level for rounding errors
  for (i in 2:length(approx[[level+1]])){
    if (isTRUE(all.equal(approx[[level+1]][i-1],approx[[level+1]][i]))){
      approx[[level+1]][i]<-approx[[level+1]][i-1]
    }
  }
  return(approx)
}

plotBlocks<-function(y,xlim,ylab="",xaxt="n",yaxt="n",las=0,RHS=FALSE,xlab=""){
  blocks<-length(y)
  x<-seq(0,xlim[2],xlim[2]/blocks)
  plot(rep(0,length(y)+2),c(y,0.05,-0.05),type="n",xlim=xlim,xaxp=c(0,xlim[2],8),ylab=ylab,xlab=xlab,yaxt=yaxt,xaxt=xaxt,las=las)
  if (RHS==TRUE){
    axis(side = 4,las=2)
  }
  for(i in 1:blocks){
    lines(x=c(x[i],x[i+1]),y=c(y[i],y[i]))
    if(i < blocks){
      lines(x=c(x[i+1],x[i+1]),y=c(y[i],y[i+1]))
    }
  }
}
plotWaveletMODecomposition<-function(DWT,Approx,fileName,graphTitle){
  jpeg(fileName,width=14,height=21,units="cm",quality=100,res=300)
  par(mar=c(0.5,0,0,0.5),oma=c(4,9.5,4,4),xpd=NA)
  layout(matrix(c(1:30), 15, 2, byrow = FALSE))
  plotBlocks(Approx[[1]],xlim=c(1,16384),yaxt="s",las=2)
  mtext(paste0("Original\nSignal"),side = 2, line = 4.5,cex=0.75,las=1)
  title(main="Detail coefficients",line=1)
  mtext(graphTitle,outer=TRUE,line=2.5)
  for (i in 1:14){
    plotBlocks(DWT$data[[i]],xlim=c(1,16384),yaxt="s",las=2,xaxt=ifelse(i==14,"s","n"),xlab=ifelse(i==14,"Kb",""))
    mtext(paste0("Scale ",2^i),side = 2, line = 4.5,cex=0.75,las=1)
    mtext(paste0("d(",i,")"),side = 2, line = 3,cex=0.75)
  }
  for (i in 0:14){
    plotBlocks(Approx[[i+1]],xlim=c(1,16384),RHS=TRUE,xaxt=ifelse(i==14,"s","n"),las=2,xlab=ifelse(i==14,"Kb",""))
    if(i==0){title(main="Approximation coefficients",line=1)}
    mtext(paste0("a(",i,")"), side = 4, line = 3,cex =0.75) #WAS line =0.5
  }
  dev.off()
}


