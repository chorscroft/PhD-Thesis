## Wavelet example for thesis

#working directory
setwd("\\\\filestore.soton.ac.uk\\users\\ch19g17\\mydocuments\\Wavelet\\Example Wavelet Work")

##Load prerequisites
#source("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wavelet/Comparison/getWaveletFile.R")
require(wmtsa)
require(biwavelet)
require(fields)

## Create Fine Scale Changes dataset
#First time generating
#dfExample1<-data.frame(x=c(1:1024),y=rnorm(1024,0,1))
#write.csv(dfExample1,"Example1.csv",row.names=FALSE)

#Read from .csv from now on
dfExample1<-read.csv("Example1.csv")

##Create Wider Scale Changes dataset
#First time generating
#<-data.frame(x=c(1:1024),y=sin(c(1:1024)*pi/32)+rnorm(1024,0,1))
#write.csv(dfExample2,"Example2.csv",row.names=FALSE)

#Read from .csv from now on
dfExample2<-read.csv("Example2.csv")

##Create Wider Scale2 Changes dataset
#First time generating
#dfExample3<-data.frame(x=c(1:1024),y=sin(c(1:1024)*pi/32)+rnorm(1024,0,1))
#write.csv(dfExample3,"Example3.csv",row.names=FALSE)

#Read from .csv from now on
dfExample3<-read.csv("Example3.csv")

## Plot example signals

jpeg("ExampleSignals.jpeg",width=13,height=17,units="cm",quality=100,res=300)
par(mfcol=c(3,1))
plot(dfExample1$x,dfExample1$y,type="l",xlab="Kb",ylab="Signal",ylim=c(-5,5),xaxp=c(0,1024,8),main="Example 1: noise") 
plot(dfExample2$x,dfExample2$y,type="l",xlab="Kb",ylab="Signal",ylim=c(-5,5),xaxp=c(0,1024,8),main="Example 2: wave") 
plot(dfExample3$x,dfExample3$y,type="l",xlab="Kb",ylab="Signal",ylim=c(-5,5),xaxp=c(0,1024,8),main="Example 3: wave") 
dev.off()

## Get DWT

Example1_DWT<-wavDWT(dfExample1$y,wavelet="haar")
Example2_DWT<-wavDWT(dfExample2$y,wavelet="haar")
Example3_DWT<-wavDWT(dfExample3$y,wavelet="haar")


## Get approximation coefficients

getApproxCoeff <- function(x){
  approx<-list(s0=x)
  n<-length(x)
  level<-0
  while (n> 1) {
    n<-n/2
    level<-level+1
    s<-rep(NA,n)
    for (i in 1:n){
      s[i]<-(approx[[level]][2*i]+approx[[level]][2*i-1])/sqrt(2)
    }
    approx[[paste0("s",level)]]<-s
  }
  return(approx)
}

Example1_Approx<-getApproxCoeff(dfExample1$y)
Example2_Approx<-getApproxCoeff(dfExample2$y)
Example3_Approx<-getApproxCoeff(dfExample3$y)

## Plot the wavelet decomposition

plotBlocks<-function(y,xlim,ylab="",xaxt="n",yaxt="n",las=0,RHS=FALSE){
  blocks<-length(y)
  x<-seq(0,xlim[2],xlim[2]/blocks)
  plot(rep(0,length(y)),y,type="n",xlim=xlim,xaxp=c(0,xlim[2],8),ylab=ylab,xlab="",yaxt=yaxt,xaxt=xaxt,las=las)
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
plotWaveletDecomposition<-function(DWT,Approx,fileName,graphTitle){
  jpeg(fileName,width=14,height=21,units="cm",quality=100,res=300)
  par(mar=c(0.5,0,0,0.5),oma=c(4,9.5,4,4),xpd=NA)
  layout(matrix(c(1:22), 11, 2, byrow = FALSE))
  plotBlocks(Approx[[1]],xlim=c(1,1024),yaxt="s",las=2)
  mtext(paste0("Original\nSignal"),side = 2, line = 4.5,cex=0.75,las=1)
  title(main="Detail coefficients",line=1)
  mtext(graphTitle,outer=TRUE,line=2.5)
  for (i in 1:9){
    plotBlocks(DWT$data[[i]],xlim=c(1,1024),yaxt="s",las=2)
    mtext(paste0("Scale ",2^i),side = 2, line = 4.5,cex=0.75,las=1)
    mtext(paste0("d(",i,")"),side = 2, line = 3,cex=0.75)
  }
  #d10
  plot(0,0,type="n",xaxp=c(0,1024,8),xlim=c(1,1024),ylab="",xlab="Kb",xaxt="s",las=2)
  lines(c(0,1024),c(DWT$data[[10]],DWT$data[[10]]))
  mtext(paste0("Scale ",2^10),side = 2, line = 4.5,cex=0.75,las=1)
  mtext("d(10)",side = 2, line = 3,cex=0.75)
  for (i in 0:9){
    plotBlocks(Approx[[i+1]],xlim=c(1,1024),RHS=TRUE)
    if(i==0){title(main="Approximation coefficients",line=1)}
    mtext(paste0("a(",i,")"), side = 4, line = 3,cex =0.75)
  }
  plot(0,Approx[[11]],type="n",xaxp=c(0,1024,8),xlim=c(1,1024),ylab="",xlab="Kb",yaxt="n",las=2)
  lines(c(0,1024),c(Approx[[11]],Approx[[11]]))
  mtext("a(10)", side = 4, line = 3,cex =0.75)
  axis(side = 4,las=2)
  dev.off()
}
plotWaveletDecomposition(Example1_DWT,Example1_Approx,"Example1Decomposition.jpeg","DWT of Example 1")
plotWaveletDecomposition(Example2_DWT,Example2_Approx,"Example2Decomposition.jpeg","DWT of Example 2")
plotWaveletDecomposition(Example3_DWT,Example3_Approx,"Example3Decomposition.jpeg","DWT of Example 3")

## Get values for example reconstruction
Example2_DWT$data$s10
Example2_DWT$data$d10
Example2_Approx$s9
(0.0722-(-0.4326))/sqrt(2)
(0.0722+(-0.4326))/sqrt(2)

getPowerSpectrum<-function(dataDWT){
  powerSpectrum<-NULL
  for(i in 1:dataDWT$dictionary$n.levels){
    powerSpectrum[i]<-sum(dataDWT$data[[i]]^2,na.rm = TRUE)
  }
  powerSpectrum<-powerSpectrum/sum(powerSpectrum,na.rm=TRUE)
  return(powerSpectrum)
}
## Power Spectrum
Example1_PS<-getPowerSpectrum(Example1_DWT)
Example2_PS<-getPowerSpectrum(Example2_DWT)
Example3_PS<-getPowerSpectrum(Example3_DWT)

jpeg("PowerSpectrumExamples.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(Example1_PS,type="l",xlim=c(1,10),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(Example1_PS,Example2_PS,Example3_PS,na.rm=TRUE)),main="Proportion of Variance")
axis(1,1:10,labels=2^(1:10),las=2)
lines(Example2_PS,col="blue",lty=2)
lines(Example3_PS,col="red",lty=3)
legend("topright",legend=c("Example 1","Example 2","Example 3"),col=c("black","blue","red"),lty=c(1,2,3))
dev.off()

## Show what happens when rotate 15

dfExample2r<-dfExample2
dfExample2r$y<-c(dfExample2$y[-c(1:15)],dfExample2$y[1:15])
Example2r_DWT<-wavDWT(dfExample2r$y,wavelet="haar")
Example2r_Approx<-getApproxCoeff(dfExample2r$y)
plotWaveletDecomposition(Example2r_DWT,Example2r_Approx,"Example2r15Decomposition.jpeg","DWT of Example 2, rotated by 15")

Example2r_PS<-getPowerSpectrum(Example2r_DWT)
jpeg("PowerSpectrumExamplesr15.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(Example2r_PS,type="l",xlim=c(1,10),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(Example2r_PS,na.rm=TRUE)),main="Proportion of Variance",col="green",lty=2)
lines(Example2_PS,col="blue")
legend("topright",legend=c("Example 2","Example 2,\nrotated by 15",""),col=c("blue","green",NA),lty=c(1,2))
axis(1,1:10,labels=2^(1:10),las=2)
dev.off()

## Get MODWT

Example1_MODWT<-wavMODWT(dfExample1$y,wavelet="haar")
Example2_MODWT<-wavMODWT(dfExample2$y,wavelet="haar")
Example3_MODWT<-wavMODWT(dfExample3$y,wavelet="haar")

getApproxCoeffMO <- function(x){
  approx<-list(s0=x)
  n<-length(x)
  level<-0
  for (i in 1:10){
    level<-level+1
    s<-rep(NA,n)
    for (t in 1:n){
      index<-t-2^(i-1)
      if (index<=0){
        index <- 1024+index
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
Example1_MOApprox<-getApproxCoeffMO(dfExample1$y)
Example2_MOApprox<-getApproxCoeffMO(dfExample2$y)
Example3_MOApprox<-getApproxCoeffMO(dfExample3$y)


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
  layout(matrix(c(1:22), 11, 2, byrow = FALSE))
  plotBlocks(Approx[[1]],xlim=c(1,1024),yaxt="s",las=2)
  mtext(paste0("Original\nSignal"),side = 2, line = 4.5,cex=0.75,las=1)
  title(main="Detail coefficients",line=1)
  mtext(graphTitle,outer=TRUE,line=2.5)
  for (i in 1:10){
    plotBlocks(DWT$data[[i]],xlim=c(1,1024),yaxt="s",las=2,xaxt=ifelse(i==10,"s","n"),xlab=ifelse(i==10,"Kb",""))
    mtext(paste0("Scale ",2^i),side = 2, line = 4.5,cex=0.75,las=1)
    mtext(paste0("d(",i,")"),side = 2, line = 3,cex=0.75)
  }
  for (i in 0:10){
    plotBlocks(Approx[[i+1]],xlim=c(1,1024),RHS=TRUE,xaxt=ifelse(i==10,"s","n"),las=2,xlab=ifelse(i==10,"Kb",""))
    if(i==0){title(main="Approximation coefficients",line=1)}
    mtext(paste0("a(",i,")"), side = 4, line = 3,cex =0.75) #WAS line =0.5
  }
  dev.off()
}

plotWaveletMODecomposition(Example1_MODWT,Example1_MOApprox,"Example1MODecomposition.jpeg","MODWT of Example 1")
plotWaveletMODecomposition(Example2_MODWT,Example2_MOApprox,"Example2MODecomposition.jpeg","MODWT of Example 2")
plotWaveletMODecomposition(Example3_MODWT,Example3_MOApprox,"Example3MODecomposition.jpeg","MODWT of Example 3")

## Power spectrum of MODWT
Example1_MOPS<-getPowerSpectrum(Example1_MODWT)
Example2_MOPS<-getPowerSpectrum(Example2_MODWT)
Example3_MOPS<-getPowerSpectrum(Example3_MODWT)

jpeg("PowerSpectrumMOExamples.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(Example1_MOPS,type="l",xlim=c(1,10),xlab="Scale Kb",ylab="Proportion of Variance",xaxt="n",ylim=c(0,max(Example1_MOPS,Example2_MOPS,Example3_MOPS,na.rm=TRUE)),main="Proportion of Variance with MODWT")
axis(1,1:10,labels=2^(1:10),las=2)
lines(Example2_MOPS,col="blue",lty=2)
lines(Example3_MOPS,col="red",lty=3)
legend("topright",legend=c("Example 1","Example 2","Example 3"),col=c("black","blue","red"),lty=c(1,2,3))
dev.off()

## Wavelet correlation
#function
getWaveletCorrelation10<-function(x,y){ ##x and y are DWT objects
  est<-NULL
  pval<-NULL
  for (i in 1:9){  ##need more than 2 observations
    rankCor<-cor.test(x$data[[i]],y$data[[i]],method="kendall")
    est[i]<-rankCor$estimate
    pval[i]<-rankCor$p.value
  }
  return(list(est=est,pval=pval))
}

#get correlations
Example12_cor<-getWaveletCorrelation10(Example1_MODWT,Example2_MODWT)
Example13_cor<-getWaveletCorrelation10(Example1_MODWT,Example3_MODWT)
Example23_cor<-getWaveletCorrelation10(Example2_MODWT,Example3_MODWT)

jpeg("CorrelationExamples.jpeg",width=15,height=15,units="cm",quality=100,res=300)
plot(Example12_cor$est[1:8],type="l",main="Correlations",ylim=c(-1,1),col="black",xaxt="n",xlab="Scale Kb",ylab="Kendall's tau")
axis(1,at=c(1:8),labels=c(2^(1:8)),las=2)
points(which(Example12_cor$pval[1:8]<0.01),Example12_cor$est[which(Example12_cor$pval[1:8]<0.01)],col="black")
lines(Example13_cor$est[1:8],col="blue",lty=2)
points(which(Example13_cor$pval[1:8]<0.01),Example13_cor$est[which(Example13_cor$pval[1:8]<0.01)],col="blue")
lines(Example23_cor$est[1:8],col="red",lty=3)
points(which(Example23_cor$pval[1:8]<0.01),Example23_cor$est[which(Example23_cor$pval[1:8]<0.01)],col="red")
legend("bottomleft",legend=c("Examples 1 and 2","Examples 1 and 3","Examples 2 and 3"),col=c("black","blue","red"),lty=c(1,2,3))
legend("bottomright",legend="p-value < 0.01",pch=1)
dev.off()




# ----------------------------------------------------------------------------------------

###CWT

plotLegendCWT<-function(x){
  x$power <- x$power.corr
  yvals <- log2(x$period)
  zvals <- log2(abs(x$power/x$sigma2))
  zlim <- range(c(-1, 1) * max(zvals))
  zvals[zvals < zlim[1]] <- zlim[1]
  #locs <- pretty(range(zlim), n = 7)
  locs<-c(-9,-6,-3,0,3,6,9)
  leg.lab<-2^locs
  leg.lab[leg.lab<1] <- paste0("1/",(1/leg.lab[leg.lab<1]))
  fill.cols <- c("#00007F", "blue", "#007FFF", "cyan", 
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  col.pal <- colorRampPalette(fill.cols)
  fill.colors <- col.pal(64)
  image.plot(x$t, yvals, t(zvals), 
             zlim = zlim, 
             ylim = rev(range(yvals)), 
             col = fill.colors, 
             horizontal = FALSE, 
             legend.only = TRUE, 
             axis.args = list(at = locs,labels = leg.lab),
             xpd = NA)
}

Example1_CWT<-wt(dfExample1,mother="morlet",param=6)
Example2_CWT<-wt(dfExample2,mother="morlet",param=6)
Example3_CWT<-wt(dfExample3,mother="morlet",param=6)

bmp("Example1CWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(Example1_CWT, plot.cb = FALSE, plot.phase = FALSE,main="Example 1",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(Example1_CWT)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(dfExample1$x,dfExample1$y,type="l",ylab="Signal",xlab="Kb",xaxs="i",xaxp=c(0,1024,8),las=1)
dev.off()
bmp("Example2CWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(Example2_CWT, plot.cb = FALSE, plot.phase = FALSE,main="Example 2",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(Example2_CWT)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(dfExample2$x,dfExample2$y,type="l",ylab="Signal",xlab="Kb",xaxs="i",xaxp=c(0,1024,8),las=1)
dev.off()
bmp("Example3CWT.bmp",res=300,height=17,width=23,units='cm')
par(mfrow=c(2,1), oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(Example3_CWT, plot.cb = FALSE, plot.phase = FALSE,main="Example 3",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(Example3_CWT)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(dfExample3$x,dfExample3$y,type="l",ylab="Signal",xlab="Kb",xaxs="i",xaxp=c(0,1024,8),las=1)
dev.off()

# Combine example 1 and 2
bmp("CWT_12combine.bmp",res=300,height=17,width=46,units='cm')
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
par(oma = c(0, 0, 0, 1), mar = c(0, 4, 4, 5) + 0.1)
plot(Example1_CWT, plot.cb = FALSE, plot.phase = FALSE,main="Example 1",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(Example1_CWT)
plot(Example2_CWT, plot.cb = FALSE, plot.phase = FALSE,main="Example 2",xlab="",xaxt="n",lwd.sig=2,las=1)
plotLegendCWT(Example2_CWT)
par(mar = c(5, 4, 0, 5) + 0.1)
plot(dfExample1$x,dfExample1$y,type="l",ylab="Signal",xlab="Kb",xaxs="i",xaxp=c(0,1024,8),las=1)
plot(dfExample2$x,dfExample2$y,type="l",ylab="Signal",xlab="Kb",xaxs="i",xaxp=c(0,1024,8),las=1)
dev.off()



## Wavelet Coherence
Example12_coh<-wtc(dfExample1,dfExample2,mother="morlet",param=6,nrands = 1000)
Example13_coh<-wtc(dfExample1,dfExample3,mother="morlet",param=6,nrands = 1000)
Example23_coh<-wtc(dfExample2,dfExample3,mother="morlet",param=6,nrands = 1000)

jpeg("Example12Coh.jpeg",width=15,height=15,units="cm",quality=100,res=300)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(Example12_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Examples 1 and 2")
dev.off()
jpeg("Example13Coh.jpeg",width=15,height=15,units="cm",quality=100,res=300)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(Example13_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Examples 1 and 3")
dev.off()
jpeg("Example23Coh.jpeg",width=15,height=15,units="cm",quality=100,res=300)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(Example23_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Wavelet Coherence: Examples 2 and 3")
dev.off()

## plot all on one graph

bmp("Coherence_combine.bmp",res=300,height=30,width=30,units='cm')
par(mfrow=c(2,2),oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(Example12_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Examples 1 and 2")
plot(Example13_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Examples 1 and 3")
plot(Example23_coh, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2,  ylab = "Scale", xlab = "KB", 
     plot.cb = TRUE, main = "Examples 2 and 3")
dev.off()

## plot Morlet 6


f <- function(x) (pi^(-1/4))*exp(-(0+1i)*6*x)*exp(-(x^2)/2)
x <- seq(-5, 5, by=0.00001)
png("Morlet6.png",width=15,height=15,units="cm",res=300)
plot(x,Re(f(x)),type="l",ylim=c(-0.8,0.8),ylab=expression(psi(s)),xlab="s")
lines(x,Im(f(x)),lty=3,col="blue")
legend("topright",legend=c("Real","Imaginary"),lty=c(1,3),col=c("black","blue"))
dev.off()



