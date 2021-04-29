

## Draws the Haar wavelet where:
##       psi(s) =  1  for 0   <  s <= 1/2
##              = -1  for 1/2 <= s <  1
##              =  0  otherwise

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wavelet")
tiff("HaarWavelet.tiff",res=300,height=8,width=17,units='cm',compression="lzw")

par(mfrow=c(1,4),mar=c(5,0,4,0),oma=c(1,5,1,1))
## original
## s=[-1,0] psi[s] = [1,-1]
plot((-2):2,(-2):2,type="n",main="a) Haar Wavelet",xlab="s",xlim=c(-2,2))
title(ylab=expression(paste(psi,"(s)")),outer=TRUE)
abline(h=0,col="grey")
abline(v=0,col="grey")
box()
segments(x0=-3,y0=0,x1=0) #across
segments(x0=0,y0=0,y1=1) #up
segments(x0=0,y0=1,x1=0.5) #across
segments(x0=0.5,y0=1,y1=-1) #down
segments(x0=0.5,y0=-1,x1=1) #across
segments(x0=1,y0=-1,y1=0) #up
segments(x0=1,y0=0,x1=3) #across

##shift
## s=[0,1] psi[s] = [1,-1]
plot((-2):2,(-2):2,type="n",main="b) Shifted",xlab="s",ylab="", yaxt="n",xlim=c(-2,2))
abline(h=0,col="grey")
abline(v=0,col="grey")
box()
segments(x0=-3,y0=0,x1=-1) #across
segments(x0=-1,y0=0,y1=1) #up
segments(x0=-1,y0=1,x1=-0.5) #across
segments(x0=-0.5,y0=1,y1=-1) #down
segments(x0=-0.5,y0=-1,x1=0) #across
segments(x0=0,y0=-1,y1=0) #up
segments(x0=0,y0=0,x1=3) #across

##stretch
## s=[-1,1] psi[s] = [-1/sqrt(2),1/sqrt(2)]
plot((-2):2,(-2):2,type="n",main="c) Stretched",xlab="s",ylab="", yaxt="n",xlim=c(-2,2))
abline(h=0,col="grey")
abline(v=0,col="grey")
box()
segments(x0=3,y0=0,x1=1) #across
segments(x0=1,y0=0,y1=(-1/sqrt(2))) #down
segments(x0=1,y0=(-1/sqrt(2)),x1=0) #across
segments(x0=0,y0=(-1/sqrt(2)),y1=(1/sqrt(2))) #up
segments(x0=0,y0=(1/sqrt(2)),x1=-1) #across
segments(x0=-1,y0=(1/sqrt(2)),y1=0) #down
segments(x0=-1,y0=0,x1=-3) #across


##stretch shift and flip
## s=[-1,1] psi[s] = [-1/sqrt(2),1/sqrt(2)]
plot((-2):2,(-2):2,type="n",main="d) Flipped",xlab="s",ylab="", yaxt="n",xlim=c(-2,2))
abline(h=0,col="grey")
abline(v=0,col="grey")
box()
segments(x0=-3,y0=0,x1=-1) #across
segments(x0=-1,y0=0,y1=(-1/sqrt(2))) #down
segments(x0=-1,y0=(-1/sqrt(2)),x1=0) #across
segments(x0=0,y0=(-1/sqrt(2)),y1=(1/sqrt(2))) #up
segments(x0=0,y0=(1/sqrt(2)),x1=1) #across
segments(x0=1,y0=(1/sqrt(2)),y1=0) #down
segments(x0=1,y0=0,x1=3) #across

par(mfrow=c(1,1))

dev.off()