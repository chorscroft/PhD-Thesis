library(grDevices)

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/African Maps Paper/Simulations/Output")

# GDmap<-read.csv("GDmap.csv")
# GSmap<-read.csv("GSmap.csv")
# 
# 
# plot(GDmap$kb,GDmap$LDU,type="l",ylim=c(0,100),xlab="Kb",ylab="LDU")
# lines(GSmap$kb,GSmap$LDU,col="blue")
# legend("topleft",legend=c("Dense","Sparse"),col=c("black","blue"),lty=1)
# 
# GDmap$kb[nrow(GDmap)]

### for all 16 maps

maps<-vector("list",16)
for (i in 5:20){
  mapname<-paste0("G",i,"map.csv")
  tempmap<-read.csv(mapname)
  maps[[i-4]]<-tempmap
  maps[[i-4]]$EBT<-maps[[i-4]]$LDU/0.1
}
map0<-read.csv("G0map.csv")
map0$EBT<-map0$LDU/0.1

colourPalette<-heat.colors(16)
bmp("kbLDU.bmp",res=300,height=14,width=14,units='cm')
plot(maps[[1]]$kb,maps[[1]]$LDU,type="l",ylim=c(0,100),xlab="Kb",ylab="LDU",col=colourPalette[1])
for(i in 2:16){
  lines(maps[[i]]$kb,maps[[i]]$LDU,col=colourPalette[i])
}
#lines(map0$kb,map0$LDU,col="black")
legend("topleft",title="Genes/Mb",legend=c(5,6,7,8,"...",17,18,19,20),col=c(colourPalette[1],colourPalette[2],colourPalette[3],colourPalette[4],"white",colourPalette[13],colourPalette[14],colourPalette[15],colourPalette[16]),lty=1)
dev.off()

## plot the EBT against the Gene density

mapsDF<-data.frame(GeneDensity=c(5:20),EBT=NA)
for (i in 1:16){
  mapsDF$EBT[i]<-maps[[i]]$EBT[nrow(maps[[i]])]
}
#mapsDF<-rbind(mapsDF,c(0,map0$EBT[nrow(map0)]))
bmp("geneDensityEBT.bmp",res=300,height=13,width=17,units='cm')
plot(mapsDF$GeneDensity,mapsDF$EBT,xlab="Gene density (genes per Mb)",ylab="Effective bottleneck time (generations)")
lmodel<-lm(mapsDF$EBT~mapsDF$GeneDensity)
abline(lmodel)
legend("topright",c(paste0("R-squared = ",round(summary(lmodel)$r.squared,2)),paste0("P-value = ",round(summary(lmodel)$coefficients[2,4],6))))
dev.off()