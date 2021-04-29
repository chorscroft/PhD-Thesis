
tiff(file="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/180412 Version after Peer Review/Version 3/figure4_v3tiff.tiff",res=300,height=17,width=17,units='cm',compression="lzw")

require(pROC)
h12rocP <- roc(h12MaxMatrixP[,1],h12MaxMatrixP[,2])
nSLrocP <- roc(nslMaxMatrixPStand[,1],nslMaxMatrixPStand[,2])
zalpharocP <- roc(zalphaMaxMatrixP[,1],zalphaMaxMatrixP[,2])
tselrocPT <- roc(tSelMaxMatrixTmed[,1],tSelMaxMatrixTmed[,2])

plot(zalpharocP,main="ROC curves")
plot.roc(h12rocP, add=TRUE, col="blue")
plot(nSLrocP, add=TRUE, col="red")
plot(tselrocPT, add=TRUE, col="green")
zauc <- auc(zalpharocP)
zauc5 <- auc(zalpharocP,partial.auc=c(0.95,1))/0.05
nauc <- auc(nSLrocP)
nauc5 <- auc(nSLrocP,partial.auc=c(0.95,1))/0.05
hauc <- auc(h12rocP)
hauc5 <- auc(h12rocP,partial.auc=c(0.95,1))/0.05
tauc <- auc(tselrocPT)
tauc5 <- auc(tselrocPT,partial.auc=c(0.95,1))/0.05

#zleg<-paste0(": AUC=",round(zauc*100),"% pAUC=",round(zauc5*100),"%")
#nleg<-paste0(": AUC=",round(nauc*100),"% pAUC=",round(nauc5*100),"%")
#hleg<-paste0("H12: AUC=",round(hauc*100),"% pAUC=",round(hauc5*100),"%")
#tleg<-paste0("TSel: AUC=",round(tauc*100),"% pAUC=",round(tauc5*100),"%")
zleg<-paste0(": ",round(zauc*100),"%, ",round(zauc5*100),"%")
nleg<-paste0(": ",round(nauc*100),"%, ",round(nauc5*100),"%")
hleg<-paste0("H12: ",round(hauc*100),"%, ",round(hauc5*100),"%")
tleg<-paste0("TSel: ",round(tauc*100),"%, ",round(tauc5*100),"%")


legend("bottomright",legend=c(as.expression(substitute("Z"[alpha]*zzleg,list(zzleg=zleg))),hleg,as.expression(substitute("nS"[L]*nnleg,list(nnleg=nleg))),tleg),col=c("black","blue","red","green"),lty=1,title="Method: AUC, pAUC")
dev.off()

#legend("bottomright",legend=c(hleg,as.expression(substitute("nS"[L]*nnleg,list(nnleg=nleg))),as.expression(substitute("Z"[alpha]*zzleg,list(zzleg=zleg))),tleg),col=c("black","red","blue","green"),lty=1)


# plot(1:10, type="n", xlab="", ylab="", main = "plot math & numbers")
# theta <- 1.23 ; mtext(bquote(hat(theta) == .(theta)), line= .25)
# for(i in 2:9)
#   text(i, i+1, substitute(list(xi, eta) == group("(",list(x,y),")"),
#                           list(x = i, y = i+1)))
# 
# 
# 
# 
# 
# 
# plot(simDF$SNPpos,simDF$h12,main="Neutral",xlab = "BP",ylab="H12")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# plot(simDF$SNPpos,simDF$h12,main="Selected",xlab = "BP",ylab="H12")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# 
# plot(simDF$SNPpos,simDF$zalpha,main="Neutral",xlab = "BP",ylab="zalpha")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# plot(simDF$SNPpos,simDF$zalpha,main="Selected",xlab = "BP",ylab="zalpha")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# abline(v=250000,col="green")
# 
# plot(simDF$SNPpos,simDF$zalpha,main="Neutral",xlab = "BP",ylab="zalpha")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# plot(simDF$SNPpos,simDF$zalpha,main="Selected",xlab = "BP",ylab="zalpha")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# 
# plot(simDF$SNPpos,simDF$SL,main="Neutral",xlab = "BP",ylab="nSL")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# plot(simDF$SNPpos,simDF$SL,main="Selected",xlab = "BP",ylab="nSL")
# abline(v=150000,col="red")
# abline(v=350000,col="red")
# 
# ### average graph of zalpha
# zalphaSelectedFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaSelected/"
# midsection <- c(150001,350000) #For analysis
# #Selected  
# sim <- 1
# mergeDF <- read.table(paste0(zalphaSelectedFolder,"zalpha_",sim,".txt"),col.names = c("SNPpos","zalpha"))
# for (sim in 2:nSims){
#   simDF <- read.table(paste0(zalphaSelectedFolder,"zalpha_",sim,".txt"),col.names = c("SNPpos","zalpha"))
#   mergeDF <- merge(mergeDF,simDF,by=1,all=TRUE)
# }
# averageZalpha <- cbind(mergeDF$SNPpos,apply(mergeDF[,c(-1)],1,sum,na.rm=TRUE),apply(mergeDF[,c(-1)],1,count))
# averageZalpha <- cbind(sapply(averageZalpha[,1],round1000),averageZalpha)
# binZalpha <- seq(0,max(averageZalpha[,1]),1000)
# binZalpha <- cbind(binZalpha,rep(0,length(binZalpha)),rep(0,length(binZalpha)))
# for (i in 1:nrow(averageZalpha)){
#   binZalpha[averageZalpha[i,1]/1000,2] <- binZalpha[averageZalpha[i,1]/1000,2]+averageZalpha[i,3]
#   binZalpha[averageZalpha[i,1]/1000,3] <- binZalpha[averageZalpha[i,1]/1000,3]+averageZalpha[i,4]
# }
# binZalpha <- cbind(binZalpha,rep(NA,nrow(binZalpha)))
# for (i in 1:nrow(binZalpha)){
#   if (binZalpha[i,3]>0){
#     binZalpha[i,4] <- binZalpha[i,2]/binZalpha[i,3]
#   }
# }
# plot(binZalpha[,1],binZalpha[,4])
# 
# plot(averageZalpha[averageZalpha[,1]>=midsection[1] & averageZalpha[,1]<=midsection[2],])
# 
