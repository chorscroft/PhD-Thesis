
## Graphs the Rho values for the four datasets across chromosome 22
analyseLDhat <- function(res,startPos,endPos,graphTitle,last){
  overallrho <- res[1,]
  res <- res[-1,]
  
  if (last==TRUE){
    plot(res$Loci,res$Mean_rho,type="n",xlim=c(16000.000,52000.000),ylim=c(0,125),las=1,ylab=graphTitle)
    title(xlab="chromosome 22 Kb",outer=TRUE)
  } else { 
    plot(res$Loci,res$Mean_rho,type="n",xlim=c(16000.000,52000.000),ylim=c(0,125),las=1,ylab=graphTitle,xaxt='n',xlab="")
  }
  for(i in 1:(nrow(res)-1)){
    #horizontal
    segments(x0=res$Loci[i]+startPos-0.001,y0=res$Mean_rho[i],x1=res$Loci[i+1]+startPos-0.001)
    #vertical
    segments(x0=res$Loci[i+1]+startPos-0.001,y0=res$Mean_rho[i],y1=res$Mean_rho[i+1])
  }
  #last horizontal
  segments(x0=res$Loci[nrow(res)]+startPos-0.001,y0=res$Mean_rho[nrow(res)],x1=endPos)
  
}

#tiff("testgraphs.tif",res = 300,compression = "lzw",height=14,width=12,units="cm")
par(mfrow=c(4,1),mar=c(0,5,0,0),oma=c(5,1,1,1)) 
### CEU ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22")
resC <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosC <- 20000.428
endPosC <- 51219.641
analyseLDhat(resC,startPosC,endPosC,"CEU",FALSE)

### Wellderly ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Wellderly")
resW <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosW <- 16051.295
endPosW <- 51218.377       ##16051.295+35167.083-0.001
analyseLDhat(resW,startPosW,endPosW,"Wellderly",FALSE)

### Baganda ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Baganda")
resB <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosB <- 16051.347
endPosB <- 51238.130       ##16051.347+35186.784-0.001
analyseLDhat(resB,startPosB,endPosB,"Baganda",FALSE)

### Zulu ###
setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/Zulu")
resZ <- read.table("res.txt",header=TRUE,colClasses = "numeric")
startPosZ <- 16051.347
endPosZ <- 51238.130       ##16051.347+35186.784-0.001
analyseLDhat(resZ,startPosZ,endPosZ,"Zulu",TRUE)
#dev.off()

max(resC$Mean_rho[-1])
max(resW$Mean_rho[-1])
max(resB$Mean_rho[-1])
max(resZ$Mean_rho[-1])

## Graphs the correlations between the four datasets (SNPs in common)

# getCorrelations<-function(common,col1,col2){
#   returnData<-list()
#   returnData$cor<-cor(common[,col1],common[,col2])
#   returnData$rsquared<-cor(common[,col1],common[,col2])^2 
#   returnData$pvalue<-cor.test(common[,col1],common[,col2])$p.value
#   return(returnData)
# }
res1<-resC[-1,]
res2<-resW[-1,]
res3<-resB[-1,]
res4<-resZ[-1,]
res1$Locus<-res1$Loci+startPosC-0.001
res2$Locus<-res2$Loci+startPosW-0.001
res3$Locus<-res3$Loci+startPosB-0.001
res4$Locus<-res4$Loci+startPosZ-0.001

common<-merge(res1[,c(2,6)],res2[,c(2,6)],by="Locus",all=FALSE)
common<-merge(common,res3[,c(2,6)],by="Locus",all=FALSE)
common<-merge(common,res4[,c(2,6)],by="Locus",all=FALSE)
colnames(common)<-c("Locus","C","W","B","Z")
# 
# corCW<-getCorrelations(common,2,3)
# corCB<-getCorrelations(common,2,4)
# corWB<-getCorrelations(common,3,4)
# corCZ<-getCorrelations(common,2,5)
# corWZ<-getCorrelations(common,3,5)
# corBZ<-getCorrelations(common,4,5)
# 
# ##plot correlations
# require(corrplot)
# Mr<-matrix(0,nrow=4,ncol=4,dimnames=list(c("CEU","Wellderly","Baganda","Zulu"),c("CEU","Wellderly","Baganda","Zulu")))
# Mr[1,2]<-corCW$rsquared
# Mr[1,3]<-corCB$rsquared
# Mr[1,4]<-corCZ$rsquared
# Mr[2,3]<-corWB$rsquared
# Mr[2,4]<-corWZ$rsquared
# Mr[3,4]<-corBZ$rsquared
# Mr[2,1]<-corCW$rsquared
# Mr[3,1]<-corCB$rsquared
# Mr[4,1]<-corCZ$rsquared
# Mr[3,2]<-corWB$rsquared
# Mr[4,2]<-corWZ$rsquared
# Mr[4,3]<-corBZ$rsquared
# Mp<-matrix(0,nrow=4,ncol=4,dimnames=list(c("CEU","Wellderly","Baganda","Zulu"),c("CEU","Wellderly","Baganda","Zulu")))
# Mp[1,2]<-corCW$pvalue
# Mp[1,3]<-corCB$pvalue
# Mp[1,4]<-corCZ$pvalue
# Mp[2,3]<-corWB$pvalue
# Mp[2,4]<-corWZ$pvalue
# Mp[3,4]<-corBZ$pvalue
# Mp[2,1]<-corCW$pvalue
# Mp[3,1]<-corCB$pvalue
# Mp[4,1]<-corCZ$pvalue
# Mp[3,2]<-corWB$pvalue
# Mp[4,2]<-corWZ$pvalue
# Mp[4,3]<-corBZ$pvalue
# par(mfrow=c(1,1))
# corrplot(Mr,method="color",type="upper",addCoef.col = "black",diag=FALSE,tl.col="black",cl.lim=c(0,1))




getKendallCorrelations<-function(common,col1,col2){
  returnData<-list()
  temp<-cor.test(common[,col1],common[,col2],method="kendall")
  returnData$cor<-temp$est
  returnData$pvalue<-temp$p.value
  return(returnData)
}
KcorCW<-getKendallCorrelations(common,2,3)
KcorCB<-getKendallCorrelations(common,2,4)
KcorWB<-getKendallCorrelations(common,3,4)
KcorCZ<-getKendallCorrelations(common,2,5)
KcorWZ<-getKendallCorrelations(common,3,5)
KcorBZ<-getKendallCorrelations(common,4,5)

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

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/LDHat/CEU22/LDMAP comparison")
bmp("kendallCompare.bmp",res=300,height=14,width=15,units='cm')
corrplot(Mr,method="color",type="upper",addCoef.col = "black",diag=FALSE,tl.col="black",p.mat=Mp,sig.level=0.01,cl.lim=c(-1,1))
dev.off()