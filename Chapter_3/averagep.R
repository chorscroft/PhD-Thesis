getMidP <- function(simFolder,prefix,nSLFlag){
  ### average graph of stat
  
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  nSims <- length(list.files(simFolder))
  midpoints<-vector("numeric",nSims)
  
  for (sim in 1:nSims){
    simDF <- read.table(paste0(simFolder,prefix,"_",sim,".txt"),header=isTRUE(nSLFlag>0))
    if (nSLFlag ==1){
      simDF <- simDF[simDF[,2]>0,] #pos
    }
    if (nSLFlag ==2){
      simDF <- simDF[simDF[,2]<0,] #neg
    }
    if (nSLFlag ==3){
      simDF[simDF[,2]<0,2] <- -simDF[simDF[,2]<0,2] #abs
    }  
    colnames(simDF) = c("SNPpos",prefix)
    #find SNP closest to midpoint 250000
    mid <- 0
    distance <- 999999
    print.default(nrow(simDF))
    for (i in 1:nrow(simDF)){
      if (abs(simDF[i,1]-250000)<distance && !is.na(simDF[i,2])){
        mid <- i
        distance <- abs(simDF[i,1]-250000)
      }
    }
    midpoints[sim]<- simDF[mid,2]
  }
  return(midpoints)
}
zalphaNeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaNeutral/"
zalphaSelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/zalpha/zalphaSelectedP/"
zalphaNeutralMid <- getMidP(zalphaNeutralFolder,"zalpha",0)
zalphaSelectMid <-getMidP(zalphaSelectedPFolder,"zalpha",0)
zalphattest <- t.test(zalphaNeutralMid,zalphaSelectMid)
zalphattest$p.value
zalphamannw <- wilcox.test(zalphaNeutralMid,zalphaSelectMid)
zalphamannw$p.value

h12NeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12Neutral/"
h12SelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/h12/h12SelectedP/"
h12NeutralMid <- getMidP(h12NeutralFolder,"h12",0)
h12SelectMid <-getMidP(h12SelectedPFolder,"h12",0)
h12ttest <- t.test(h12NeutralMid,h12SelectMid)
h12ttest$p.value
h12mannw<- wilcox.test(h12NeutralMid,h12SelectMid)
h12mannw$p.value

nslNeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslNeutralStand/"
nslSelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/nSL/nslSelectedPStand/"
nslNeutralMid <- getMidP(nslNeutralFolder,"nslNeutral",3)
nslSelectMid <-getMidP(nslSelectedPFolder,"nslSelect",3)
nslttest <- t.test(nslNeutralMid,nslSelectMid)
nslttest$p.value
nslmannw <-wilcox.test(nslNeutralMid,nslSelectMid)
nslmannw$p.value

tSelNeutralFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutputT/"
tSelSelectedPFolder<- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputPT/"
tSelNeutralMid <- getMidP(tSelNeutralFolder,"tSelOutput",0)
tSelSelectMid <-getMidP(tSelSelectedPFolder,"tSelOutput",0)
tSelttest <- t.test(tSelNeutralMid,tSelSelectMid)
tSelttest$p.value
tSelmannw <-wilcox.test(tSelNeutralMid,tSelSelectMid)
tSelmannw$p.value

hist(zalphaNeutralMid)
hist(zalphaSelectMid)
hist(h12NeutralMid)
hist(h12SelectMid)
hist(nslNeutralMid)
hist(nslSelectMid)
hist(tSelNeutralMid)
hist(tSelSelectMid)
