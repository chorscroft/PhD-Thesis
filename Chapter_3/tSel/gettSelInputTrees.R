## Gets the tree files and the bp files
## and combines the to make the input files
## for the tSel software
## Inputs: Location of bp and tree files, and output folder, and total bps

gettSelInputTrees <- function(simFolder, treeFolder, outputFolder, nbps){
 
  nSims <- length(list.files(treeFolder))
  
  ##add "/" to end of folder name if not already
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  if (substr(treeFolder,nchar(treeFolder),nchar(treeFolder))!="/"){
    treeFolder <- paste0(treeFolder,"/")
  }  
  if (substr(outputFolder,nchar(outputFolder),nchar(outputFolder))!="/"){
    outputFolder <- paste0(outputFolder,"/")
  } 
  
  for (sim in 1:nSims){
    
    simBPs <- unlist(read.table(paste0(simFolder,"simBPs_",sim,".txt")),use.names = FALSE)  
    simTree <-read.table(paste0(treeFolder,"simTree_",sim,".txt"), colClasses = "character")
    colnames(simTree)<-"Tree"
    for (i in 1:nrow(simTree)){
      b <- 1
      while(substr(simTree$Tree[i],b,b)!="]"){
        b<-b+1
      }
      simTree$snps[i] <- as.integer(substr(simTree$Tree[i],2,b-1))
      simTree$Tree[i] <- substr(simTree$Tree[i],b+1,nchar(simTree$Tree[i]))
    }
    nSNPs <- length(simBPs)
    tselDF <- data.frame(chr=integer(),position=integer(),averagePairs=double(),maxPairs=double(),medPairs=double(),varPairs=double(),fracMax=double(),quar1Pairs=double(),quar3Pairs=double())
    
    treeStep <- 0
    treeSnp <- 0
    
    for (snp in seq(1,nSNPs,20)){
          
      tselDF[snp,1] <- 1
      tselDF[snp,2] <- simBPs[snp]

      if (simBPs[snp] > treeSnp){
        while(simBPs[snp] > treeSnp){
          treeStep <- treeStep+ 1
          treeSnp <- treeSnp + simTree[treeStep,2]
        }
        #get tree features
        treeFeat <- getTreeFeatures(simTree[treeStep,1])
      }
      tselDF[snp,3] <- treeFeat$average
      tselDF[snp,4] <- treeFeat$max
      tselDF[snp,5] <- treeFeat$med
      tselDF[snp,6] <- treeFeat$var
      tselDF[snp,7] <- treeFeat$fracMax
      tselDF[snp,8] <- treeFeat$quar1
      tselDF[snp,9] <- treeFeat$quar3
    }
    write.table(tselDF,file=paste0(outputFolder,"tselInput_",sim,".txt"),row.names = FALSE)
  }
}
require(ape)
source("getTreeFeatures.R")

simNeutralFolder<-  "/home/ch19g17/msms/NeutralSimulations"
#simSelectFolder<-  "/home/ch19g17/msms/SelectionSimulationsP"
treeNeutralFolder <- "/home/ch19g17/msms/NeutralTreeT/"
#treeSelectFolder <- "/home/ch19g17/msms/SelectionTreePT/"
outputNeutralFolder <- "/home/ch19g17/tSel/tselNeutralInputT/"
#outputSelectFolder <- "/home/ch19g17/tSel/tselSelectInputPT/"
nbps <- 500000

gettSelInputTrees(simNeutralFolder,treeNeutralFolder,outputNeutralFolder,nbps)
#gettSelInputTrees(simSelectFolder,treeSelectFolder,outputSelectFolder,nbps)
