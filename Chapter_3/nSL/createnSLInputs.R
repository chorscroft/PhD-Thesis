## creates the nSL input files
## requires the msms output folder (after split)
## requires the output foler location

createnSLInputs <- function(simFolder, outputFolder){
  
  ## add "/" to end of folder name if not already
  if (substr(simFolder,nchar(simFolder),nchar(simFolder))!="/"){
    simFolder <- paste0(simFolder,"/")
  }
  ## add "/" to end of folder name if not already
  if (substr(outputFolder,nchar(outputFolder),nchar(outputFolder))!="/"){
    outputFolder <- paste0(outputFolder,"/")
  }
  
  ## count files in folder
  nSims <- length(list.files(simFolder))/2
  ## Loop over each simulation in turn
  for (sim in 1:nSims){
    ## Read in bp locations
    simBPs <- unlist(read.table(paste0(simFolder,"simBPs_",sim,".txt")),use.names = FALSE)  
    ## Read in matrix of snp status:
    ## rows are chromosomes
    ## columns are each snp
    ## 0 = ancenstral; 1 = derived
    simMatrix <-data.matrix(read.table(paste0(simFolder,"simMatrix_",sim,".txt")))
    colnames(simMatrix)<-NULL
    
    ## get the number of snps
    nSNPs <- length(simBPs) 
    nChrm <- nrow(simMatrix)
    
    ## Write samples file
    write(paste0("CH",formatC(1:nChrm,width=5,format="d",flag="0")), file=paste0(outputFolder,"samfile_",sim,".txt"))
    
    ## Write adfile
    write.table(cbind(paste0("snp",1:nSNPs),rep("A",nSNPs),rep("T",nSNPs)),file=paste0(outputFolder,"adfile_",sim,".txt"),sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
    
    ## Write Haplotype file
    textVector <- c("sim",sim,"snps:")
    snpsVector <- vector(mode="character",nSNPs)
    for (i in 1:nSNPs){
      snpsVector[i] <- paste0("  - snp",i,": ",simBPs[i])
    }
    hapsVector <- vector(mode="character",nChrm)
    for (i in 1:nChrm){
      tempHap <- simMatrix[i,]
      tempHap[tempHap==0] <- "A"  #use A for ancestral
      tempHap[tempHap==1] <- "T"  #use T for derived
      tempHap <- paste(tempHap,collapse = "")
      hapsVector[i] <- paste0("  - CH",formatC(i,width=5,format="d",flag="0"),": ",tempHap)
    }
    write(c(textVector,snpsVector,"phased_haplotypes:",hapsVector),file=paste0(outputFolder,"hapfile_",sim,".txt"))
  }
}
