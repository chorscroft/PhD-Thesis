## Reads in the output from MSMS 
## Stores in a folder with each simulation over two files
## The first stores the bp positions
## The second stores the simulations as matrices:
##     - Row = each chromosome
##     - Col = each SNP (0 if ancestral, 1 if derived)s
## Input should be 
##   the simulation file name
##   the location to save the output


splitSimulationOutput <- function(filename,outFolder){
  txtFile <- file(filename, "r")
  flag <- FALSE  ## Looks for two empty rows
  countNoSims <- 0
  if (substr(outFolder,nchar(outFolder)-1,nchar(outFolder))!="/"){
    outFolder <- paste0(outFolder,"/")
  }
  
  while ( TRUE ) {
    line <- readLines(txtFile, n = 1)
    if ( length(line) == 0 || line=="" ) {
      if(flag == TRUE){
        break
      }
      flag <- TRUE
      #Insert matrix into list
      if (countNoSims>=1){
        write.table(tempMatrix,paste0(outFolder,"simMatrix_",countNoSims,".txt"),row.names = FALSE, col.names = FALSE)
      }
    } else {
      flag <- FALSE
      if (substr(line,1,2)=="ms"){
        bp <- as.numeric(strsplit(line," ")[[1]][8])
      } else if (line=="//"){
        countNoSims <- countNoSims + 1
        matrixrow <- 0
        matrixWidth <- 0
        tempMatrix <- NA
      } else if (substr(line,1,8) == "segsites"){
        matrixWidth <- as.numeric(substr(line,11,nchar(line)))
      } else if (substr(line,1,9) == "positions"){
        posVec <- round(as.numeric(strsplit(line," ")[[1]][c(-1)])*bp)
        write.table(posVec,paste0(outFolder,"simBPs_",countNoSims,".txt"),row.names = FALSE, col.names = FALSE)
      } else if (countNoSims >= 1){
        if (matrixrow==0){
          matrixrow <- matrixrow + 1
          tempMatrix<- as.numeric(strsplit(line,"")[[1]])
        } else {
          matrixrow <- matrixrow + 1
          tempMatrix <- rbind(tempMatrix,as.numeric(strsplit(line,"")[[1]]))
        }
      }
    }
  }
  close(txtFile)
  return("Complete")
}


msmsSelectFile<-  "/home/ch19g17/msms/msmsSelectP.txt"
simSelectFolder<-  "/home/ch19g17/msms/SelectionSimulationsP"

## Split out the simulation files into folders
splitSimulationOutput(msmsSelectFile,simSelectFolder)
