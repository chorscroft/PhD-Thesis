## Reads in the output from MSMS Tree run
## Input should be 
##   the simulation file name
##   the location to save the output


splitTreeOutputTrees <- function(filename,outFolder){
  txtFile <- file(filename, "r")
  blankflag <- FALSE  ## Looks for two empty rows
  dataFlag <- FALSE ##is true when the lines contain data
  countNoSims <- 0
  if (substr(outFolder,nchar(outFolder)-1,nchar(outFolder))!="/"){
    outFolder <- paste0(outFolder,"/")
  }
  
  while ( TRUE ) {
    line <- readLines(txtFile, n = 1)
    if ( length(line) == 0 || line=="" ) {
      if(blankflag == TRUE){
        break
      }
      blankflag <- TRUE
    } else {
      blankflag <- FALSE
      if (substr(line,1,2)=="//"){
        countNoSims <- countNoSims + 1
        matrixrow <- 0
        tempMatrix <- NA
        dataFlag <- TRUE
      } else if (dataFlag == TRUE){
        if (substr(line,1,1) != "["){
          write.table(tempMatrix,paste0(outFolder,"simTree_",countNoSims,".txt"),row.names = FALSE, col.names = FALSE)      
          dataFlag <- FALSE
        } else{
          if (matrixrow==0){
            matrixrow <- matrixrow + 1
            tempMatrix<- line
          } else {
            matrixrow <- matrixrow + 1
            tempMatrix <- c(tempMatrix,line)
          }
        }
      }
    }
  }
  close(txtFile)
  return("Complete")
}


msmsNeutralFile <- "/home/ch19g17/msms/msmsNeutralTreeT.txt"
msmsSelectFile<-  "/home/ch19g17/msms/msmsSelectTreePT.txt"
simNeutralFolder<-  "/home/ch19g17/msms/NeutralTreeT"
simSelectFolder<-  "/home/ch19g17/msms/SelectionTreePT"

## Split out the simulation files into folders
splitTreeOutputTrees(msmsNeutralFile,simNeutralFolder)
splitTreeOutputTrees(msmsSelectFile,simSelectFolder)
