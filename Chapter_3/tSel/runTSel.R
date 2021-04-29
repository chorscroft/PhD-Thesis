## read in tSel inputs and run tsel

require(tsel)
runTSel <- function(inputFolder, outputFolder){
  ## add "/" to end of folder name if not already
  if (substr(inputFolder,nchar(inputFolder),nchar(inputFolder))!="/"){
    inputFolder <- paste0(inputFolder,"/")
  }
  ## add "/" to end of folder name if not already
  if (substr(outputFolder,nchar(outputFolder),nchar(outputFolder))!="/"){
    outputFolder <- paste0(outputFolder,"/")
  }
  
  ##count files in folder
  nSims <- length(list.files(inputFolder))
  for (sim in 1:nSims){
    tselInput<-getTselData(paste0(inputFolder,"tselInput_",sim,".txt"))
    #clean up
    tselInput$pos <- tselInput$pos[!is.na(tselInput$pos[,1]),]
    tselInput$features <- tselInput$features[!is.na(tselInput$features[,1]),]
    tseloutput<-getTselScore(tselInput,quantile = 100)
    write.table(tseloutput$score[,c(2,3)],file=paste0(outputFolder,"tselOutput_",sim,".txt"),col.names = FALSE, row.names = FALSE)
  }
}


#inputNFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralInput/"
#outputNFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutput/"
#inputSFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedInput/"
#outputSFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutput/"
#inputSFolder2 <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedInput2/"
#outputSFolder2 <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutput2/"
#inputSFolderP3 <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedInputP3/"
#outputSFolderP3 <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputP3/"
#runTSel(inputNFolder,outputNFolder)
#runTSel(inputSFolder,outputSFolder)
#runTSel(inputSFolder2,outputSFolder2)

inputNFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralInputT/"
outputNFolder <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/NeutralOutputT/"
inputSFolderPT <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedInputPT/"
outputSFolderPT <- "//filestore.soton.ac.uk/users/ch19g17/mydocuments/Review Paper/Simulations/tSel/SelectedOutputPT/"

runTSel(inputNFolder,outputNFolder)
runTSel(inputSFolderPT,outputSFolderPT)
