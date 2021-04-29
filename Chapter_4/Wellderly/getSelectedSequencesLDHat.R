getRandomlySelectedSequences<-function(fileName,outputFileName){

  con <- file(fileName, "r")
  sequences<-NULL
  while (TRUE){
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    } else if (substr(line,1,1)==">") {
      sequences<-c(sequences,substr(line,2,nchar(line)))
    }
  }
  close(con)
  write.table(sequences,outputFileName,quote=FALSE,row.names = FALSE,col.names = FALSE)
  
}

setwd("/temp/hgig/EXOME_DATA/Clare/WellderlyChr22/LDhat")
getRandomlySelectedSequences("sites.txt","selectedSequences.txt")
