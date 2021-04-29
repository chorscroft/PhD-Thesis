
getSNPCountByChromosome<-function(tped,output){
  myTped<-read.table(tped,header=FALSE)
  myTped<-myTped[,1]
  df<-as.data.frame(table(myTped))
  write.table(df,output,quote=FALSE,row.names=FALSE)
}


tped<-"/temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/FinalCleanDataset/cleanDogDataset.tped"
output<-"chromosomeSummary.txt"

getSNPCountByChromosome(tped,output)