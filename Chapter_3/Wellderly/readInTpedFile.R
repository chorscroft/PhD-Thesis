#Reads tped file

gettPed <- function(tPed){
  #Do some checks
  myTpedFile<-read.table(tPed)
  colnames(myTpedFile) <- c("Chromosome","Variant","Position","bp",paste0("Indiv",rep(1:((ncol(myTpedFile)-4)/2),each=2),"_",1:2))
  return(myTpedFile)
}

