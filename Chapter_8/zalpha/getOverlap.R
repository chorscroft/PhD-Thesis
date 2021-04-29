require(readxl)

SNPs<-read_excel("H:/Dogs/Candidate Regions/Candidate Regions.xlsx","Overlap","A1:C231",col_types = c("numeric","numeric","numeric"))
PrevPub<-read_excel("H:/Dogs/Candidate Regions/Candidate Regions.xlsx","PrevPublished","A1:E535",col_types = c("numeric","numeric","numeric","text","text"))

addOverlap<-function(author,curr){
  if (curr==""){
    return(author)
  } else if (grepl(author,curr,fixed=TRUE)){
    return(curr)
  } else {
    return(paste0(curr,", ",author))
  }
}
overlap<-rep("",nrow(SNPs))
for (i in 1:nrow(SNPs)){
  for (j in 1:nrow(PrevPub)){
    if (PrevPub$Build[j]=="CanFam2"){
      if (is.na(PrevPub$to[j])){
        if (PrevPub$Chr[j]==SNPs$CHR[i] & PrevPub$Region[j]==SNPs$Canfam2[i]){
          overlap[i]<-addOverlap(PrevPub$Study[j],overlap[i])
        }
      } else {
        if (PrevPub$Chr[j]==SNPs$CHR[i] & PrevPub$Region[j]<=SNPs$Canfam2[i] & PrevPub$to[j]>=SNPs$Canfam2[i]){
          overlap[i]<-addOverlap(PrevPub$Study[j],overlap[i])
          if (i == 200){
            print.default(paste("j is",j))
          }
        }        
      }
    } else {
      if (is.na(PrevPub$to[j])){
        if (PrevPub$Chr[j]==SNPs$CHR[i] & PrevPub$Region[j]==SNPs$BP[i]){
          overlap[i]<-addOverlap(PrevPub$Study[j],overlap[i])
        }
      } else {
        if (PrevPub$Chr[j]==SNPs$CHR[i] & PrevPub$Region[j]<=SNPs$BP[i] & PrevPub$to[j]>=SNPs$BP[i]){
          overlap[i]<-addOverlap(PrevPub$Study[j],overlap[i])
        }        
      }      
    }
  }  
  
}
write.table(overlap,"H:/Dogs/Candidate Regions/overlap.txt",col.names = FALSE,quote=FALSE,sep = "\t")
