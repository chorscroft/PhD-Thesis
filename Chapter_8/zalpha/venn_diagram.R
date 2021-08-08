library(VennDiagram)

set1<-1:143
set2<-88:230
venn.diagram(
  x = list(set1, set2),
  filename="//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/venn.png",
  category.names = c(expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"))) , expression(paste("Z"[alpha]^BetaCDF))),
  col=c("red","blue"),
  fill = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
  cex=2,
  cat.cex=2,
  cat.pos=c(-90,90),
  cat.dist=c(0.10,0.12),
  margin=0.1
)

