rates<-read.table("~/VariableRecombinationRates/recombinationRates.txt")
rates$V1<-rates$V1/1000


## multiply rates$v2 by 10000000 just for the y axis!
rates$V2<-rates$V2*10000000

png("~/VariableRecombinationRates/recombinationRatePlot.png",res=300,height=12,width=20,units='cm')
par(mar=c(5, 5, 4, 2) + 0.1)
plot(rates$V1,rates$V2,type="n",xlab="Kb",ylab=expression(paste("Recombination rate per bp per generation x10"^"-7")))

lines(x=c(0,rates$V1[1]),y=c(rates$V2[1],rates$V2[1]))
for (i in 1:(nrow(rates)-1)){
  lines(x=c(rates$V1[i],rates$V1[i]),y=c(rates$V2[i],rates$V2[i+1]))
  lines(x=c(rates$V1[i],rates$V1[i+1]),y=c(rates$V2[i+1],rates$V2[i+1]))
}
dev.off()
