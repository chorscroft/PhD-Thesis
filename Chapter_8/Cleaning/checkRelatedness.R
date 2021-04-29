
temp<-read.table("plink.genome",header=TRUE,stringsAsFactors=FALSE)

max(temp$PI_HAT)
sum(temp$PI_HAT>0.1875)
sum(temp$PI_HAT>0.2)
histcounts<-hist(temp$PI_HAT)
df<-data.frame(PI_HAT_interval=paste(histcounts$breaks[-21],"-",histcounts$breaks[-1]),count=histcounts$counts)

potentRel<-temp[temp$PI_HAT>0.2,]
identical<-temp[temp$PI_HAT==1,]
rm(temp)

relIDs<-c(potentRel$IID1,potentRel$IID2)
tableDogs<-table(relIDs)
row.names(tableDogs)

identIDs<-c(identical$IID1,identical$IID2)
tableIdent<-table(identIDs)


write.csv(potentRel,"potentiallyRelated.csv",row.names=FALSE,quote=FALSE)
