merged<-read.table("H:/Dogs/Candidate Regions/merged.txt",header = TRUE)

## get recombination map
chrMap<-list()

for (i in 1:38){
  chrMap[[i]]<-read.table(paste0("\\\\filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/LDhat/Maps/cM_map_chr",i,".txt"),stringsAsFactors = FALSE,header=TRUE)
}

## get overlaps with previously published data
require(readxl)
PrevPub<-read_excel("H:/Dogs/PreviouslyPublished/PreviouslyPublished.xlsx","CanFam3.1","A1:E521",col_types = c("numeric","numeric","numeric","text","numeric"))
PrevPub$size[PrevPub$size==0]<-1

## get diversity stat
L_plus_R<-read.table("H:/Dogs/Candidate Regions/L_plus_R.txt",stringsAsFactors = FALSE,header=TRUE)


## Gvis
require(Gviz)
require(ensembldb)
require(AnnotationHub)
ah <- AnnotationHub()
## Query for all available EnsDb databases
#bob<-query(ah, "EnsDb")
#d<-display(ah)

dogdb <- ah["AH79907"]
## Create a EnsDb database file from this.
DbFile <- ensDbFromAH(dogdb)
## We can either generate a database package, or directly load the data
edb <- EnsDb(DbFile)
seqlevelsStyle(edb) <- "UCSC"




plotRegion<-function(chrom,bplow,bphigh,folder){
  #set colour of title panels
  scheme<-getScheme()
  scheme$GdObject$background.title="darkgrey"
  addScheme(scheme, "myScheme")
  options(Gviz.scheme="myScheme")
  getOption("Gviz.scheme")
  
  #get y limits
  Zroe_ylim<-c(min(merged$Zroe),max(merged$Zroe))
  Zb_ylim<-c(min(merged$Zb),max(merged$Zb))
  ZroeMZbroe_ylim<-c(min(merged$ZroeMZbroe,na.rm = TRUE),max(merged$ZroeMZbroe,na.rm = TRUE))
  ZbOZbb_ylim<-c(min(merged$ZbOZbb,na.rm = TRUE),max(merged$ZbOZbb,na.rm = TRUE))
  
  recomb_ylim<-c(0,max(sapply(1:38,function(x) max(chrMap[[x]]$cMMb))))
  div_ylim<-c(0,3.5)
  
  #get median
  Zroe_median<-median(merged$Zroe)
  Zb_median<-median(merged$Zb)
  ZroeMZbroe_median<-median(merged$ZroeMZbroe,na.rm=TRUE)
  ZbOZbb_median<-median(merged$ZbOZbb,na.rm=TRUE)
  
  #get outliers
  Zroe_top01<-quantile(merged$Zroe,0.999)
  Zb_top01<-quantile(merged$Zb,0.999)
  div_bot01<-quantile(log10(L_plus_R$LplusR+1),0.001)
  
  grtrack<-getGeneRegionTrackForGviz(edb,chromosome=chrom,start=bplow,end=bphigh,featureIs="tx_biotype")
  grtrack@elementMetadata@listData$symbol[is.na(grtrack@elementMetadata@listData$symbol)]<-grtrack@elementMetadata@listData$gene[is.na(grtrack@elementMetadata@listData$symbol)]
  genesInRegion<-unique(grtrack@elementMetadata@listData$symbol)
  grtrack<-GeneRegionTrack(grtrack, collapseTranscripts = "meta", shape = "box", fill="grey", 
                           transcriptAnnotation = "symbol",from=bplow,to=bphigh,name="Genes")
  
  ideoTrack <- IdeogramTrack(genome = "canFam3", chromosome = paste0("chr",chrom))
  gtrack <- GenomeAxisTrack()
  
  region<-merged$CHR==chrom & merged$BP>=bplow & merged$BP<=bphigh
  zroetrack<-DataTrack(start=merged$BP[region],end=merged$BP[region],data=merged$Zroe[region],name="Za_r/E[r]",genome="canFam3",chromosome = paste0("chr",chrom),ylim=Zroe_ylim)
  zroetoptrack<-DataTrack(start=c(bplow,bphigh),end=c(bplow,bphigh),data=c(Zroe_top01,Zroe_top01),name="Za_r/E[r]",genome="canFam3",chromosome = paste0("chr",chrom),ylim=Zroe_ylim,type="l",col="red",lty=2)
  zroecombtrack<-OverlayTrack(trackList=list(zroetrack, zroetoptrack))
  
  zbtrack<-DataTrack(start=merged$BP[region],end=merged$BP[region],data=merged$Zb[region],name="Za_BetaC",genome="canFam3",chromosome = paste0("chr",chrom),ylim=Zb_ylim,legend=TRUE)
  zbtoptrack<-DataTrack(start=c(bplow,bphigh),end=c(bplow,bphigh),data=c(Zb_top01,Zb_top01),name="Za_BetaC",genome="canFam3",chromosome = paste0("chr",chrom),ylim=Zb_ylim,type="l",col="red",lty=2,groups=c("Top 0.1%"),legend=TRUE)
  zbcombtrack<-OverlayTrack(trackList=list(zbtrack, zbtoptrack))
  
  ZroeMZbroetrack<-DataTrack(start=merged$BP[region],end=merged$BP[region],data=merged$ZroeMZbroe[region],name="Za_r/E[r]\n-Zb_r/E[r]",genome="canFam3",chromosome = paste0("chr",chrom),ylim=ZroeMZbroe_ylim,col="darkgreen")
  ZroeMZbroemedtrack<-DataTrack(start=c(bplow,bphigh),end=c(bplow,bphigh),data=c(ZroeMZbroe_median,ZroeMZbroe_median),name="Za_r/E[r]\n-Zb_r/E[r]",genome="canFam3",chromosome = paste0("chr",chrom),ylim=ZroeMZbroe_ylim,type="l",col="black",lty=2)
  ZroeMZbroecombtrack<-OverlayTrack(trackList=list(ZroeMZbroetrack, ZroeMZbroemedtrack))
  
  ZbOZbbtrack<-DataTrack(start=merged$BP[region],end=merged$BP[region],data=merged$ZbOZbb[region],name="Za_BetaC\n/Zb_BetaC",genome="canFam3",chromosome = paste0("chr",chrom),ylim=ZbOZbb_ylim,col="darkgreen")
  ZbOZbbmedtrack<-DataTrack(start=c(bplow,bphigh),end=c(bplow,bphigh),data=c(ZbOZbb_median,ZbOZbb_median),name="Za_BetaC\n/Zb_BetaC",genome="canFam3",chromosome = paste0("chr",chrom),ylim=ZbOZbb_ylim,type="l",col="black",lty=2)
  ZbOZbbcombtrack<-OverlayTrack(trackList=list(ZbOZbbtrack, ZbOZbbmedtrack))
  
  
  recombRegion<-chrMap[[chrom]]$loci>=bplow & chrMap[[chrom]]$loci<=bphigh
  recombtrack<-DataTrack(start=chrMap[[chrom]]$loci[recombRegion],end=chrMap[[chrom]]$loci[recombRegion],data=chrMap[[chrom]]$cMMb[recombRegion],name="cM/Mb",genome="canFam3",chromosome = paste0("chr",chrom),type="l",ylim=recomb_ylim)
  
  diverseRegion<-L_plus_R$CHR==chrom & L_plus_R$BP>=bplow & L_plus_R$BP<=bphigh
  divtrack<-DataTrack(start=L_plus_R$BP[diverseRegion],end=L_plus_R$BP[diverseRegion],data=log10(L_plus_R$LplusR[diverseRegion]+1),name="L_plus_R",genome="canFam3",chromosome = paste0("chr",chrom),ylim=div_ylim)
  divbottrack<-DataTrack(start=c(bplow,bphigh),end=c(bplow,bphigh),data=c(div_bot01,div_bot01),name="L_plus_R",genome="canFam3",chromosome = paste0("chr",chrom),ylim=div_ylim,type="l",col="red",lty=2,groups=c("Bottom 0.1%"),legend=TRUE)
  divcombtrack<-OverlayTrack(trackList=list(divtrack, divbottrack))
  
  
  tempprevpub<-PrevPub[PrevPub$Chr==chrom,]
  prevpubtrack<-AnnotationTrack(start=tempprevpub$Region,width=tempprevpub$size,chromosome = chrom,strand="*",id=tempprevpub$Study,genome="canFam3",name="Pub",fill="pink",featureAnnotation = "id",fontsize.feature=6,fontcolor.feature="black")
  
  png(paste0("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/",folder,"/chr",chrom,"_",bplow,"_",bphigh,".png"),res=300,height=20,width=14,units='cm')
  plotTracks(list(ideoTrack,gtrack,grtrack,zroecombtrack,zbcombtrack,ZroeMZbroecombtrack,ZbOZbbcombtrack,recombtrack,divcombtrack,prevpubtrack), from = bplow, to = bphigh)
  dev.off()
  

}


options(scipen=999)

## Set region to plot
buffer<-1000000

## 11 Regions validated both ways
plotRegion(2,61876498-buffer,61901702+buffer,"Graphs")
plotRegion(4,57366377-buffer,57366377+buffer,"Graphs")
plotRegion(5,4064061-buffer,4093514+buffer,"Graphs")
plotRegion(6,33510473-buffer,33510473+buffer,"Graphs")
plotRegion(7,24652821-buffer,24664438+buffer,"Graphs")
plotRegion(10,44372549-buffer,44388924+buffer,"Graphs")
plotRegion(11,54324689-buffer,54391443+buffer,"Graphs")
plotRegion(15,20317533-buffer,20317533+buffer,"Graphs")
plotRegion(16,7462818-buffer,7462818+buffer,"Graphs")
plotRegion(19,4813917-buffer,6590666+buffer,"Graphs")  
plotRegion(30,1558195-buffer,1732646+buffer,"Graphs")


## Regions unique to this study
plotRegion(1,43001368-buffer,43001368+buffer,"GraphsNoOverlap")
plotRegion(1,96115461-buffer,96115461+buffer,"GraphsNoOverlap")
plotRegion(2,71434345-buffer,71434345+buffer,"GraphsNoOverlap")
plotRegion(3,17490492-buffer,17516194+buffer,"GraphsNoOverlap")
plotRegion(3,72708942-buffer,72708942+buffer,"GraphsNoOverlap")
plotRegion(4,17518453-buffer,17518453+buffer,"GraphsNoOverlap")
plotRegion(4,57345395-buffer,57345395+buffer,"GraphsNoOverlap")
plotRegion(5,6838932-buffer,6859691+buffer,"GraphsNoOverlap")
plotRegion(5,40202215-buffer,40202215+buffer,"GraphsNoOverlap")
plotRegion(8,7735497-buffer,7735497+buffer,"GraphsNoOverlap")
plotRegion(9,29752455-buffer,29752455+buffer,"GraphsNoOverlap")
plotRegion(10,46053118-buffer,46053118+buffer,"GraphsNoOverlap")
plotRegion(12,26284264-buffer,26284264+buffer,"GraphsNoOverlap")
plotRegion(12,31691990-buffer,31835704+buffer,"GraphsNoOverlap")
plotRegion(14,8117811-buffer,8117811+buffer,"GraphsNoOverlap")
plotRegion(17,3753156-buffer,3753156+buffer,"GraphsNoOverlap")
plotRegion(18,29595073-buffer,29595073+buffer,"GraphsNoOverlap")
plotRegion(19,7095253-buffer,7122489+buffer,"GraphsNoOverlap")
plotRegion(20,13387022-buffer,13387022+buffer,"GraphsNoOverlap")
plotRegion(22,11073667-buffer,12039716+buffer,"GraphsNoOverlap")
plotRegion(22,18774821-buffer,19925395+buffer,"GraphsNoOverlap")
plotRegion(22,31194138-buffer,31347124+buffer,"GraphsNoOverlap")
plotRegion(26,22151015-buffer,22156289+buffer,"GraphsNoOverlap")
plotRegion(27,44328723-buffer,44328723+buffer,"GraphsNoOverlap")
plotRegion(30,4822803-buffer,4822803+buffer,"GraphsNoOverlap")
plotRegion(32,24657487-buffer,25070561+buffer,"GraphsNoOverlap")

### genes
#IGF1 15: 41,202,518-41,275,794 
plotRegion(15,41202518-buffer,41275794+buffer,"Graphs_test")
#MBP 1: 2,846,589-2,951,860 
plotRegion(1,2846589-buffer,2951860+buffer,"Graphs_test")
#SEMA3D 18: 24,262,125-24,394,155 
plotRegion(18,24262125-buffer,24394155+buffer,"Graphs_test")


### interesting regions
#10 2000000-9000000 overlapped a lot
plotRegion(10,2000000,9000000,"Graphs_test")
#18 0-2000000 overlapped a lot
plotRegion(18,0,2000000+buffer,"Graphs_test")
#31	3000000-4000000 overlapped a lot
plotRegion(31,3000000-buffer,4000000+buffer,"Graphs_test")


#### get genes

getGenesInRegion<-function(chrom,bplow,bphigh,folder){
  
  grtrack<-getGeneRegionTrackForGviz(edb,chromosome=chrom,start=bplow,end=bphigh,featureIs="tx_biotype")
  grtrack@elementMetadata@listData$symbol[is.na(grtrack@elementMetadata@listData$symbol)]<-grtrack@elementMetadata@listData$gene[is.na(grtrack@elementMetadata@listData$symbol)]
  genesInRegion<-unique(grtrack@elementMetadata@listData$symbol)
  geneEnsembl<-grtrack@elementMetadata@listData$gene
  
  write.table(genesInRegion,paste0("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/",folder,"/genes/chr",chrom,"_",bplow,"_",bphigh,".txt"),quote = FALSE,col.names = FALSE,row.names = FALSE)
  named<-genesInRegion[substr(genesInRegion,1,7)!="ENSCAFG"]
  named<-paste(named,collapse =", ")
  write.table(named,paste0("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/",folder,"/genes/chr",chrom,"_",bplow,"_",bphigh,".txt"),quote = FALSE,col.names = FALSE,row.names = FALSE,append = TRUE)

  write.table(geneEnsembl,paste0("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Candidate Regions/ensemblgenes.txt"),quote = FALSE,col.names = FALSE,row.names = FALSE,append = TRUE)
  
  
}


## Set region to get genes for
buffer<-250000

## 11 Regions validated both ways
getGenesInRegion(2,61876498-buffer,61901702+buffer,"Graphs")
getGenesInRegion(4,57366377-buffer,57366377+buffer,"Graphs")
getGenesInRegion(5,4064061-buffer,4093514+buffer,"Graphs")
getGenesInRegion(6,33510473-buffer,33510473+buffer,"Graphs")
getGenesInRegion(7,24652821-buffer,24664438+buffer,"Graphs")
getGenesInRegion(10,44372549-buffer,44388924+buffer,"Graphs")
getGenesInRegion(11,54324689-buffer,54391443+buffer,"Graphs")
getGenesInRegion(15,20317533-buffer,20317533+buffer,"Graphs")
getGenesInRegion(16,7462818-buffer,7462818+buffer,"Graphs")
getGenesInRegion(19,4813917-buffer,6590666+buffer,"Graphs")  
getGenesInRegion(30,1558195-buffer,1732646+buffer,"Graphs")


## Regions unique to this study
getGenesInRegion(1,43001368-buffer,43001368+buffer,"GraphsNoOverlap")
getGenesInRegion(1,96115461-buffer,96115461+buffer,"GraphsNoOverlap")
getGenesInRegion(2,71434345-buffer,71434345+buffer,"GraphsNoOverlap")
getGenesInRegion(3,17490492-buffer,17516194+buffer,"GraphsNoOverlap")
getGenesInRegion(3,72708942-buffer,72708942+buffer,"GraphsNoOverlap")
getGenesInRegion(4,17518453-buffer,17518453+buffer,"GraphsNoOverlap")
getGenesInRegion(4,57345395-buffer,57345395+buffer,"GraphsNoOverlap")
getGenesInRegion(5,6838932-buffer,6859691+buffer,"GraphsNoOverlap")
getGenesInRegion(5,40202215-buffer,40202215+buffer,"GraphsNoOverlap")
getGenesInRegion(8,7735497-buffer,7735497+buffer,"GraphsNoOverlap")
getGenesInRegion(9,29752455-buffer,29752455+buffer,"GraphsNoOverlap")
getGenesInRegion(10,46053118-buffer,46053118+buffer,"GraphsNoOverlap")
getGenesInRegion(12,26284264-buffer,26284264+buffer,"GraphsNoOverlap")
getGenesInRegion(12,31691990-buffer,31835704+buffer,"GraphsNoOverlap")
getGenesInRegion(14,8117811-buffer,8117811+buffer,"GraphsNoOverlap")
getGenesInRegion(17,3753156-buffer,3753156+buffer,"GraphsNoOverlap")
getGenesInRegion(18,29595073-buffer,29595073+buffer,"GraphsNoOverlap")
getGenesInRegion(19,7095253-buffer,7122489+buffer,"GraphsNoOverlap")
getGenesInRegion(20,13387022-buffer,13387022+buffer,"GraphsNoOverlap")
getGenesInRegion(22,11073667-buffer,12039716+buffer,"GraphsNoOverlap")
getGenesInRegion(22,18774821-buffer,19925395+buffer,"GraphsNoOverlap")
getGenesInRegion(22,31194138-buffer,31347124+buffer,"GraphsNoOverlap")
getGenesInRegion(26,22151015-buffer,22156289+buffer,"GraphsNoOverlap")
getGenesInRegion(27,44328723-buffer,44328723+buffer,"GraphsNoOverlap")
getGenesInRegion(30,4822803-buffer,4822803+buffer,"GraphsNoOverlap")
getGenesInRegion(32,24657487-buffer,25070561+buffer,"GraphsNoOverlap")


