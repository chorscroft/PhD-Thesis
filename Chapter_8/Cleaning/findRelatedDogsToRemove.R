## This piece of code finds related dogs in the plink.genome output
## and decides which ones to remove
## First, remove dogs with an identical pair
##   Find the dog with the most missing
##   Of these pick the top one
##   Remove dog
##   Repeat until all identicals are gone
## Next, remove the rest of the pairings
##   Find the dog with the most pairs
##   Of these find the dog with the most missing
##   Of these, pick the top one
##   Remove dog
##   Repeat until all related pairs are gone



## Read in Genome file
genome <- read.table("plink.genome",header=TRUE,stringsAsFactors=FALSE)
# Only keep required columns, and where PI_HAT > 0.2
genome<-genome[genome$PI_HAT>0.2,c(2,4,10)]

## Read in missingness file
missingness<- read.table("plink.imiss",header=TRUE,stringsAsFactors=FALSE)
# Only keep required columns
missingness<-missingness[,c(2,4)]

## Vector of dogs to remove
dogsToRemove<-NULL

##remove identical dogs first
while(sum(genome$PI_HAT==1)>0){
  ## Find unique dogs with a PI_HAT = 1 for at least one pairing
  identical<-data.frame(IID=unique(c(genome$IID1[genome$PI_HAT==1],genome$IID2[genome$PI_HAT==1])),stringsAsFactors=FALSE)
  ##merge on missing values
  identical<-merge(identical,missingness,by="IID",)
  ##merge on a count of pairs
  identical$pairs<-sapply(identical$IID,function(x){sum(genome$IID1==x)+sum(genome$IID2==x)})
  ##sort by missing
  identical<-identical[order(-identical$N_MISS),]
  ##write the top dog to the remove list
  removeDog<-identical$IID[1]
  dogsToRemove<-c(dogsToRemove,removeDog)
  ##remove the top dog's pairings from the genome file
  genome<-genome[genome$IID1 != removeDog & genome$IID2 != removeDog,]
}

##remove the rest
while(nrow(genome)>0){
  ## Find unique dogs
  uniqueDogs<-data.frame(IID=unique(c(genome$IID1,genome$IID2)),stringsAsFactors=FALSE)
  ##merge on missing values
  uniqueDogs<-merge(uniqueDogs,missingness,by="IID",)
  ##merge on a count of pairs
  uniqueDogs$pairs<-sapply(uniqueDogs$IID,function(x){sum(genome$IID1==x)+sum(genome$IID2==x)})
  ##sort by pairs count then by missing
  uniqueDogs<-uniqueDogs[order(-uniqueDogs$pairs,-uniqueDogs$N_MISS),]
  ##write the top dog to the remove list
  removeDog<-uniqueDogs$IID[1]
  dogsToRemove<-c(dogsToRemove,removeDog)
  ##remove the top dog's pairings from the genome file
  genome<-genome[genome$IID1 != removeDog & genome$IID2 != removeDog,]
}

## write dogs to remove out to a .txt file
write.table(cbind(dogsToRemove,dogsToRemove),"dogstoremove.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
