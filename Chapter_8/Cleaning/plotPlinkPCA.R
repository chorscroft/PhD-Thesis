### Reads in the PCA Analysis from Plink ###
### And then performs k-means clustering ###
#require(patchwork)
#require(tidygraph)
#require(ggplot2)
#require(PCDimension)
require(RColorBrewer)
require(Rtsne)

### PCA ###

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Dogs/Related Dogs/PCA")
pca <-read.table("plink.eigenvec")
pca[,1:2] <- lapply(pca[,1:2], as.character)

# Look at pca graph in two PCs
plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2")

# get sum of diagonals for total variance
rel <-read.table("plink.rel.diag")
sum(rel)

pcaval<-read.table("plink.eigenval")
plot(pcaval$V1,type="o",xlab="PC",ylab="Variances",xlim=c(1,20))

#check brokenStick method
#bs<-brokenStick(1:3700,3700)
#lines(bs*sum(rel),col="blue")
#dims<-bsDimension(pcaval$V1)
#AuerG<-AuerGervini(pcaval,dd=c(3700,111314))
#agDims<-agDimension(AuerG)

### k-means clustering ###

noPC<-7

#variance
sum(pcaval$V1[1:noPC])/sum(rel)

##choose number of clusters, given # PC as above:
wss<-(nrow(pca[,3:(2+noPC)])-1)*sum(apply(pca[3:(2+noPC)],2,var))
for (i in 2:15) {    #15 just a random high number - as max clusters
  wss[i] <- sum(kmeans(pca[3:(2+noPC)],i,1000,25)$withinss)
}

plot(1:15,wss,type="o",xlab="Clusters",ylab="Within groups Sum of Squares")

bmp("variancesclusters.bmp",res=300,height=10,width=17,units='cm')#,compression="lzw")
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(cex=0.75)
plot(pcaval$V1,type="o",xlab="PC",ylab="Variances",xlim=c(1,20))
title("A",adj=0)
plot(1:15,wss,type="o",xlab="Clusters",ylab="Within groups Sum of Squares")
title("B",adj=0)
dev.off()
par(mfrow=c(1,1))




noClust<-9
## pcakmeans with # PCs and # clusters
pcakmeans<-kmeans(pca[,3:(2+noPC)],noClust,1000,25)
table(pcakmeans$cluster)

#sort clusters and reassign
clustrank<-noClust+1-rank(table(pcakmeans$cluster))
newCluster<-clustrank[pcakmeans$cluster]
n<-table(newCluster)

bmp("PCAclusters.bmp",res=300,height=16,width=16,units='cm')#,compression="lzw")
plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2",col=brewer.pal(noClust,"Set1")[newCluster])
legend("bottomright",paste(c(1:noClust),": ",n),pch=1,col=brewer.pal(noClust,"Set1"),title = "cluster: n")
dev.off()
bmp("PCAclustersAll.bmp",res=300,height=16,width=16,units='cm')
pairs(pca[,3:(noPC+2)],col=brewer.pal(noClust,"Set1")[newCluster],labels=paste0("PC",c(1:noClust)),lower.panel = NULL)
legend("bottomleft",paste(c(1:noClust),": ",n),pch=1,col=brewer.pal(noClust,"Set1"),title = "cluster: n",xpd=TRUE)
dev.off()

indivToKeep<- pca$V1[newCluster==1]
write.table(data.frame(indivToKeep,indivToKeep),"indivToKeep.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

# get all of the individuals in each cluster
for (clust in 1:9){
  write.table(data.frame(pca$V1[newCluster==clust],pca$V1[newCluster==clust]),paste0("cluster_",clust,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE)
}

## Look at the three datasets
origDataset<-rep(1,3701)
origDataset[substr(pca$V1,1,9)=="breed_dog"]<-2
origDataset[substr(pca$V1,1,7)=="at_risk"]<-3
table(origDataset)

plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2",col=brewer.pal(noClust,"Set1")[newCluster],pch=origDataset)
plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2",col=origDataset)
legend("bottomright",c("discovery","breed_dog","at_risk"),col=c(1:3),pch=1)

bmp("PCAoriginaldata.bmp",res=300,height=16,width=16,units='cm')#,compression="lzw")
barplot(table(origDataset,newCluster),col=c("black","red","blue"),xlab="Cluster",ylab="Number of Dogs")
legend("topright",c("discovery","breed_dog","at_risk"),fill=c("black","red","blue"))
dev.off()
