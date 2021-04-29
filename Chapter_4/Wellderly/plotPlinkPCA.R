### Reads in the PCA Analysis from Plink ###
### And then performs k-means clustering ###

### PCA ###

setwd("//filestore.soton.ac.uk/users/ch19g17/mydocuments/Wellderly Filter")
pca <-read.table("plink.eigenvec")
pca[,1:2] <- lapply(pca[,1:2], as.character)

plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2")

pcaval<-read.table("plink.eigenval")
bmp("variances.bmp",res=300,height=14,width=14,units='cm')#,compression="lzw")
plot(pcaval$V1,type="o",xlab="PC",ylab="Variances")
dev.off()

### k-means clustering ###

##choose number of clusters, given 4 PC as above:
wss<-(nrow(pca[,3:6])-1)*sum(apply(pca[3:6],2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(pca[3:6],i,1000,25)$withinss)
}
bmp("calcclusters.bmp",res=300,height=14,width=14,units='cm')#,compression="lzw")
plot(1:15,wss,type="o",xlab="Clusters",ylab="Within groups Sum of Squares")
dev.off()

## pcakmeans with 4 PCs and 5 clusters
pcakmeans<-kmeans(pca[,3:6],5,1000,25)
table(pcakmeans$cluster)

#sort clusters and reassign
clustrank<-6-rank(table(pcakmeans$cluster))
newCluster<-clustrank[pcakmeans$cluster]
n<-table(newCluster)

bmp("PCAclusters.bmp",res=300,height=16,width=16,units='cm')#,compression="lzw")
plot(pca[,3],pca[,4],xlab="PC1",ylab="PC2",col=newCluster)
legend("bottomright",paste(c(1:5),": ",n),pch=1,col=c(1:5),title = "cluster: n")
dev.off()
indivToKeep<- pca$V1[newCluster==1]
write(indivToKeep,"indivToKeep.txt")
