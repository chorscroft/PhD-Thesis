#require(ggraph)
#require(tidygraph)
#require(patchwork)
require(qgraph)

setwd("\\\\filestore.soton.ac.uk\\users\\ch19g17\\mydocuments\\Dogs\\Related Dogs")
potRel<-read.csv("potentiallyRelated.csv",stringsAsFactors = FALSE)

# ## graph the data nicely
# ggraphData<-function(dataset,title){
#   gGraph<-as_tbl_graph(data.frame(from=dataset$IID1,to=dataset$IID2,Identical=ifelse(dataset$PI_HAT==1,TRUE,FALSE)))
#   ggraph(gGraph,layout="mds")+
#     geom_edge_link(aes(col=Identical),edge_width=1)+
#     geom_node_point(size=0.5)+
#     ggtitle(title)+
#     scale_edge_color_manual(values=c("TRUE"="blue","FALSE"="red"))
# }
# 
# 
# ## Just Identical
 identical<-potRel[potRel$PI_HAT==1,]
# identIDs<-table(c(identical$IID1,identical$IID2))
# sum(identIDs==2)
# 
# 
# identGraph<-ggraphData(identical,"Identical Dogs")
# 
# ## All potentially related
# 
# allPotRelGraph<-ggraphData(potRel,"All Potentially Related")
# 
# 
# ## Only non-identical
# nonIdent<-potRel[potRel$PI_HAT!=1,]
# nonIdentGraph<-ggraphData(nonIdent,"Only Non-Identical")
# 
# jpeg("AllPotentiallyRelated.jpg",res=300,height=14,width=16,units='cm',quality=100)
# allPotRelGraph + theme_graph(plot_margin = margin(0, 0, 0, 0))
# dev.off()
# 
# jpeg("IdenticalDogs.jpg",res=300,height=14,width=16,units='cm',quality=100)
# identGraph + theme_graph(plot_margin = margin(0, 0, 0, 0))
# dev.off()
# # #graph with scale for PI_HAT
# #   gGraph<-as_tbl_graph(data.frame(from=nonIdent$IID1,to=nonIdent$IID2,PI_HAT=nonIdent$PI_HAT))
# #   ggraph(gGraph,layout="mds")+
# #     geom_edge_link(aes(col=PI_HAT),edge_width=1)+
# #     geom_node_point()+
# #     ggtitle("Only Non-Identical")+
# #     scale_edge_color_gradient(low="red",high="blue")
# 
# 
# 
# 
# 
# ## Only within the at_risk dataset
 at_risk<-potRel[substr(potRel$IID1,1,7)=="at_risk" & substr(potRel$IID2,1,7)=="at_risk",]
# risk<-ggraphData(at_risk,"at_risk")
# 
# ## Only within the breed_dog dataset
 breed_dog<-potRel[substr(potRel$IID1,1,9)=="breed_dog" & substr(potRel$IID2,1,9)=="breed_dog",]
# breed<-ggraphData(breed_dog,"breed_dog")
# 
# ##Only within the discovery dataset
 discovery<-potRel[nchar(as.vector(potRel$IID1))==4 & nchar(as.vector(potRel$IID2))==4,]
# disc<-ggraphData(discovery,"discovery")
# 
# disc + risk + breed & theme(legend.position = "bottom")
# 
# ## between at_risk and breed_dog
 at_risk_breed_dog<-potRel[(substr(potRel$IID1,1,7)=="at_risk" & substr(potRel$IID2,1,9)=="breed_dog") | (substr(potRel$IID2,1,7)=="at_risk" & substr(potRel$IID1,1,9)=="breed_dog"),]
# riskbreed<-ggraphData(at_risk_breed_dog,"at_risk and breed_dog")
# 
# ## between at_risk and discovery
 at_risk_discovery<-potRel[(substr(potRel$IID1,1,7)=="at_risk" & nchar(as.vector(potRel$IID2))==4) | (substr(potRel$IID2,1,7)=="at_risk" & nchar(as.vector(potRel$IID1))==4),]
# riskdisc<-ggraphData(at_risk_discovery,"discovery and at_risk")
# 
# ## between breed_dog and discovery
 breed_dog_discovery<-potRel[(substr(potRel$IID1,1,9)=="breed_dog" & nchar(as.vector(potRel$IID2))==4) | (substr(potRel$IID2,1,9)=="breed_dog" & nchar(as.vector(potRel$IID1))==4),]
# breeddisc<-ggraphData(breed_dog_discovery,"discovery and breed_dog")
# 
# riskdisc + breeddisc + riskbreed & theme(legend.position = "bottom")
# 
# 
# disc<-disc+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0)) + theme(legend.position = "none")
# risk<-risk+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0)) + theme(legend.position = "none")
# breed<-breed+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0)) + theme(legend.position = "none")
# riskdisc<-riskdisc+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0))
# breeddisc<-breeddisc+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0))
# riskbreed<-riskbreed+theme_graph(title_size = 10,plot_margin = margin(0, 0, 0, 0))
# jpeg("relatedWithinAndBetween.jpg",res=300,height=14,width=16,units='cm',quality=100)
# (disc + risk + breed) / ((riskdisc +breeddisc + riskbreed) + plot_layout(guides = 'collect')  & theme(legend.position = 'bottom')) + plot_annotation(tag_levels = 'A')
# dev.off()
# 
# 
## find clusters

clusters<-list()
dogs<-data.frame(IID=unique(c(potRel$IID1,potRel$IID2)),cluster=0)
noClusters <- 0
for (i in 1:nrow(potRel)){
  dog1<-potRel$IID1[i]
  dog2<-potRel$IID2[i]
  if (dogs$cluster[dogs$IID==dog1] == 0 & dogs$cluster[dogs$IID==dog2] == 0){
    noClusters<-noClusters+1
    dogs$cluster[dogs$IID==dog1] <- noClusters
    dogs$cluster[dogs$IID==dog2] <- noClusters
    clusters[[noClusters]]<-c(dog1,dog2)
  } else if (dogs$cluster[dogs$IID==dog1] == 0){
    dogs$cluster[dogs$IID==dog1]<-dogs$cluster[dogs$IID==dog2]
    clusters[[dogs$cluster[dogs$IID==dog1]]]<-c(clusters[[dogs$cluster[dogs$IID==dog1]]],dog1)
  } else if (dogs$cluster[dogs$IID==dog2] == 0){
    dogs$cluster[dogs$IID==dog2]<-dogs$cluster[dogs$IID==dog1]
    clusters[[dogs$cluster[dogs$IID==dog2]]]<-c(clusters[[dogs$cluster[dogs$IID==dog2]]],dog2)
  } else if (dogs$cluster[dogs$IID==dog1] != 0 & dogs$cluster[dogs$IID==dog2] != 0 & dogs$cluster[dogs$IID==dog1] != dogs$cluster[dogs$IID==dog2]){
    # merge clusters
    if (dogs$cluster[dogs$IID==dog1] < dogs$cluster[dogs$IID==dog2]){
      clusters[[dogs$cluster[dogs$IID==dog1]]]<-c(clusters[[dogs$cluster[dogs$IID==dog1]]],clusters[[dogs$cluster[dogs$IID==dog2]]])
      clusters[[dogs$cluster[dogs$IID==dog2]]]<-NA
      dogs$cluster[dogs$cluster==dogs$cluster[dogs$IID==dog2]]<-dogs$cluster[dogs$IID==dog1]
    } else {
      clusters[[dogs$cluster[dogs$IID==dog2]]]<-c(clusters[[dogs$cluster[dogs$IID==dog2]]],clusters[[dogs$cluster[dogs$IID==dog1]]])
      clusters[[dogs$cluster[dogs$IID==dog1]]]<-NA
      dogs$cluster[dogs$cluster==dogs$cluster[dogs$IID==dog1]]<-dogs$cluster[dogs$IID==dog2]
    }
  }
}
sum(table(table(dogs$cluster)))
table(table(dogs$cluster))
# 
# ## filter by cluster 57: n=88
# clusters[[57]]
# cluster57rel<- subset(potRel, potRel$IID1 %in% clusters[[57]])
# cluster57Graph<-ggraphData(cluster57rel,"Dog in largest cluster")
# bob<-unique(c(cluster57rel$IID1,cluster57rel$IID2))
# 
# 
# ####################
# # Create graph of highschool friendships
# graph <- as_tbl_graph(highschool) %>% 
#   mutate(Popularity = centrality_degree(mode = 'in'))
# 
# # plot using ggraph
# ggraph(graph, layout = 'kk') + 
#   geom_edge_link(aes(col=year)) + 
#   geom_node_point()+
#   ggtitle("highschool")
# 
# 
# qgraphedgemat<-as.matrix(cbind(highschool$from,highschool$to))
# qgraph(qgraphedgemat,labels=FALSE,directed=FALSE)


#######################
# Dogs with qgraph

edgeMat<-cbind(potRel$IID1,potRel$IID2)

#qgraph(edgeMat,labels=FALSE,directed=FALSE)


## function to create graph
create_qgraph_legend<-function(dataset,title){
  edgeMat<-cbind(dataset$IID1,dataset$IID2)
  qgraph(edgeMat,layout="spring",repulsion=0.9,title=title,labels=FALSE,directed=FALSE,edge.color=ifelse(dataset$PI_HAT==1,"blue","red"),edge.width=2,node.width=0.5,mar=c(3,3,3,5))
  legend("bottomright",c("TRUE","FALSE"),col = c("blue","red"),lty=1,lwd=2,title = "Identical",bty="n")
}
create_qgraph<-function(dataset,title){
  edgeMat<-cbind(dataset$IID1,dataset$IID2)
  qgraph(edgeMat,layout="spring",repulsion=0.9,title=title,labels=FALSE,directed=FALSE,edge.color=ifelse(dataset$PI_HAT==1,"blue","red"),edge.width=2,node.width=0.5,mar=c(3,3,3,3))
}


# All potentially related dogs
#allPotRelGraph<-create_qgraph(potRel,"All Potentially Related Dogs")
jpeg("AllPotentiallyRelated.jpg",res=300,height=14,width=16,units='cm',quality=100)
create_qgraph_legend(potRel,"All Potentially Related Dogs") #+ theme_graph(plot_margin = margin(0, 0, 0, 0))
dev.off()

# Only identical
create_qgraph(identical,"Only Identical")

# Only in breed_dog
breed<-create_qgraph(breed_dog,"breed_dog")
#Only in at_risk
risk<-create_qgraph(at_risk,"at_risk")
#Only in discovery
disc<-create_qgraph(discovery,"discovery")
# Between at_risk and discovery
riskdisc<-create_qgraph(at_risk_discovery,"discovery and at_risk")
# Between breed_dog and discovery
breeddisc<-create_qgraph(breed_dog_discovery,"discovery and breed_dog")
# Between at_risk and breed_dog
riskbreed<-create_qgraph(at_risk_breed_dog,"at_risk and breed_dog")

jpeg("relatedWithinAndBetween.jpg",res=300,height=14,width=20,units='cm',quality=100)
layout(matrix(c(1,1,2,2,3,3,7,4,4,5,5,6,6,7),2,7,byrow = TRUE))
create_qgraph(discovery,"A discovery")
create_qgraph(at_risk,"B at_risk")
create_qgraph(breed_dog,"C breed_dog")
create_qgraph(at_risk_discovery,"D discovery and at_risk")
create_qgraph(breed_dog_discovery,"E discovery and breed_dog")
create_qgraph(at_risk_breed_dog,"F at_risk and breed_dog")
par(mar=c(0,0,0,0))
plot.new()
legend("center",c("TRUE","FALSE"),col = c("blue","red"),lty=1,lwd=2,title = "Identical",bty="n")
dev.off()