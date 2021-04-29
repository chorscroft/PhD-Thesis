
bob<-read.table("plink.hwe",header=TRUE)

max(bob$P)
sum(bob$P>=0.001)

sum(bob$P>=0.001)/nrow(bob)
