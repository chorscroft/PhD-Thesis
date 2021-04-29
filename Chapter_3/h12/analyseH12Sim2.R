
## Set parameters
#simNeutralFolder<-  "/home/ch19g17/msms/NeutralSimulations"
simSelectFolder<-  "/home/ch19g17/msms/SelectionSimulations2"
#simNeutralOut <- "h12NeutralOut/"
simSelectOut <- "h12SelectOut2/"

windowSize <- 400 #200 snps either side
stepSize <- 20 #snps

midsection <- c(150001,350000)

## Get R functions
source("H12Sim.R")


## Get h12 values for simulations
#h12Neutral <- h12Sim(simNeutralFolder, windowSize,stepSize,simNeutralOut)
h12Select <- h12Sim(simSelectFolder, windowSize,stepSize,simSelectOut)

## Find max values for h12 in the middle 200kb section
## for each simulation and record in matrix where
## first column is 0 if from neutral simulation and
## 1 where it is from the selection simulation
#h12MaxMatrix <-geth12MaxMatrix(h12Neutral,h12Select,midsection)
#write.table(h12MaxMatrix,file="h12MaxMatrix.txt",col.names = FALSE,row.names = FALSE)










