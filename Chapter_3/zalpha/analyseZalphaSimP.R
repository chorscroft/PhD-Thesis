
## Set parameters
#simNeutralFolder<-  "/home/ch19g17/msms/NeutralSimulations"
simSelectFolder<-  "/home/ch19g17/msms/SelectionSimulationsP"
#outputNeutral <- "zalphaNeutral/"
outputSelected <- "zalphaSelectedP/"
windowSize <- 200 #100kb either side
stepSize <- 1 #snps
midsection <- c(150001,350000) #For analysis

## Get R functions
source("r2Matrix.R")
source("zalphaSimNew.R")


## Get zalpha values for simulations
#zalphaNeutral <- zalphaSim(simNeutralFolder, windowSize, stepSize, outputNeutral)
zalphaSelect <- zalphaSim(simSelectFolder, windowSize, stepSize, outputSelected)

## Find max values for zalpha in the middle 200kb section
## for each simulation and record in matrix where
## first column is 0 if from neutral simulation and
## 1 where it is from the selection simulation
#zalphaMaxMatrix <-getZalphaMaxMatrix(zalphaNeutral,zalphaSelect,midsection)
#write.table(zalphaMaxMatrix,file="zalphaMaxMatrix.txt",col.names = FALSE,row.names = FALSE)






