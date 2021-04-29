##Initialises a zAlpha run

tPed <- "wellderly012.tped"
windowSize <- 200 #100kb either side
source("readInTpedFile.R")
source("removeRandomSNPs.R")
source("create3By3Matrix012.R")
source("HillAlgorithm.R")
source("r2Matrix.R")
source("zalpha.R")


myZalpha <- zalpha(tPed,windowSize)
write(t(myZalpha),file="zalphaoutput4425.txt",ncolumns=2)

