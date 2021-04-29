## This code takes a VCF file and converts it into
## the two files needed to run LDHat:
## 1) ".locs" which is the position of the SNPs
##     in relation to the window size and 
## 2) ".sites" which is the file containing the data by 
##    sample, where 0 and 1 are homozygotes and 2 is a
##    heterozygote (? if the data is unknown)
##
## The inputs for the functions are assigned at the
## bottom of this code.
## 1) Change the working directory
## 2) give the name of the vcfFile
## 3) give the prefix to use for the output files
## 4) convert Bp to Kb (TRUE or FALSE)

#################################################

## Find the header line in a VCF File
findHeaderLine = function(vcfName) {
  con = file(vcfName, "r")
  headerLine <- 0
  while (TRUE){
    headerLine <- headerLine +1
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      print.default("No header line in vcf file")
      headerLine <- 0
      break
    } else if (substr(line,1,6)=="#CHROM") {
      break
    }
  }
  close(con)
  return(headerLine)
}

## Get .site and .loc files
vcfToLDHatInput <- function(vcfName,prefix="vcfdata",kbFlag=FALSE){
  
  ## Get the header row number
  headerLine <- findHeaderLine(vcfName)
  
  ## stop code if no header line
  if (headerLine==0){
    quit()
  }
  
  ## read in data table
  dataW<-read.delim(vcfName,skip=headerLine-1)
  
  ## Keep only unique lines
  dataW<-dataW[!duplicated(dataW$POS),]
  
  ## calculates window size and 
  ## changes snp positions so they are 
  ## relative to the window
  
  windowSize <- dataW$POS[nrow(dataW)]-dataW$POS[1]+1
  locations <- dataW$POS - dataW$POS[1]+1
  
  ## If required, convert Bp to Kb
  if (kbFlag == TRUE){
    windowSize <- format(windowSize/1000,digits = 8)
    locations <- format(locations/1000,digits=8)
  }
  
  ## Creates Locs file
  write(c(nrow(dataW),as.character(windowSize),"L"),paste0(prefix,".locs"),ncolumns=3,sep="\t")
  write(locations,paste0(prefix,".locs"),ncolumns=1,append =TRUE)
  
  missingvar<-NULL
  ## Creates Sites file
  write(c(ncol(dataW)-9,nrow(dataW),2),paste0(prefix,".sites"),ncolumns = 3,sep="\t")
  for (i in 10:ncol(dataW)){
    write(paste0(">",colnames(dataW)[i]),paste0(prefix,".sites"),ncolumns = 1,append=TRUE)
    genotypeData <- NULL
    for (j in 1:nrow(dataW)){
      var01 <- substr(dataW[j,i],1,3)
      if(var01 == "0/0" || var01 == "0|0"){
        genotypeData<-paste0(genotypeData,"0")
      } else if (var01 == "0/1" || var01 == "1/0" || var01 == "0|1" || var01 == "1|0"){
        genotypeData<-paste0(genotypeData,"2")
      } else if (var01 == "1/1" || var01 == "1|1"){
        genotypeData<-paste0(genotypeData,"1")
      } else {
        genotypeData<-paste0(genotypeData,"?")
        missingvar<-c(missingvar,var01)
      }
    }
    write(genotypeData,paste0(prefix,".sites"),ncolumns = 1,append=TRUE)
  }
  return(TRUE)
  write(missingvar,"missingvar.txt",ncolumns = 1)
}

#####################################################

## Set working directory
setwd("/temp/hgig/EXOME_DATA/Clare/LDhat/Baganda")

## Name of the VCF file
vcfName <- "chr22_Baganda.vcf"

## Prefix for output files
prefix <- "Baganda22"

## Runs the code
vcfToLDHatInput(vcfName,prefix,TRUE)


