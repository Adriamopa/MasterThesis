#!/usr/bin/Rscript

input <- commandArgs(T)

wd <- input[1]
pop <- input[2]

## For testing
#wd <- "/home/elias/Desktop/ImputationBenchmarking/scoreInvHap/1_ObtainInputs/Inputs/HsInv0072"

VCF <- read.table(paste0(wd, "/Body1"), sep = "\t", header = T, check.names = F, comment.char = "@")

### New IDs
for(i in seq(dim(VCF)[1])){
  if(VCF[i,3] == "."){
    VCF[i,3] <- paste(VCF[i,1], VCF[i,2], VCF[i,4], VCF[i,5], sep = "_")
  }
}

write.table(VCF, paste0(wd, "/Body2"), sep = "\t", quote = F, row.names = F, col.names = T)

## Only for chr X
if(mean(nchar(VCF[1,10:dim(VCF)[2]])) < 3){
  ## Sex
  Sex <- as.data.frame(colnames(VCF)[10:dim(VCF)[2]])
  colnames(Sex) <- "Sample"
  rownames(Sex) <- Sex$Sample
  Sex$Gender <- sapply(VCF[1,10:dim(VCF)[2]], function(x){
    if(nchar(x) == 3){
      "Female"
    } else {
      "Male"
    }
  })
  
  write.table(Sex, paste0(wd, "/Sex_", pop), sep = "\t", quote = F, row.names = F, col.names = T)
  
  ## Also make males diploid
  VCF[VCF == "0"] <- "0|0"
  VCF[VCF == "1"] <- "1|1"
  
  write.table(VCF, paste0(wd, "/Body2_dip"), sep = "\t", quote = F, row.names = F, col.names = T)
}



