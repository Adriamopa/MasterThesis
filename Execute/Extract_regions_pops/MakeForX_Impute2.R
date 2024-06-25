#!/usr/bin/Rscript

input <- commandArgs(T)

## For testing
# input <- c("HsInv0072", "EUR", "NA06985")

inv <- input[1]

Samples <- read.table(paste0("../VCFs/Inversion_Regions/Impute2_Phase/", inv, "/Samples"), comment.char = "@")
Samples <- as.character(Samples)

PopFile <- read.table("../VCFs/30X/Panel30x", header = T)
rownames(PopFile) <- PopFile$sample
PopFile <- PopFile[Samples,]

All <- as.data.frame(Samples)
colnames(All) <- "ID_1"
All$ID_2 <- All$ID_1
All$missing <- rep("0.0", dim(All)[1])
All$sex <- sapply(PopFile$gender, function(x){
  if(x == "male"){
    1
  } else {
    2
  }
})

Base <- as.data.frame(t(c("0", "0", "0", "D")))
colnames(Base) <- c("ID_1", "ID_2", "missing", "sex")

Ref <- rbind(Base, All)

write.table(Ref, paste0("../VCFs/Inversion_Regions/Impute2_Phase/", inv,"/Reference.sample"), quote = F, sep = "\t", row.names = F, col.names = T)
