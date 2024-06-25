#!/usr/bin/Rscript

input <- commandArgs(T)

## For testing
# input <- c("HsInv0015", "0")

inv <- input[1]
rep <- input[2]

## Reading data

ImpGT <- read.table(paste0("Results_ScoreInvHap/", rep, "/Outputs/",inv,"/sIH_output.csv"), sep = "\t", header = T)
rownames(ImpGT) <- ImpGT$Samples


ExpGT <- read.table("../VCFs/30X/GTypesINVs.csv", sep = "\t", header = T)
rownames(ExpGT) <- ExpGT$Sample.ID
ExpGT <- ExpGT[rownames(ImpGT), c("Sample.ID","SuperPop",inv)]
# ExpGT[,inv] <- sapply(ExpGT[,inv], function(x){
#   y <- gsub("Std", "0", x)
#   y <- gsub("Inv", "1", y)
#   y
# })

## Joining
ImpGT$ImpGenos <- sapply(ImpGT$Result, function(x){
  if(x == "NN"){
    "0/0"
  } else if(x == "NI"){
    "0/1"
  } else if(x == "II"){
    "1/1"
  }
})

## VCF
dummyVCF <- matrix(c(1, 1000, "Exp", "Std", "Inv", ".", "PASS", ".", "GT",
                     1, 1001, "Imp", "Std", "Inv", ".", "PASS", ".", "GT"), 
                   byrow = T, nrow = 2)
colnames(dummyVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
AFR <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Experimental", "ImpGenos")]))
EUR <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Experimental", "ImpGenos")]))
EAS <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Experimental", "ImpGenos")]))

write.table(AFR, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/AFR.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EUR, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EUR.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EAS, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EAS.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)







