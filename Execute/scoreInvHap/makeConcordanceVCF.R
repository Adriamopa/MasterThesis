#!/usr/bin/Rscript

library(stringr)
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
VCF_Exp <- matrix(c(1, 1000, "T1", "Std", "Inv", ".", "PASS", ".", "GT"), 
                  byrow = T, nrow = 1)
VCF_Imp <- matrix(c(1, 1000, "T1", "Std", "Inv", ".", "PASS", ".", "GT"), 
                  byrow = T, nrow = 1)

AFR_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Experimental")]))
AFR_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("ImpGenos")]))

colnames(AFR_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))
colnames(AFR_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))

EUR_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Experimental")]))
EUR_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("ImpGenos")]))

colnames(EUR_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))
colnames(EUR_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))

EAS_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Experimental")]))
EAS_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("ImpGenos")]))

colnames(EAS_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))
colnames(EAS_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Samples")]))


write.table(EUR_Exp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EUR_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EUR_Imp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EUR_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EAS_Exp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EAS_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EAS_Imp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/EAS_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(AFR_Exp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/AFR_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(AFR_Imp, paste0("Results_ScoreInvHap/", rep, "/R2/", inv, "/AFR_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)







