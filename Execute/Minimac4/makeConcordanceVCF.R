#!/usr/bin/Rscript

library(stringr)

input <- commandArgs(T)

## For testing
# input <- c("HsInv0015", "0.50")

inv <- input[1]
filter <- as.numeric(input[2])
fChar <- input[2]
condition <- input[3]
rep <- input[4]

## Reading data

ImpGT <- read.table(paste0("Results_Minimac4/",condition, "/", rep ,"/Inversions/",inv,"/ImpGT"), sep = "\t", header = F, fill = T)
colnames(ImpGT) <- c("Sample", "GT", "GP")
rownames(ImpGT) <- ImpGT$Sample
##Unphase
ImpGT$GT <- sapply(ImpGT$GT, function(x){
  y = gsub("\\|", "/", x)
  y = gsub("1/0", "0/1", y)
  y
})


ExpGT <- read.table("../VCFs/30X/GTypesINVs.csv", sep = "\t", header = T)
rownames(ExpGT) <- ExpGT$Sample.ID
ExpGT <- ExpGT[rownames(ImpGT), c("Sample.ID","SuperPop",inv)]
ExpGT[,inv] <- sapply(ExpGT[,inv], function(x){
  y <- gsub("Std", "0", x)
  y <- gsub("Inv", "1", y)
  y
})

## Joining
ImpGT$Exp <- ExpGT[rownames(ImpGT), inv]

## Filtering
#### The filtering goes from 0 to 1, being 0 no filter and 1 exact filter
#### Filter measures the deviation from exact values 0, 1 and 2
#### For example, filter = 0.2 for an genotype 0/1 accept the value if:
####  probability is ranging 0.6-1.4
#### For filter = 0.8 the range must be 0.9-1.1 for heterozygotes

ImpGT$Filt <- sapply(ImpGT$Sample, function(x){
  
  if(ImpGT[x, "GP"] == ""){
    if(ImpGT[x, "DS"] < 0.5){
      1 - ImpGT[x, "DS"]
    } else {
      ImpGT[x, "DS"]
    }
    
  } else {
    max(as.numeric(str_split(string = ImpGT[x, "GP"], pattern = ",", simplify = T)))
  }
  
})
ImpGT$Pass <- sapply(ImpGT$Filt, function(x){
  if(x >= filter){
    T
  } else {
    F
  }
})

ImpGT <- ImpGT[ImpGT$Pass,]


## VCF
VCF_Exp <- matrix(c(1, 1000, "T1", "Std", "Inv", ".", "PASS", ".", "GT"), 
                  byrow = T, nrow = 1)
VCF_Imp <- matrix(c(1, 1000, "T1", "Std", "Inv", ".", "PASS", ".", "GT"), 
                  byrow = T, nrow = 1)

AFR_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Exp")]))
AFR_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("GT")]))

colnames(AFR_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))
colnames(AFR_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))

EUR_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Exp")]))
EUR_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("GT")]))

colnames(EUR_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))
colnames(EUR_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))

EAS_Exp <- cbind(VCF_Exp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("Exp")]))
EAS_Imp <- cbind(VCF_Imp, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                        "Sample.ID"],c("GT")]))

colnames(EAS_Exp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))
colnames(EAS_Imp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Sample")]))


write.table(EUR_Exp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/EUR_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EUR_Imp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/EUR_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EAS_Exp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/EAS_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(EAS_Imp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/EAS_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(AFR_Exp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/AFR_Exp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)
write.table(AFR_Imp, paste0("Results_Minimac4/",condition, "/", rep,"/R2/", inv, "/", fChar, "/AFR_Imp.vcf"), sep = "\t", row.names = F, col.names = T,
            quote = F)







