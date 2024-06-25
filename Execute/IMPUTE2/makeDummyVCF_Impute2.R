#!/usr/bin/Rscript

input <- commandArgs(T)

## For testing
# input <- c("HsInv0114", "0.90")

inv <- input[1]
filter <- as.numeric(input[2])
fChar <- input[2]
condition <- input[3]
rep <- input[4]

## Reading data

ImpGT <- read.table(paste0("Results_Impute2/",condition,"/", rep, "/Inversions/", inv,"/ImpGT"), sep = "\t", header = F)
colnames(ImpGT) <- c("Sample", "SS", "SI", "II")
rownames(ImpGT) <- ImpGT$Sample
##Unphase
# ImpGT$GT <- sapply(ImpGT$GT, function(x){
#   y = gsub("\\|", "/", x)
#   y = gsub("1/0", "0/1", y)
#   y
# })


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
#### Filter measeres the deviation from exact values 0, 1 and 2
#### For example, filter = 0.2 for an genotype 0/1 accept the value if:
####  probability is ranging 0.6-1.4
#### For filter = 0.8 the range must be 0.9-1.1 for heterozygotes

ImpGT$Max <- sapply(ImpGT$Sample, function(x){
  max(ImpGT[ImpGT$Sample == x, c("SS", "SI", "II")])
})
ImpGT$Pass <- sapply(ImpGT$Max, function(x){
  if(x > filter){
    TRUE
  } else {
    FALSE
  }
})
ImpGT <- ImpGT[ImpGT$Pass,]

ImpGT$Imp <- sapply(ImpGT$Sample, function(x){
  posibs <- ImpGT[x, c("SS", "SI", "II")]
  c("0/0", "0/1", "1/1")[which.max(posibs)]
})

## VCF
dummyVCF <- matrix(c(1, 1000, "Exp", "Std", "Inv", ".", "PASS", ".", "GT",
                     1, 1001, "Imp", "Std", "Inv", ".", "PASS", ".", "GT"), 
                   byrow = T, nrow = 2)
colnames(dummyVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

AFR <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "AFR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Exp", "Imp")]))
EUR <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EUR","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Exp", "Imp")]))
EAS <- cbind(dummyVCF, t(ImpGT[ExpGT[intersect(ExpGT[ExpGT$SuperPop == "EAS","Sample.ID"], rownames(ImpGT)),
                                     "Sample.ID"],c("Exp", "Imp")]))

write.table(AFR, paste0("Results_Impute2/",condition,"/",rep,"/R2/", inv, "/", fChar, "/AFR.vcf"), sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = F)
write.table(EUR, paste0("Results_Impute2/",condition,"/",rep,"/R2/", inv, "/", fChar, "/EUR.vcf"), sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = F)
write.table(EAS, paste0("Results_Impute2/",condition,"/",rep,"/R2/", inv, "/", fChar, "/EAS.vcf"), sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = F)






