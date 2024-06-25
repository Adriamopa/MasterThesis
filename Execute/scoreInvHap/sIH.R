#!/usr/bin/Rscript
suppressMessages(library(scoreInvHap))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(snpStats))
opt <- options(show.error.messages = FALSE)
on.exit(options(opt))

input <- commandArgs(T)

inv <- input[1]
pop <- input[2]
rep <- input[3]

## For testing
##inv <- "HsInv0072"
##pop <- "EUR"

#R2input
R2input <- read.table(paste0("Results_ScoreInvHap/Inputs/", inv,"/R2_", pop, ".input"), sep = "\t", header = T)
if(dim(R2input)[1]==0){
  write.table(data.frame("Samples" = c(""), "Result" = c(""), "ImpGenos" = c(""), "Experimental" = c("")), paste0("Results_ScoreInvHap/",rep, "/Outputs/", inv,"/sIH_output_", pop, ".csv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  print(paste0("Monomorphic variant for pop: ", pop, " and inversion: ", inv))
  stop()
}
R2input <- R2input[R2input$SNP_B != inv,]

dups <- R2input$SNP_B[duplicated(R2input$SNP_B)]
'%!in%' <- Negate('%in%')
R2aux <- R2input[R2input$SNP_B %!in% dups,]

for(i in dups){
  ss <- R2input[R2input$SNP_B == i,]
  R2aux <- rbind(R2aux, ss[which.max(ss$R2),])
}

R2input <- R2aux
rm(R2aux)

R2 <- as.numeric(R2input$R2)
attr(R2, "names") <- R2input$SNP_B

#VCF

vcf <- read.table(paste0("Results_ScoreInvHap/Inputs/", inv, "/Region_", pop, ".csv"), sep = "\t", header = T,
                  comment.char = "@", check.names = F)
InvLine <- vcf[vcf$ID == inv,]

if(mean(nchar(InvLine[10:dim(InvLine)[2]])) < 3){
  modeX = T
} else {
  modeX = F
}

vcf_woinv <- vcf[!grepl("HsInv*", vcf$ID),]
vcf_woinv <- vcf_woinv[vcf_woinv$ID %in% R2input$SNP_B,]

InvGen <- InvLine[,10:dim(InvLine)[2]]
InvGen[InvGen == "0"] <- 4
InvGen[InvGen == "1"] <- 3
InvGen[InvGen == "0/0"] <- 1
InvGen[InvGen == "0/1"] <- 2
InvGen[InvGen == "1/1"] <- 3
InvGen[InvGen == 4] <- 1

  
#hRef

hRef <- paste0(vcf_woinv$REF, vcf_woinv$ALT)
attr(hRef, "names") <- vcf_woinv$ID
hRef[nchar(hRef) !=2] <- "AT"

# InvGen <- vcf[grepl("HsInv*", vcf$ID),10:dim(vcf)[2]]



RefInput <- list()

from <- c("0|0", "1|0", "0|1", "1|1", "0", "1")
to <- c(1,2,2,3,1,3)

combos <- c()

for(i in seq(length(vcf_woinv$ID))){
  rs <- vcf_woinv$ID[i]
  comb_snp <- c(paste(vcf_woinv$REF[i], vcf_woinv$REF[i], sep = ""),
                paste(vcf_woinv$REF[i], vcf_woinv$ALT[i], sep = ""),
                paste(vcf_woinv$ALT[i], vcf_woinv$ALT[i], sep = "")
  )
  if(nchar(comb_snp[1]) != 2 || nchar(comb_snp[2]) != 2 || nchar(comb_snp[3]) != 2){
    comb_snp <- c("AA", "AT", "TT")
    
  }
  comb_haplo <- c("NN", "NI", "II")
  
  snp_gen <- as.character(vcf_woinv[i,10:dim(vcf_woinv)[2]])
  # snp_gen[snp_gen == "0|0"] <- 1
  # snp_gen[snp_gen == "0|1"] <- 2
  # snp_gen[snp_gen == "1|0"] <- 2
  # snp_gen[snp_gen == "1|1"] <- 3
  snp_gen <- as.numeric(unlist(sapply(snp_gen, function(gt){
    to[grepl(gt, from, fixed = T)]
  })))
  
  comp <- t(rbind(InvGen, snp_gen))
  cm <- as.matrix(table(as.data.frame(comp)))

  ocm <- matrix(rep(0,9), nrow = 3, ncol = 3)

  for(i in seq(dim(cm)[1])){
    for(j in seq(dim(cm)[2])){
      f <- as.integer(rownames(cm)[i])
      c <- as.integer(rownames(cm)[j])
      
      ocm[f,c] <- cm[i,j]
    }
  }

  cm <- as.table(as.matrix(ocm))

  attr(cm, "dimnames") <- list(haplos = comb_haplo, comb_snp)

  for(r in seq(dim(cm)[1])){
    cm[r,] <- cm[r,]/sum(cm[r,])
  }

  RefInput[[rs]] <- cm
  
}

bed <- read.plink(paste0("Results_ScoreInvHap/Inputs/", inv, "/Region_", pop, ".bed"))
bedmap <- bed[["map"]]
allele_join <- paste(bedmap$allele.1, bedmap$allele.2, sep = "")
a1 <- bedmap$allele.1
a2 <- bedmap$allele.2

for(i in seq(length(allele_join))){
  if(nchar(allele_join[i]) > 2){
    a1[i] <- "A"
    a2[i] <- "T"
  }
}

bedmap$allele.1 <- a1
bedmap$allele.2 <- a2

bed[["map"]] <- bedmap


# vcfInput <- readVcf("data/HsInv0501.vcf", "hg19")
results <- scoreInvHap(bed, SNPsR2 = R2, hetRefs = hRef, Refs = RefInput, 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 16))
# results <- scoreInvHap(bed, inv = "inv8_001")

out <- as.data.frame(names(results@classification))
colnames(out) <- "Samples"
rownames(out) <- out$Samples
out$Result <- as.character(results@classification)

if(modeX){
  Sex <- read.table(paste0("Results_ScoreInvHap/Inputs/", inv,"/Sex_", pop), sep = "\t", header = T)
  rownames(Sex) <- Sex$Sample
  
  out$ImpGenos <- sapply(out$Samples, function(x){
    gender <- Sex[x,"Gender"]
    
    if(gender == "Male"){
      from <- c("NN", "NI", "II")
      to <- c("0", "Unk", "1")
      gt <- out[x, "Result"]
      
      to[grepl(gt, from, fixed = T)]
    } else {
      from <- c("NN", "NI", "II")
      to <- c("0/0", "0/1", "1/1")
      gt <- out[x, "Result"]
      
      to[grepl(gt, from, fixed = T)]
    }
    
  })
} else {
  out$ImpGenos <- out$Result
  out[out$ImpGenos == "NN", "ImpGenos"] <- "0/0"
  out[out$ImpGenos == "NI", "ImpGenos"] <- "0/1"
  out[out$ImpGenos == "II", "ImpGenos"] <- "1/1"
}

out$Experimental <- as.character(InvLine[,rownames(out)])

write.table(out, paste0("Results_ScoreInvHap/",rep, "/Outputs/", inv,"/sIH_output_", pop, ".csv"),
            quote = F, sep = "\t", row.names = F, col.names = T)




