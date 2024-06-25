#!/usr/bin/Rscript

suppressPackageStartupMessages(library(stringr))
input <- commandArgs(T)

impute_method <- input[1] #Impute2 Impute5 Beagle Minimac4
#impute_method <- "Impute5"
pop <- input[2] # pops <- c("AFR", "EUR", "EAS")
#pop <- "EAS"
condition <- input[3]
#condition <- "Standard"
filter <- input[4] # filters <- c(0.50, 0.70, 0.80, 0.90, 0.95)
#filter <- "0.90"
filt <- as.numeric(filter)
wd <- getwd()
#wd <- "/home/aftimo/Documents/practiques/Documents_Elias/IntroMasterAdria/Execute"
invs <- readLines(paste0(wd,"/../VCFs/30X/Invs2Imp"))

## Declare Useful functions
get_num_indivs <- function(ImpGT) {
  
  if (impute_method == "Impute2") {
    colnames(ImpGT) <- c("Sample", "SS", "SI", "II")
    ImpGT$Filt <- sapply(ImpGT$Sample, function(x){
      max(ImpGT[ImpGT$Sample == x, c("SS", "SI", "II")])
    })
    
  } else if (impute_method == "Minimac4") {
    colnames(ImpGT) <- c("Sample", "GT", "GP")
    ImpGT$Filt <- sapply(ImpGT$GP, function(x){
      max(as.numeric(str_split(string = x, pattern = ",", simplify = TRUE)))
    })
    
  } else {
    colnames(ImpGT) <- c("Sample", "GT", "DS", "GP")
    ImpGT$Filt <- sapply(ImpGT$GP, function(x){
      max(as.numeric(str_split(string = x, pattern = ",", simplify = TRUE)))
    })
  }
  
  ImpGT$Pass <- sapply(ImpGT$Filt, function(x){
    if(x >= filt){
      TRUE
    } else {
      FALSE
    }
  })
  num_indiv <- sum(ImpGT$Pass)
  return(num_indiv)
}


## Get the input data
if (is.na(condition)) {
  Repetitions <- dir(paste0(wd,"/Results_", impute_method), pattern = '^Rep_[0-9]$', include.dirs = T, full.names = T)
} else { 
  Repetitions <- dir(paste0(wd,"/Results_", impute_method, "/", condition), pattern = '^Rep_[0-9]$', include.dirs = T, full.names = T)
}
num_reps <- length(Repetitions)
for (i in 1:num_reps) {
  assign(paste0("Rep_",as.character(i)), dir(paste0(Repetitions[i],"/Inversions"), include.dirs = T, full.names = T))
  assign(paste0("Rep_",as.character(i)),lapply(get(paste0("Rep_",as.character(i))), function(x) paste0(x,"/ImpGT_", pop)))
}

for(i in 1:num_reps) {
  current <- paste0("Rep_",as.character(i))
  for (j in 1:length(invs)) {
    tryCatch({
      assign(paste0(current,".",invs[j]),read.csv(file = as.character(get(current)[j]), sep = "\t",header = F))
    }, error = function(e) {
      print(paste0("File for imp_method: ", impute_method, ", condition: ", condition, ", inv: ", invs[j], " population: ", pop," filter: ", filter, " and rep: ", current," is empty!!"))
      assign(paste0(current,".",invs[j]), "")
    })
  }
}



## Do the work
Output_data <- data.frame()
for (i in 1:length(invs)) {
  rep_indivs <- c()
  for(j in 1:num_reps) {
    current <- paste0("Rep_",as.character(j))
    assign(paste0(current,".",invs[i]),read.csv(file = as.character(get(current)[i]), sep = "\t",header = F))
    rep_indivs <- c(rep_indivs,get_num_indivs(get(paste0(current,".",invs[i]))))
  }
  final_num <- round(mean(rep_indivs))
  Output_data <- rbind(Output_data, c(invs[i], final_num))
}
names(Output_data) <- c("INV", paste0("Indvs.",impute_method,".",condition,".",filter,".",pop))


## Save the Data
tryCatch({
  Previous_data <- read.csv(paste0(wd,"/Results_Stats/Indiv_per_Filt.txt"), sep=",", header = T)
  Bind_data <- cbind(Previous_data,as.double(Output_data[,2:ncol(Output_data)]))
  names(Bind_data)[length(Bind_data)] <- paste0("Indvs.",impute_method,".",condition,".",filter,".",pop)
  write.csv(Bind_data,paste0(wd,"/Results_Stats/Indiv_per_Filt.txt"), row.names = FALSE)
  print("Adding to Output File")
}, error = function(e) {
  write.csv(Output_data,paste0(wd,"/Results_Stats/Indiv_per_Filt.txt"), row.names = FALSE)
  print("Initializing Output File")
})
