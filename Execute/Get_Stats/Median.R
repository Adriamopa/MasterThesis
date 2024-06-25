## Load required libraries
#library(dplyr)
#library(tidyverse)
suppressPackageStartupMessages(library(readr))


## Get the input agruments
input <- commandArgs(T)
impute_method <- input[1] #Impute2 Impute5 Beagle ScoreInvHap Minimac4
#impute_method <- "Impute5"
pop <- input[2] # pops <- c("AFR", "EUR", "EAS")
#pop <- "EAS"
condition <- input[3]
#condition <- "Standard"
filter <- input[4] # filters <- c(0.50, 0.70, 0.80, 0.90, 0.95)
#filter <- "0.90"
wd <- getwd()
#wd <- "/home/aftimo/Documents/practiques/Documents_Elias/IntroMasterAdria/Execute"


## Get the input data
if (is.na(condition)) {
  Repetitions <- dir(paste0(wd,"/Results_", impute_method), pattern = '^Rep_[0-9]$', include.dirs = T, full.names = T)
} else { 
  Repetitions <- dir(paste0(wd,"/Results_", impute_method, "/", condition), pattern = '^Rep_[0-9]$', include.dirs = T, full.names = T)
}
num_reps <- length(Repetitions)
if (is.na(filter)) {
  final_links <- lapply(Repetitions, function(x) paste0(x,"/Imp_", pop, ".r2"))
} else {
  list_links <- lapply(Repetitions, function(x) paste0(x,"/Lists/"))
  filter_links <- paste0(list_links, rep(filter, length.out = num_reps))
  final_links <- lapply(filter_links, function(x) paste0(x,"/Imp_", pop, ".r2"))
}

for(i in 1:num_reps) {
  assign(paste0("Rep_",as.character(i)),read.csv(file = as.character(final_links[i]), sep = "\t",header = F, na.strings = c("","R2")))
}



## Join the tables to work with
full_data <- data.frame("Rep_1" = Rep_1[,2], row.names = Rep_1[,1] )
if (num_reps >= 2){
  for (i in 2:num_reps){
    full_data[paste0("Rep_",i)] <- get(paste0("Rep_",i))[,2]
  }
}


## Perform the average of the repetition of the imputation
Summary_data <- data.frame(row.names = row.names(full_data), "Mean" = rowMeans(full_data, na.rm = T), "Standard Deviation" = apply(full_data, 1, sd,na.rm = T))
Summary_data$Mean[is.nan(Summary_data$Mean)] <- NA


## Check the imputability of the inversions, prepare output data
Summary_data["Very.Imputable"] <- Summary_data$Mean >= 0.9
Summary_data["Quite.Imputable"] <- (Summary_data$Mean >= 0.8) & (Summary_data$Mean < 0.9)
Summary_data["Just.Imputable"] <- Summary_data$Mean >= 0.8
Summary_data["Barely.Imputable"] <- (Summary_data$Mean >= 0.6) & (Summary_data$Mean < 0.8)
Summary_data["Not.Imputable"] <- Summary_data$Mean < 0.6
Summary_data["Is.NA"] <- is.na(Summary_data$Mean)

Output_data <- Summary_data
Output_data$INV <- row.names(Output_data)
Output_data <- Output_data[,c(ncol(Output_data),c(1:ncol(Output_data))-1)]
row.names(Output_data) <- c()
names(Output_data) <- c(paste0("INV"),
                        paste0("Mean.",impute_method,".",condition,".",filter,".",pop),
                        paste0("SD.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Very.Imputable.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Quite.Imputable.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Just.Imputable.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Barely.Imputable.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Not.Imputable.",impute_method,".",condition,".",filter,".",pop),
                        paste0("Is.NA.",impute_method,".",condition,".",filter,".",pop))


## Get the available Results File and append new results. If not, create the Results File
tryCatch({
  Previous_data <- read.csv(paste0(wd,"/Results_Stats/FullData_Imputation.txt"), sep=",", header = T)
  Bind_data <- cbind(Previous_data,Output_data[,2:ncol(Output_data)])
  write.csv(Bind_data,paste0(wd,"/Results_Stats/FullData_Imputation.txt"), row.names = FALSE)
  print("Adding to Output File")
}, error = function(e) {
  write.csv(Output_data,paste0(wd,"/Results_Stats/FullData_Imputation.txt"), row.names = FALSE)
  print("Initializing Output File")
})


## Output the summary results in a Summary File
write_lines(paste0("IMPUTATION METHOD: ", impute_method, "\tCONDITION: ", condition,"\tGP FILTER: ", filter, "\tPOPULATION: ",pop),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
write_lines(paste0("\tVERY IMPUTABLE:\t",sum(Summary_data$Very.Imputable, na.rm = T)),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
write_lines(paste0("\tQUITE IMPUTABLE:\t",sum(Summary_data$Quite.Imputable, na.rm = T)),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
write_lines(paste0("\tBARELY IMPUTABLE:\t",sum(Summary_data$Barely.Imputable, na.rm = T)),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
write_lines(paste0("\tNOT IMPUTABLE:\t",sum(Summary_data$Not.Imputable, na.rm = T)),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
write_lines(paste0("\tNA:\t",sum(Summary_data$Is.NA, na.rm = T),"\n","-------------------------------------------------------------------------"),paste0(wd,"/Results_Stats/Summary_Stats.txt"), sep="\n", append = TRUE)
print("Adding to Summary File")


