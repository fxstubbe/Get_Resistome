###
# Get the resistome of each strain using CARD
###

### 1) Load packages###
# ----------------------------------------------------------------------------------------------------------- #

suppressPackageStartupMessages(library(data.table))

### 1) paths and inputs ###
# ----------------------------------------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)
f.root <- getwd()
f.contigs <- as.character(args[1])

#Get contigs ID
setwd(f.contigs)
my_files <- list.files(pattern = ".final.contigs.fasta")


#Getting the resistome

setwd(f.root)
f.resistome <- "find_resistome.sh"

for (i in 1:length( my_files ) ){
  cat( paste( "echo   Getting resistome of : ", gsub(".final.contigs.fasta", "",my_files[i]), sep = "" ), sep = "\n", file = f.resistome, append = T )
  cat( paste( "echo Task : ", as.character(i), "/", as.character(length(my_files)), sep = ""), sep = "\n", file = f.resistome, append = T )
  cat( paste( "rgi main",
              " --input_sequence /Users/stubbf02/Fx_Stubbe/ressources/genomes/staph/MGE_project/CC8/contigs/", my_files[i],
              " --output_file ./Resistome/", gsub("final.contigs.fasta", "", my_files[i]), "_resistome",
              " --input_type contig",
              " --clean",
              " --local", sep = ""), sep = "\n" , file = f.resistome, append = T )
}

cat(paste("echo ... DONE with Resistome Finding", sep = ""), sep = "\n", file = f.resistome, append = T )
