
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)
library(futile.logger)
library(VennDiagram)

### 1) Paths & Inputs ###
# ----------------------------------------------------------------------------------------------------------- #

#Get Orthodox list
setwd("/Users/stubbf02/Fx_Stubbe/ressources/genomes/staph/All/orthodox/fastq/contigs")
orthodox <- list.files()
orthodox <- sapply(c(1:length(orthodox)), function(x) { gsub(".final.contigs.fasta", "", orthodox[x] ) } ) 

#Get resistome of each strain
setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD/Resistome")
my_files <- list.files(pattern = ".txt")
Isolates <- sapply(c(1:length(my_files)), function(x) { gsub("._resistome.txt", "", my_files[x] ) } ) 
my_tables <- lapply(my_files, function(i) { fread( i, select=c(9,15) ) } )

### 2) Data Wrangling ###
# ----------------------------------------------------------------------------------------------------------- #

if (Isolates[x] %in% orthodox){
  for (resistance in "list de resistance") {}
}

#Initialising a dataframe containing all unique"Best_Hit_ARO" & "Drug Class"

df <- bind_rows(my_tables) 
df <- trial[!duplicated(trial$Best_Hit_ARO),]
df$nOrtho <- 0
df$nHospi <- 0


#Counting How many orthodox or non orthodox hit for a given resistance 

for( k in 1:length(Isolates) ){
  
  m <- as.data.frame(matrix(NA,nrow = nrow(my_tables[[k]]), ncol =  ncol(my_tables[[k]])))
  m <- cbind( m, my_tables[[k]])
  m <- m[,-c(1,2)]
  
  if (Isolates[k] %in% orthodox) { 
    df$nOrtho <- sapply( c( 1:length( df$Best_Hit_ARO ) ), function( x )  { 
      if( df$Best_Hit_ARO[x] %in% m$Best_Hit_ARO ) { 
        as.numeric(df$nOrtho[x]) + 1 
      } else { as.numeric(df$nOrtho[x]) + 0 }
    } )
  } else {
    df$nHospi <- sapply( c( 1:length( df$Best_Hit_ARO ) ), function( x )  { 
      if( df$Best_Hit_ARO[x] %in% m$Best_Hit_ARO ) { 
        as.numeric(df$nHospi[x]) + 1 
      } else { as.numeric(df$nHospi[x]) + 0 }
    } )
  }
}

#Adding Columns containing frequencies of the resistance in either orthodox or non orthodox

df$Ortho_freq <- sapply( c( 1:nrow(df) ), function(x) {
  ( df$nOrtho[x]/length(orthodox) )*100
} ) 

df$Hospi_freq <- sapply( c( 1:nrow(df) ), function(x) {
  ( df$nHospi[x]/(length(Isolates) - length(orthodox) ) )*100
} ) 

df$Count <- sapply( c( 1:nrow(df) ), function(x) {
  df$nOrtho[x] + df$nHospi[x]
} ) 

write.csv(df, "resistome_summary.csv")