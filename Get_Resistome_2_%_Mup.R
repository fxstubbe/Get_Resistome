#Make a summary file of resistomes BUT only from strains matching in MupA

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

#Get isolates matching in mupirocin

setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC8/CC8/Tables/Summary")
matching <- fread("iTOL_feat_70_full_match_independent_sig.txt", header = F )
colnames(matching) <- c("Isolate", "plasmid", "mupA", "qacA")
matching <- matching[which(matching$mupA == 1),]
matching_ID <- matching$Isolate
matching_ID <- sapply( c(1:length(matching_ID) ), function(x){
  trans <- gsub(" ", "_", matching_ID[x])
  #trans <- str_c(c(trans, "._resistome.txt"), collapse = "")
  return(trans) } )
matching_ID <- matching_ID [! matching_ID %in% 'USA300_FPR3757' ]

#Get resistome of each strain
setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD/Resistome")
my_files <- list.files(pattern = ".txt")
Isolates <- sapply(c(1:length(my_files)), function(x) { gsub("._resistome.txt", "", my_files[x] ) } ) 
my_tables <- lapply(my_files, function(i) { fread( i, select=c(9,15) ) } )


### 2) Data Wrangling ###
# ----------------------------------------------------------------------------------------------------------- #
for (i in 1:length(Isolates)){
  if (Isolates[i] %in% matching_ID){print("yep")} else {print("nop")} 
}


#Initialising a dataframe containing all unique"Best_Hit_ARO" & "Drug Class"

df <- bind_rows(my_tables) 
df <- df[!duplicated(df$Best_Hit_ARO),]
df$nMupA <- 0
df$nNotMupA <- 0


#Counting How many orthodox or non orthodox hit for a given resistance 

for( k in 1:length(Isolates) ){
  
  m <- as.data.frame(matrix(NA,nrow = nrow(my_tables[[k]]), ncol =  ncol(my_tables[[k]])))
  m <- cbind( m, my_tables[[k]])
  m <- m[,-c(1,2)]
  
  if (Isolates[k] %in% matching_ID) { 
    df$nMupA <- sapply( c( 1:length( df$Best_Hit_ARO ) ), function( x )  { 
      if( df$Best_Hit_ARO[x] %in% m$Best_Hit_ARO ) { 
        as.numeric(df$nMupA[x]) + 1 
      } else { as.numeric(df$nMupA[x]) + 0 }
    } )
  } else {
    df$nNotMupA <- sapply( c( 1:length( df$Best_Hit_ARO ) ), function( x )  { 
      if( df$Best_Hit_ARO[x] %in% m$Best_Hit_ARO ) { 
        as.numeric(df$nNotMupA[x]) + 1 
      } else { as.numeric(df$nNotMupA[x]) + 0 }
    } )
  }
}

#Adding Columns containing frequencies of the resistance in either orthodox or non orthodox

df$Ortho_freq <- sapply( c( 1:nrow(df) ), function(x) {
  ( df$nMupA[x] / length(matching_ID) )*100
} ) 

df$Hospi_freq <- sapply( c( 1:nrow(df) ), function(x) {
  ( df$nNotMupA[x] / (length(Isolates) - length(matching_ID) ) )*100
} ) 

df$Count <- sapply( c( 1:nrow(df) ), function(x) {
  df$nMupA[x] + df$nNotMupA[x]
} ) 

setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD")
write.csv(df, "resistome_%mupA_.csv")
