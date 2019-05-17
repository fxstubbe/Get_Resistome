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
   trans <- str_c(c(trans, "._resistome.txt"), collapse = "")
   return(trans) } )
matching_ID <- matching_ID [! matching_ID %in% 'USA300_FPR3757._resistome.txt' ]

#Get resistome of each strain
setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD/Resistome")
Isolates <- sapply(c(1:length(matching_ID)), function(x) { gsub("._resistome.txt", "", matching_ID[x] ) } ) 
my_tables <- lapply(matching_ID, function(i) { fread( i, select=c(9,15) ) } )


### 2) Data Wrangling ###
# ----------------------------------------------------------------------------------------------------------- #

if (Isolates[x] %in% orthodox){
  for (resistance in "list de resistance") {}
}

#Initialising a dataframe containing all unique"Best_Hit_ARO" & "Drug Class"

df <- bind_rows(my_tables) 
df <- df[!duplicated(df$Best_Hit_ARO),]
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
  ( df$nOrtho[x]/99 )*100
} ) 

df$Hospi_freq <- sapply( c( 1:nrow(df) ), function(x) {
  ( df$nHospi[x]/(length(Isolates) - 99 ) )*100
} ) 

df$Count <- sapply( c( 1:nrow(df) ), function(x) {
  df$nOrtho[x] + df$nHospi[x]
} ) 

setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD")
write.csv(df, "resistome_mupA_resistant_summary.csv")
