
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)
library(futile.logger)
library(VennDiagram)

### 3) Plotting - Bubble Graph ###
# ----------------------------------------------------------------------------------------------------------- #

setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/MRSA-CARD/")
df <- fread("resistome_%mupA_.csv", header = T)

theme_set(
  theme_bw() + 
    theme(legend.position = "bottom")
)

ggplot(df, aes(x = Hospi_freq, y = Ortho_freq)) + 
  geom_point(aes(color = Drug, size = Count), alpha = 0.5) +
  scale_size(range = c(0.5, 22)) +
  geom_label_repel(aes(x = Hospi_freq, y = Ortho_freq,
                       label = ifelse(Ortho_freq > 40 & Hospi_freq < 50, as.character(Best_Hit_ARO),'') ),     
                   nudge_x      = -0.3,
                   direction    = "y",
                   hjust        = 0,
                   segment.size = 0.2) +
  labs(x = "Non-MupA (%)", y = "MupA (%)", color = "", size = "") +
  geom_smooth(method = "lm", se = F, color = "black", linetype = 2, size = 0.4)

### 3) Plotting - Bubble Graph ###
# ----------------------------------------------------------------------------------------------------------- #