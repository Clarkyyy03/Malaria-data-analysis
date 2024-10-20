install.packages("pheatmap")
install.packages("data.table")
install.packages("readxl")
install.packages("ggvenn")
install.packages("dplyr")
install.packages("tidyr")

library(ggvenn)
library(readxl)
library(data.table)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(tidyr)


#Venn diagrams
d1 <- read_excel("~/Desktop/data analysis/ATII 24-25 RNA-Seq.xlsx")
 list = list('7dvs.24h.DOWN' = d1$`Gene ID`[d1$Regulation.7vs24== "DOWN"],
             '7dvs.1h.DOWN' = d1$`Gene ID`[d1$Regulation.7vs1== "DOWN"])
ggvenn(list) 

# dl1 = list(`24vs1.UP`=na.omit(d$`Gene ID`[d$`log2ratio:24h v 1h` > 1 & d$`p-val:24h v 1h` < 0.05]),
#            `24vs1.DOWN`=na.omit(d$`Gene ID`[d$`log2ratio:24h v 1h` < -1 & d$`p-val:24h v 1h` < 0.05]),
#            `7vs24.UP`=na.omit(d$`Gene ID`[d$`log2ratio:7d v 24h`>1 & d$`p-val:7d v 24h`<0.05]),
#            `7vs24h.DOWN`=na.omit(d$`Gene ID`[d$`log2ratio:7d v 24h`< -1 & d$`p-val:7d v 24h`<0.05]),
#            `7vs1.UP`=na.omit(d$`Gene ID`[d$`log2ratio:7d v 1h`> 1 & d$`p-val:7d v 1h`<0.05]),
#            `7vs1.DOWN`=na.omit(d$`Gene ID`[d$`log2ratio:7d v 1h`< -1 & d$`p-val:7d v 1h`<0.05]) )
# summary(dl1)

#List those overlapped genes
library(ggvenn)
gene_sets <- list(
  `7d vs 24h DOWN` = d1$`Gene ID`[d1$Regulation.7vs24== "DOWN"],
  `7d vs 1h DOWN` = d1$`Gene ID`[d1$Regulation.7vs1== "DOWN"]
)
overlap <- Reduce(intersect, gene_sets)
write.table(overlap, file = "overlapping_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#Further GO enrichment analysis
