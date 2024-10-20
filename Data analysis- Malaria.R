install.packages("pheatmap")
install.packages("data.table")
install.packages("readxl")
install.packages("ggvenn")

library(ggvenn)
library(readxl)
library(data.table)
library(tidyverse)
library(pheatmap)

#Venn diagrams
d1 <- read_excel("~/Desktop/data analysis/ATII 24-25 RNA-Seq.xlsx")
 list = list('7dvs.24h.DOWN' = d1$`Gene ID`[d1$Regulation.7vs24== "DOWN"],
             '7dvs.1h.DOWN' = d1$`Gene ID`[d1$Regulation.7vs1== "DOWN"])
ggvenn(list) 

#List those overlapped genes
library(ggvenn)
gene_sets <- list(
  `7d vs 24h DOWN` = d1$`Gene ID`[d1$Regulation.7vs24== "DOWN"],
  `7d vs 1h DOWN` = d1$`Gene ID`[d1$Regulation.7vs1== "DOWN"]
)
overlap <- Reduce(intersect, gene_sets)
write.table(overlap, file = "overlapping_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Further GO enrichment analysis
getwd()
go.enrich = read.csv("/Users/clarkeyin/Desktop/data analysis/P.falciparum/7h.24,1h down.tsv",  sep = "\t", header = TRUE)

ggplot(go.enrich, aes(x = -log10(adjusted_p_value), y = term_name, size = interaction_size, color = highlighted)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "GO Enrichment Analysis",
    x = "-log10(Adjusted P-Value)",
    y = "GO Terms",
    size = "Interaction Size",
    color = "Highlighted"
  )
a = go.enrich$
tmp = go.enrich %>% filter(hightlighted== TRUE) %>% separate_rows("interactions", seq = ",")
genelist = unique()
