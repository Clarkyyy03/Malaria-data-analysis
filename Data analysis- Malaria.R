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
getwd()
library(ggplot2)
library(data.table)
go.enrich = fread("~/Desktop/data analysis/P.falciparum/gProfiler_pfalciparum_intersections.csv")
go.enrich$term_name <- factor(go.enrich$term_name, 
                              levels = c(
                                setdiff(go.enrich$term_name, c("adhesion of symbiont to microvasculature", 
                                                               "cell adhesion molecule binding", 
                                                               "host cell cytoplasm", 
                                                               "membrane")), 
                                "adhesion of symbiont to microvasculature", 
                                "cell adhesion molecule binding", 
                                "host cell cytoplasm", 
                                "membrane"))
ggplot(go.enrich, aes(x = -log10(adjusted_p_value), y = term_name, size = intersection_size, color = highlighted)) + 
  geom_point() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +  # Set colors for TRUE and FALSE
  theme_bw() +
  labs(
    title = "GO Enrichment Analysis",
    x = "-log10(Adjusted P-value)",
    y = "GO Terms",
    size = "Intersection Size",
    color = "Highlighted"
  )


## Heatmap - 105 overlapped genes
tmp = go.enrich %>% filter(highlighted==TRUE) %>% separate_rows("intersections", sep = ",")
genelist = unique(tmp$intersections)

mat1 = d1 %>% filter(`Gene ID` %in% genelist) %>% column_to_rownames(var="Gene ID") %>% select("percent abundance 1h","percent abundance 24h","percent abundance 7d")
pheatmap(mat1, scale = "row", cluster_cols = F)
