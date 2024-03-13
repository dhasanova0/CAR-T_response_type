#Author: Dzhansu Hasanova
#Comparison of Stator results and scRNAseq data analysis with standard methods

set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(enrichR)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import data
DE_standard <- read.csv(paste0(wd, "figures/baseline/analysis/DE_RvsNR_log.csv"))
DE_stator <- read.csv(paste0(wd, "figures/stator_run1/thesis/valid/DE_RvsNR_entire_valid.csv"))

#Check shared genes
DE_standard <- DE_standard[order(DE_standard$avg_log2FC),]
DE_stator <- DE_stator[order(DE_stator$avg_log2FC),]

stator_up_R <- rev(DE_stator$X[(nrow(DE_stator)-19):nrow(DE_stator)])
DE_up_R <- rev(DE_standard$X[(nrow(DE_standard)-19):nrow(DE_standard)])

stator_up_NR <- DE_stator$X[1:20]
DE_up_NR <- DE_standard$X[1:20]

intersect(stator_up_NR, DE_up_NR)
intersect(stator_up_R, DE_up_R)

#Compare FC
rownames(DE_standard) <- DE_standard$X
rownames(DE_stator) <- DE_stator$X

DE_20 <- DE_standard[c(intersect(stator_up_NR, DE_up_NR), intersect(stator_up_R, DE_up_R)),]
stator_20 <- DE_stator[c(intersect(stator_up_NR, DE_up_NR), intersect(stator_up_R, DE_up_R)),]
stator_20 <- stator_20[, -(7:8)]

total <- data.frame(cbind(DE_20, stator_20, by = "X"))
total <- select(total, avg_log2FC, avg_log2FC.1)

total$avg_log2FC <- abs(total$avg_log2FC)
total$avg_log2FC.1 <- abs(total$avg_log2FC.1)

colnames(total) <- c("DE_cells", "DE_States")

total_long <- total %>% 
  pivot_longer(
    cols = `DE_cells`:`DE_States`, 
    names_to = "stator_seurat",
    values_to = "avg_log2FC"
  )


ggplot(total_long, aes(y=avg_log2FC, x=stator_seurat, fill = stator_seurat)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.05), aes(color = stator_seurat)) +
  scale_fill_manual(name = "", values=c("deepskyblue4", "#E69F00"), labels = c("DE Cells", "DE States")) +
  scale_color_manual(name = "", values=c("deepskyblue4", "#E69F00"), labels = c("DE Cells", "DE States")) +
  theme_bw()  +
  xlab("") +
  ylab("average log2 FC") + stat_compare_means() +
  scale_x_discrete(labels = c("DE Cells", "DE States")) 
ggsave("violin_comparison.png", path =  paste0(wd, "data/stator_results/run3/results/comp/"), width = 5, height = 5)


#Perform GO for the top 40 DE genes

DE_standard <- DE_standard[order(DE_standard$avg_log2FC),]
DE_stator <- DE_stator[order(DE_stator$avg_log2FC),]


dbs <- c("GO_Biological_Process_2023")

enriched_standard_NR <- enrichr(c(DE_standard$X[1:20]), dbs)
plotEnrich(enriched_standard_NR[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", 
           title = "Functional Enrichment of Biological Processes in NR")+
  theme(axis.text=element_text(size=18), plot.title = element_text(face="bold", size = 20), axis.title= element_text(size = 18),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=11))+ ylab("Gene count") 
ggsave("GO_standard_NR.png", path = paste0(wd, "data/stator_results/run3/results/comp/"), width = 15, height = 5)

enriched_stator_NR <- enrichr(c(DE_stator$X[1:20]), dbs)
plotEnrich(enriched_stator_NR[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", 
           title = "Functional Enrichment of Biological Processes in NR states")+
  theme(axis.text=element_text(size=18), plot.title = element_text(face="bold", size = 20), axis.title= element_text(size = 18),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=11))+ ylab("") + ylab("Gene count")
ggsave("GO_stator_NR.png", path = paste0(wd, "data/stator_results/run3/results/comp/"), width = 18, height = 5)

enriched_standard_R <- enrichr(rev(c(DE_standard$X[(nrow(DE_standard)-19):nrow(DE_standard)])), dbs)
plotEnrich(enriched_standard_R[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", 
           title = "Functional Enrichment of Biological Processes in R")+
  theme(axis.text=element_text(size=18), plot.title = element_text(face="bold", size = 20), axis.title= element_text(size = 18),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=11))+ ylab("") + ylab("Gene count")
ggsave("GO_standard_R.png", path = paste0(wd, "data/stator_results/run3/results/comp/"), width = 19, height = 5)

enriched_stator_R <- enrichr(rev(c(DE_stator$X[(nrow(DE_stator)-19):nrow(DE_stator)])), dbs)
plotEnrich(enriched_stator_R[[1]], showTerms = 10, numChar = 100, y = "Count", orderBy = "P.value", 
           title = "Functional Enrichment of Biological Processes in R states")+
  theme(axis.text=element_text(size=18), plot.title = element_text(face="bold", size = 20), axis.title= element_text(size = 18),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=11))+ ylab("") + ylab("Gene count")
ggsave("GO_stator_R.png", path = paste0(wd, "data/stator_results/run3/results/comp/"), width = 16, height = 5)


setdiff(c(DE_stator$X[1:20], rev(c(DE_stator$X[(nrow(DE_stator)-19):nrow(DE_stator)]))), c(DE_standard$X[1:20], rev(c(DE_standard$X[(nrow(DE_standard)-19):nrow(DE_standard)]))))





