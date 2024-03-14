#Author: Dzhansu Hasanova
#Perform pseudobulk DE analysis for R and NR baseline samples

set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

wd = "/Users/dhasanova/Documents/ETH/HS23/"

seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state.rds"))


#Function for pseudobulk DE with DESeq2
DE_pseudo <- function(seu_obj, group){
  # pseudobulk the counts based on donor-condition-celltype
  pseudo_obj <- AggregateExpression(seu_obj, assays = "RNA", return.seurat = T, 
                                    group.by = group, slot="counts")
  
  
  df <- as.matrix(pseudo_obj@assays[["RNA"]]@counts)
  coldata <- pseudo_obj@meta.data
  
  dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = coldata,
                                design = ~ orig.ident)
  
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  res <- na.omit(res)
  res<-res[res$padj<0.05,]
  
  return(res)
}

DE_RvsNR <- DE_pseudo(seu_obj, c("label", "orig.ident"))
write.csv(DE_RvsNR, paste0(wd, "figures/baseline/analysis/pseudo_RvsNR.csv"))


EnhancedVolcano(DE_RvsNR , 
                rownames(DE_RvsNR),
                x ="log2FoldChange",
                y ="padj",
                title = "Pseudobulk R vs NR samples",
                subtitle = 'Baseline',
                drawConnectors = TRUE,
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_pseudo.png", path = paste0(wd, "figures/baseline/analysis/pseudo"), width = 12, height = 10)





