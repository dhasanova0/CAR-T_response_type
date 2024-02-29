#Run on eddie (TO DO: Get script from eddie)

set.seed(12342)

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(irlba)
library(Matrix)

s1_md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_corr.csv")
states <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/run1/State_Table-2023-11-17.csv")

genes <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/counts/genes.csv")


seu_obj_corr <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/baseline_doublet_filtered_md_corrannot_norm.rds")
seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/baseline_doublet_filtered_md_annot_umap.rds")

seu_obj <- RenameCells(seu_obj, new.names = gsub(".*_","",colnames(seu_obj)))
seu_obj@meta.data$cell_type_corr <- seu_obj_corr@meta.data$cell_type_corr



DefaultAssay(seu_obj) <- "RNA"

#Log normalize
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj), assay = 'RNA')
seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
seu_obj <- RunUMAP(seu_obj, dims = 1:20)

states$Genes <- strsplit(states$Genes, "_")
states$D.tuple <- gsub("\\D", "", states$D.tuple)

raw_counts <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/counts/baseline_counts.csv")

# Get genes in interaction
for(l in 1:nrow(states)){
  list_1 <- c()
  list_0 <- c()
  MFI <- c()
  for (i in 1:length(states$Genes[[l]])){
    
    if (substr(states$D.tuple[l], i, i) == '1'){
      list_1 <- c(list_1,states$Genes[[l]][i])
    } else if (substr(states$D.tuple[l], i, i) == '0'){
      list_0 <- c(list_0,states$Genes[[l]][i])
    }
  }
  if (length(list_1) > 0){
    MFI <- raw_counts[, apply(raw_counts[list_1, ], 2, function(x) all(x > 0))]
  }
  if (length(list_0) > 0){
    MFI <- colnames(MFI[, apply(MFI[list_0, ], 2, function(x) all(x == 0))])
  }
  
  seu_obj@meta.data[ , ncol(seu_obj@meta.data) + 1] <- ifelse(rownames(seu_obj@meta.data) %in% MFI, "yes", "no")
  colnames(seu_obj@meta.data)[ncol(seu_obj@meta.data)] <- paste0('MFI_', l)
}




