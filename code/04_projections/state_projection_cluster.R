#Author: Dzhansu Hasanova
#Projection of interactions (d-tuples) on remaining cells. If a given expression pattern
#of an interaction is in a cell then the cell is labeled with the respective interaction
#This script was run on the university cluster eddie (HPC)
#Files in data/output_baseline/states were generated with this script with according state and Seurat objects as input

set.seed(12342)

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

wd = "/Users/dhasanova/Documents/ETH/HS23/"

#states <- read.csv("/exports/igmm/eddie/ponting-lab/dzhansu/data/output/counts/State_Table-2023-11-17.csv")
#seu_obj <- readRDS("/exports/igmm/eddie/ponting-lab/dzhansu/data/output/counts/baseline_doublet_filtered_md_annot.rds")

states <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))
seu_obj <- readRDS(paste0(wd, "data/output_baseline/baseline_doublet_filtered_md_annot.rds"))

seu_obj <- RenameCells(seu_obj, new.names = gsub(".*_","",colnames(seu_obj)))

#Re-normalise
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj), assay = 'RNA')
seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
seu_obj <- RunUMAP(seu_obj, dims = 1:20)

#Get states
states$Genes <- strsplit(states$Genes, "_")
states$D.tuple <- gsub("\\D", "", states$D.tuple)


raw_counts <-as.data.frame(seu_obj@assays[["RNA"]]@counts)

# Get genes in interaction
for(l in 1:nrow(states)){
  print(l)
  list_1 <- c()
  list_0 <- c()
  MFI <- c()
  for (i in 1:length(states$Genes[[l]])){
    print(i)
    if (substr(states$D.tuple[l], i, i) == '1'){
      list_1 <- c(list_1,states$Genes[[l]][i])
    } else if (substr(states$D.tuple[l], i, i) == '0'){
      list_0 <- c(list_0,states$Genes[[l]][i])
    }
  }
  print(list_1)
  print(list_0)	
  if (length(list_1) > 0){
    MFI <- raw_counts[, apply(raw_counts[list_1, ], 2, function(x) all(x > 0))]
  }else {
   print("else")
   MFI <- raw_counts
  }
  if (length(list_0) > 0){
    MFI <- MFI[, apply(MFI[list_0, ], 2, function(x) all(x == 0))]
  }

  MFI <- colnames(MFI)
  
  seu_obj@meta.data[ , ncol(seu_obj@meta.data) + 1] <- ifelse(rownames(seu_obj@meta.data) %in% MFI, "yes", "no")
  colnames(seu_obj@meta.data)[ncol(seu_obj@meta.data)] <- paste0('MFI_', l)
}


#saveRDS(seu_obj, "/exports/igmm/eddie/ponting-lab/dzhansu/data/output/counts/baseline_state.rds")

saveRDS(seu_obj, paste0(wd, "data/output_baseline/states/baseline_state.rds"))

