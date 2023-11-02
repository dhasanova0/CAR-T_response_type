#Author: Dzhansu Hasanova 10162023
#Integration and TCR/BCR removal and btach effect correction

# Filter BCR/TCR genes

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seurat_obj <- readRDS(paste0(wd, "/data/output/baseline_doublet_filtered_md_annot.rds"))
seurat_obj_no_rbc <- readRDS(paste0(wd, "data/output/baseline_doublet_no_rbc_filtered_md_annot.rds"))

filter_tcr_bcr <- function(seurat_object){
  
  #BCR filtering
  seurat_object <- seurat_object[!grepl("^IG[HKL]V", rownames(seurat_object)), ]
  seurat_object <- seurat_object[!grepl("^IG[HKL]J", rownames(seurat_object)), ]
  seurat_object <- seurat_object[!grepl("^IG[HKL]C", rownames(seurat_object)), ]
  seurat_object <- seurat_object[!grepl("^IGH[ADEGM]", rownames(seurat_object)), ]
  
  #TCR filtering
  seurat_object <- seurat_object[!grepl("^TR[ABDG][VJC]", rownames(seurat_object)), ]
  
  return(seurat_object)
}

seurat_obj <- filter_tcr_bcr(seurat_obj)
seurat_obj_no_rbc <- filter_tcr_bcr(seurat_obj_no_rbc)

#Normalize data

normalisation <- function(seurat_object){
  

  seurat_object <- SCTransform(seurat_object)
  seurat_object <- RunPCA(object = seurat_object)
  seurat_object <- FindNeighbors(object = seurat_object, dims = 1:20)
  seurat_object <- FindClusters(object = seurat_object)
  seurat_object <- RunUMAP(object = seurat_object, dims = 1:20)
  
  return(seurat_object)
  
} 


#Perform normalization
seurat_obj <- normalisation(seurat_obj)
seurat_obj_no_rbc <- normalisation(seurat_obj_no_rbc)

output1 <- paste0(wd, "data/output/baseline_doublet_filtered_md_annot_norm.rds")
output2 <- paste0(wd, "data/output/baseline_doublet_no_rbc_filtered_md_annot_norm.rds")
if (file.exists(output1)) {cat("The folder already exists")} else {dir.create(output1)}

saveRDS(seurat_obj, output1)
saveRDS(seurat_obj_no_rbc, output2)


integration <- function(seurat_object){
  # Perform integration on sample level
  
  
  ifnb.list <- SplitObject(seurat_object, split.by = "orig.ident")
  ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
  
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                           anchor.features = features)
  immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
  
  immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
  immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:20)
  
  immune.combined.sct <- FindNeighbors(object = immune.combined.sct, dims = 1:20)
  immune.combined.sct <- FindClusters(object = immune.combined.sct)
  
  return(immune.combined.sct)
  
}

#Perform integration
seurat_obj_batch_corr <- integration(seurat_obj)
seurat_obj_no_rbc_batch_corr <- integration(seurat_obj_no_rbc)

output1 <- paste0(wd, "data/output/baseline_doublet_filtered_md_annot_norm_batch.rds")
output2 <- paste0(wd, "data/output/baseline_doublet_no_rbc_filtered_md_annot_norm_batch.rds")

saveRDS(seurat_obj, output1)
saveRDS(seurat_obj_no_rbc, output2)

