#Author: Dzhansu Hasanova 10042023
#Subset data into 10k R, 10k NR

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Define functions for downsampling to 20k (10k R, 10k NR)
#Make sure that number of cells are distributed as in full dataset

subset_10k <- function(seurat_object){
  
  # Calculate the target number of cells for each patient (assuming a list of patient IDs in the 'patient_id' column)
  target_cells_per_patient <- 10000 * table(seurat_object@meta.data$orig.ident) / ncol(seurat_object)
  target_cells_per_patient <- round(target_cells_per_patient)
  
  ident.list <- SplitObject(seurat_object, split.by = "orig.ident")
  
  for (i in 1:length(ident.list)){
    
    ident.list[[i]] <- ident.list[[i]][, sample(colnames(ident.list[[i]]), size =target_cells_per_patient[[i]], replace=F)]
    
  }
  
  subset_1 <- merge(ident.list[[1]], y = ident.list[-c(1)], add.cell.ids = ls(ident.list), project = "baseline")
  
  return(subset_1)
  
}

subsample_seurat <- function(seurat_object, output_path){
  
  seurat_NR <- subset(seurat_object, subset = label == "NR")
  seurat_R <- subset(seurat_object, subset = label == "R")
  
  n_splits <- min(length(colnames(seurat_NR)), length(colnames(seurat_R))) / 10000

  
  for (i in 1:round(n_splits)){
    gc()
    
    if (length(colnames(seurat_NR)) & length(colnames(seurat_R)) > 10000){
      
      
      subset1_NR <- subset_10k(seurat_NR)
      subset1_R <- subset_10k(seurat_R)
      
      subset_1 <- merge(subset1_NR, y = subset1_R, add.cell.ids = c("NR", "R"),  project = paste0("subset_",i))
      subset_1 <- RenameCells(subset_1, new.names = gsub(".*_","",colnames(subset_1)))
      
      saveRDS(subset_1, paste0(output_path, "subset_",i,".rds"))
      rm(subset_1)
      gc()
      subset1_NR <- RenameCells(subset1_NR, new.names = gsub(".*_","",colnames(subset1_NR)))
      subset1_R <- RenameCells(subset1_R, new.names = gsub(".*_","",colnames(subset1_R)))
      seurat_R <- RenameCells(seurat_R, new.names = gsub(".*_","",colnames(seurat_R)))
      seurat_NR <- RenameCells(seurat_NR, new.names = gsub(".*_","",colnames(seurat_NR)))
      
      seurat_NR <- seurat_NR[,!colnames(seurat_NR) %in% colnames(subset1_NR)]
      seurat_R <- seurat_R[,!colnames(seurat_R) %in% colnames(subset1_R)]
      
      
      gc()
      
    }
  }
  
  subset1_NR <- RenameCells(subset1_NR, new.names = gsub(".*_","",colnames(subset1_NR)))
  subset1_R <- RenameCells(subset1_R, new.names = gsub(".*_","",colnames(subset1_R)))
  seurat_R <- RenameCells(seurat_R, new.names = gsub(".*_","",colnames(seurat_R)))
  seurat_NR <- RenameCells(seurat_NR, new.names = gsub(".*_","",colnames(seurat_NR)))
  
  seurat_NR <- seurat_NR[,!colnames(seurat_NR) %in% colnames(subset1_NR)]
  seurat_R <- seurat_R[,!colnames(seurat_R) %in% colnames(subset1_R)]
  
  subset_4 <- merge(seurat_NR, y = seurat_R, add.cell.ids = c("NR", "R"),  project = paste0("subset_",i+1))
  saveRDS(subset_4, paste0(output_path, "susbet", i+1 , "_not_balanced.rds"))
  
}

#Define output paths and create if doesn't exist
output_path1 <- paste0(wd, "data/output/stator_input/rds_new/")
output_path2 <- paste0(wd, "data/output/stator_input_no_rbc/rds_new/")

if (file.exists(output_path1)) {cat("The folder already exists")} else {dir.create(output_path1)}
if (file.exists(output_path2)) {cat("The folder already exists")} else {dir.create(output_path2)}

#Load data for downsampling and save downsampled files
seurat_obj <- readRDS(paste0(wd, "data/output/baseline_doublet_filtered_md_annot.rds"))
subsample_seurat(seurat_obj, output_path1)
rm(seurat_obj)
gc()

#Load data for downsampling and save downsampled files
seurat_obj_no_rbc <- readRDS(paste0(wd, "data/output/baseline_doublet_no_rbc_filtered_annot_md.rds"))
subsample_seurat(seurat_obj_no_rbc, output_path2)
rm(seurat_obj_no_rbc)
gc()






