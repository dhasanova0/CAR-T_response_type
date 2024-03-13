#Author: Dzhansu Hasanova 
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
#Make sure that number of cells are distributed as in original dataset

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
  #Subsample 10k R and 10k NR
  
  #Get R and NR cells
  seurat_NR <- subset(seurat_object, subset = label == "NR")
  seurat_R <- subset(seurat_object, subset = label == "R")
  
  n_splits <- min(length(colnames(seurat_NR)), length(colnames(seurat_R))) / 10000

  #Split data in batches to approx. 20k cells total preserving cell type distribution in patients
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
  
  #Rename cells - only barcode remains
  subset1_NR <- RenameCells(subset1_NR, new.names = gsub(".*_","",colnames(subset1_NR)))
  subset1_R <- RenameCells(subset1_R, new.names = gsub(".*_","",colnames(subset1_R)))
  seurat_R <- RenameCells(seurat_R, new.names = gsub(".*_","",colnames(seurat_R)))
  seurat_NR <- RenameCells(seurat_NR, new.names = gsub(".*_","",colnames(seurat_NR)))
  
  #Create the last subset with the remianing cells
  seurat_NR <- seurat_NR[,!colnames(seurat_NR) %in% colnames(subset1_NR)]
  seurat_R <- seurat_R[,!colnames(seurat_R) %in% colnames(subset1_R)]
  
  subset_4 <- merge(seurat_NR, y = seurat_R, add.cell.ids = c("NR", "R"),  project = paste0("subset_",i+1))
  saveRDS(subset_4, paste0(output_path, "susbet", i+1 , "_not_balanced.rds"))
  
}




#Load data for downsampling and save downsampled files baseline
#Define output paths and create if doesn't exist
output_path1_baseline <- paste0(wd, "data/output_baseline/stator_input/rds_new/")


if (file.exists(output_path1_baseline)) {cat("The folder already exists")} else {dir.create(output_path1_baseline)}


seurat_obj <- readRDS(paste0(wd, "data/output_baseline/baseline_doublet_filtered_md_annot.rds"))
subsample_seurat(seurat_obj, output_path1_baseline)
rm(seurat_obj)
gc()

#Load data for downsampling and save downsampled files post-treatment

output_path_posttreatment <- paste0(wd, "data/output_posttreatment/stator_input/combined/rds/")

output_path_post_axi_D7 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/rds_D7/")
output_path_post_tisa_D7 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/rds_D7/")
output_path_post_axi_D14 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/rds_D14/")
output_path_post_tisa_D14 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/rds_D14/")

if (file.exists(output_path_posttreatment)) {cat("The folder already exists")} else {dir.create(output_path_posttreatment)}
if (file.exists(output_path_post_axi)) {cat("The folder already exists")} else {dir.create(output_path_post_axi)}
if (file.exists(output_path_post_tisa)) {cat("The folder already exists")} else {dir.create(output_path_post_tisa)}

seurat_obj_post <- readRDS(paste0(wd, "data/output_posttreatment/post_doublet_filtered_md_annot.rds"))
colnames(seurat_obj_post@meta.data)[26] <- "label"

#Get the 2 axi and tisa-cel samples
D14_axi_cel <- subset(seurat_obj_post, subset = orig.ident == "Patient12-D14" | orig.ident == "Patient14-D14")
D14_tisa_cel <- subset(seurat_obj_post, subset = orig.ident == "Patient20-D14" | orig.ident == "Patient21-D14")

D7_axi_cel <- subset(seurat_obj_post, subset = Product == "Axi-cel")
D7_axi_cel <- subset(D7_axi_cel, subset = orig.ident == "Patient12-D14" | orig.ident == "Patient14-D14", invert = TRUE) #Remova D14 samples
D7_tisa_cel <- subset(seurat_obj_post, subset = Product == "Tisa-cel")
D7_tisa_cel <- subset(D7_tisa_cel, subset = orig.ident == "Patient20-D14" | orig.ident == "Patient21-D14", invert = TRUE) #Remova D14 samples

subsample_seurat(seurat_obj_post, output_path_posttreatment)
gc()


subsample_seurat(D7_axi_cel, output_path_post_axi_D7)
subsample_seurat(D7_tisa_cel, output_path_post_tisa_D7)

#No need to subsample D14 samples, since they already contain around 20k cells
saveRDS(D14_tisa_cel,paste0(output_path_post_tisa_D14, "D14_tisa.rds"))
saveRDS(D14_axi_cel, paste0(output_path_post_axi_D14, "D14_axi.rds"))

#T cells baseline
if (file.exists(output_path_baseline_T)) {cat("The folder already exists")} else {dir.create(output_path_baseline_T)}
output_path_baseline_T <- paste0(wd, "data/output_baseline/stator_input_T/rds/")
seu_obj_base_T <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/seu_T.rds")
subsample_seurat(seu_obj_base_T, output_path_baseline_T)
rm(seu_obj_base_T)
gc()




