#Author: Dzhansu Hasanova
#Merge Baseline Samples into Seurat object

#Import libraries
library(Seurat)
library(stringr)
library(SeuratDisk)
library(dplyr)
library("readxl")

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"
samples_baseline <- 'data/GSE197268_RAW/baseline'
samples_post <- 'data/GSE197268_RAW/posttreatment'


merge_data <- function(path, project){
  #Get Data
  data_dir <- path
  samples <- list.dirs(data_dir)
  samples <- samples[-c(1)]
  
  #Create lists to save names of samples and create Seurat objects
  names <- c()
  object_list <- c()
  
  #Create a seurat object for each sample
  for (i in samples){
    data <- Read10X(data.dir = i)
    object_list <- c(object_list, CreateSeuratObject(counts = data, project = sub(".*/", "", i)))
    names <- c(names, sub(".*/", "", i))
  }
  
  # Merge seurat objects of baseline samples
  seurat_obj <- merge(object_list[[1]], y = object_list[-c(1)], add.cell.ids = names, project = project)
  rm(object_list)
  gc()
  
  #Rename cells remove the name added when merging only barcodes left
  seurat_obj <- RenameCells(seurat_obj, new.names = gsub(".*_","",colnames(seurat_obj)))
  return(seurat_obj)
}

baseline_seurat <- merge_data(paste0(wd, samples_baseline), "baseline")
post_seurat <- merge_data(paste0(wd, samples_post), "post")

#Add metadata
labels <- read_excel(paste0(wd, "data/metadata/labels.xlsx"))

md_post <- read_excel(paste0(wd, "data/metadata/metadata.xlsx"))
colnames(md_post) <- gsub(" ", "_", colnames(md_post))

md_post$patient_no <- sub(".*-", "", md_post$Sample_ID)
post_seurat@meta.data <- post_seurat@meta.data[,1:3]
post_seurat@meta.data$patient_no <- gsub(".*?(\\d+)-(D7|D14)", "\\1", post_seurat@meta.data$orig.ident)
post_seurat@meta.data$patient_no <- ifelse(nchar(post_seurat@meta.data$patient_no) == 1, 
                                           paste0("0", post_seurat@meta.data$patient_no), post_seurat@meta.data$patient_no)


add_md <- function(md, seurat_obj, merge_by){
  
  orig_rownames <- rownames(seurat_obj@meta.data)
  seurat_obj@meta.data$r_names <- rownames(seurat_obj@meta.data)
  seurat_obj@meta.data <- merge(seurat_obj@meta.data, md, by = merge_by, all.x = TRUE)
  seurat_obj@meta.data <- seurat_obj@meta.data[match(orig_rownames, seurat_obj@meta.data$r_names),]
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$r_names
  seurat_obj@meta.data <- subset(seurat_obj@meta.data, select=-c(r_names))
  
  return(seurat_obj)
  
}


baseline_seurat <- add_md(labels, baseline_seurat, "orig.ident")
post_seurat <- add_md(md_post, post_seurat, "patient_no")


#Save Seurat Object
output_baseline <- paste0(wd, "data/output_baseline/")
if (file.exists(output_baseline)) {cat("The folder already exists")} else {dir.create(output_baseline)}
saveRDS(baseline_seurat, paste0(wd, "data/output_baseline/baseline_raw_md.rds"))

#Save Seurat Object
output_post <- paste0(wd, "data/output_posttreatment/")
if (file.exists(output_post)) {cat("The folder already exists")} else {dir.create(output_post)}
saveRDS(post_seurat, paste0(wd, "data/output_posttreatment/post_raw_md.rds"))

