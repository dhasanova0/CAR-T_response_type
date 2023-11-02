#Author: Dzhansu Hasanova
#Merge Baseline Samples into Seurat object

#Import libraries
library(Seurat)
library(stringr)
library(SeuratDisk)
library(dplyr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Get Data
data_dir <- paste0(wd, 'data/GSE197268_RAW/baseline')
baseline_samples <- list.dirs(data_dir)
baseline_samples <- baseline_samples[-c(1)]

#Create lists to save names of samples and create Seurat objects
names <- c()
object_list <- c()

#Create a seurat object for each sample
for (i in baseline_samples){
  data <- Read10X(data.dir = i)
  object_list <- c(object_list, CreateSeuratObject(counts = data, project = str_sub(i, 65)))
  names <- c(names, str_sub(i, 65))
}

# Merge seurat objects of baseline samples
baseline_seurat <- merge(object_list[[1]], y = object_list[-c(1)], add.cell.ids = names, project = "baseline")
rm(object_list)
gc()

#Rename cells remove the name added when merging only batcodes left
baseline_seurat <- RenameCells(baseline_seurat, new.names = gsub(".*_","",colnames(baseline_seurat)))

#Add metadata
labels <- read_excel(paste0(wd, "data/metadata/labels.xlsx"))

baseline_seurat@meta.data$r_names <- rownames(baseline_seurat@meta.data)
baseline_seurat@meta.data <- merge(baseline_seurat@meta.data, labels, by ="orig.ident", all.x = TRUE)
rownames(baseline_seurat@meta.data) <- baseline_seurat@meta.data$r_names
baseline_seurat@meta.data <- subset(baseline_seurat@meta.data, select=-c(r_names))

#Save Seurat Object
output <- paste0(wd, "data/output/")
if (file.exists(output)) {cat("The folder already exists")} else {dir.create(output)}
saveRDS(baseline_seurat, paste0(wd, "data/output/baseline_raw_md.rds"))
