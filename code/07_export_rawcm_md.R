#Sanity checks for Stator inputs

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)

set.seed(100)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

input_path1 <- paste0(wd, "data/output/stator_input/rds_new/")
input_path2 <- paste0(wd, "data/output/stator_input_no_rbc/rds_new/")

output_path1 <- (paste0(wd, "data/output/stator_input/raw_md/"))
output_path2 <- paste0(wd, "data/output/stator_input_no_rbc/raw_md/")

if (file.exists(output_path1)) {cat("The folder already exists")} else {dir.create(output_path1)}
if (file.exists(output_path1)) {cat("The folder already exists")} else {dir.create(output_path2)}


save_csv <- function(input_path, output_path){
  
  ls <- list.files(input_path, full.names=TRUE)
  seurat_object <- lapply(ls, readRDS)
  
  for (i in 1:length(seurat_object)){
    
    write.csv(t(as.data.frame(seurat_object[[i]]@assays[["RNA"]]@counts)), file=paste0(output_path, "subset",i,"_raw.csv"), quote=FALSE)
    write.csv(seurat_object[[i]]@meta.data, file=paste0(output_path, "subset",i,"_md.csv"), quote=FALSE)
    gc()
  }
  rm(seurat_object)
  gc()
}

save_csv(input_path1, output_path1)
save_csv(input_path2, output_path2)




