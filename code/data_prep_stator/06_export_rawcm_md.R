#Author: Dzhansu Hasanova 10042023
#Save csv files of count matrices for Stator input

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)

set.seed(100)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Define input paths
input_path1 <- paste0(wd, "data/output_baseline/stator_input/rds_new/")

input_T <- paste0(wd, "data/output_baseline/stator_input_T/rds")

input_path_post_axi_D7 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/rds_D7/")
input_path_post_tisa_D7 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/rds_D7/")

input_path_post_axi_D14 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/rds_D14/")
input_path_post_tisa_D14 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/rds_D14/")

#Define output paths
output_path1 <- (paste0(wd, "data/output_baseline/stator_input/all_cells/"))

output_T <- paste0(wd, "data/output_baseline/stator_input_T/")

output_path_post_axi_D7 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/raw_md_D7/")
output_path_post_tisa_D7 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/raw_md_D7/")

output_path_post_axi_D14 <- paste0(wd, "data/output_posttreatment/stator_input/axi_cel/raw_md_D14/")
output_path_post_tisa_D14 <- paste0(wd, "data/output_posttreatment/stator_input/tisa_cel/raw_md_D14/")



save_csv <- function(input_path, output_path){
  #Input: the subsampled rds files
  #Save csv files of transposed count matrices for stator (cell x gene)
  #Save csv files of metadata
  
  ls <- list.files(input_path, full.names=TRUE)
  seurat_object <- lapply(ls, readRDS)
  
  for (i in 1:length(seurat_object)){
    
    write.csv(t(as.data.frame(seurat_object[[i]]@assays[["RNA"]]@counts)), file=paste0(output_path, "subset",i,"_raw.csv"), quote=FALSE)
    write.csv(as.data.frame(seurat_object[[i]]@meta.data), file=paste0(output_path, "subset",i,"_md.csv"), quote=TRUE)
    gc()
  }
  rm(seurat_object)
  gc()
}

save_csv(input_path1, output_path1)
save_csv(input_path2, output_path2)

save_csv(input_T, output_T)

save_csv(input_path_post_axi_D7, output_path_post_axi_D7)
save_csv(input_path_post_tisa_D7, output_path_post_tisa_D7)

save_csv(input_path_post_axi_D14, output_path_post_axi_D14)
save_csv(input_path_post_tisa_D14, output_path_post_tisa_D14)




