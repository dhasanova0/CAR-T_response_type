#Author: Dzhansu Hasanova 10232023
#Annotation on QC data (no normalization and integration)
#Cell annotation with Haradhvala et.al. 

library(Seurat)
library("stringr")
library(SeuratDisk)
library(dplyr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output/baseline_raw_md.rds"))

#Import annotations from Haradhvala et al.
annotations_haradhvala <- read.csv(paste0(wd, "data/metadata/Fig1_CART_global_obs.csv"), row.names=1)
#Subset baseline samples
annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_baseline$barcodes <- rownames(annotations_baseline)
rownames(annotations_baseline) <- NULL
annotations_baseline <- annotations_baseline[ , c("cell_type", "barcodes")]


orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcodes <- rownames(seurat_object@meta.data)

annotations_baseline <- annotations_baseline %>%
  filter(annotations_baseline$barcodes %in% orig_rownames)

md <- merge(annotations_baseline, seurat_object@meta.data, by ="barcodes", all.y = TRUE)
md <- md[match(orig_rownames, md$barcodes),]


seurat_object <- AddMetaData(seurat_object, md$cell_type , col.name = "cell_type")

saveRDS(seurat_object, paste0(wd, "data/output/baseline_raw_md_annot.rds"))

