#Author: Dzhansu Hasanova 10232023
#Annotation on QC data (no normalization and integration)
#Cell annotation with Haradhvala et.al. 

library(Seurat)
library("stringr")
library(SeuratDisk)
library(dplyr)
library(readr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output/baseline_raw_md.rds"))

annotations_haradhvala_sub <- read_delim("Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
annotations_baseline_sub <- annotations_haradhvala_sub[grep("Baseline", annotations_haradhvala_sub$timepoint), ]
annotations_baseline_sub <- annotations_baseline_sub[c("...1", "subtype", "class")]
annotations_baseline_sub$class <- paste0(annotations_baseline_sub$subtype, "_", annotations_baseline_sub$class)

#Import annotations from Haradhvala et al.
annotations_haradhvala <- read.csv(paste0(wd, "data/metadata/Fig1_CART_global_obs.csv"), row.names=1)
#Subset baseline samples
annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_baseline$barcodes <- rownames(annotations_baseline)
rownames(annotations_baseline) <- NULL

annotations_baseline <- left_join(annotations_baseline, annotations_baseline_sub, by = join_by(barcodes == ...1))
annotations_baseline$subtype <- ifelse(is.na(annotations_baseline$subtype), annotations_baseline$cell_type, annotations_baseline$subtype)
annotations_baseline$class <- ifelse(is.na(annotations_baseline$class), annotations_baseline$cell_type, annotations_baseline$class)

annotations_baseline <- annotations_baseline[ , c("subtype", "class","barcodes")]


orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcodes <- rownames(seurat_object@meta.data)

annotations_baseline <- annotations_baseline %>%
  filter(annotations_baseline$barcodes %in% orig_rownames)

md <- merge(annotations_baseline, seurat_object@meta.data, by ="barcodes", all.y = TRUE)
md <- md[match(orig_rownames, md$barcodes),]


seurat_object <- AddMetaData(seurat_object, md$subtype , col.name = "cell_type")
seurat_object <- AddMetaData(seurat_object, md$class , col.name = "subtype")


saveRDS(seurat_object, paste0(wd, "data/output/baseline_raw_md_annot_corr.rds"))

