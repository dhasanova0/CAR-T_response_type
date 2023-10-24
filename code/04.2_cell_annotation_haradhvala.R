#Author: Dzhansu Hasanova 10232023
#Cell annotation with Haradhvala et.al. 


library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(DoubletFinder)
library(celldex)
library(SingleR)

norm <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_normalised_TCR_BCR_filtered.rds")
norm_no_rbc <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_no_rbc_normalised_TCR_BCR_filtered.rds")

batchcorr <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_normalised_TCR_BCR_filtered_batchcorr.rds")
batchcorr_no_rbc <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_no_rbc_normalised_TCR_BCR_filtered_batchcorr.rds")


#Import annotations from Haradhvala et al.
annotations_haradhvala <- read.csv("~/Documents/ETH/HS23/data/Fig1_CART_global_obs.csv", row.names=1)

#Subset baseline samples
annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_baseline$barcodes <- rownames(annotations_baseline)
rownames(annotations_baseline) <- NULL

annotate <- function(seurat_object, ref_annotation){
  
  #Match metadata file to annotation file
  seurat_object@meta.data$barcodes <- gsub(".*_","",rownames(seurat_object@meta.data))
  
  #Keep rownames after merge
  seurat_object@meta.data$r_names <- rownames(seurat_object@meta.data)

  #Combine metadata and cell annotations
  seurat_object@meta.data <- merge(seurat_object@meta.data, ref_annotation, by ="barcodes", all.x = TRUE)

  rownames(seurat_object@meta.data) <- seurat_object@meta.data$r_names
  seurat_object@meta.data <- subset(seurat_object@meta.data, select=-c(r_names))
  
  return(seurat_object)
  
}

norm_annot <- annotate(norm, annotations_baseline)
norm_no_rbc_annot <- annotate(norm_no_rbc, annotations_baseline)

batchcorr_annot <- annotate(batchcorr, annotations_baseline)
batchcorr_no_rbc_annot <- annotate(batchcorr_no_rbc, annotations_baseline)

if (!dir.exists("/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/")) dir.create("/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/")
saveRDS(norm_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/norm_refann.rds")
saveRDS(norm_no_rbc_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/norm_norbc_refann.rds")

saveRDS(batchcorr_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/batchcorr_refann.rds")
saveRDS(batchcorr_no_rbc_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/batchcorr_norbc_refann.rds")


DimPlot(batchcorr_annot, reduction = "umap", group.by = "cell_type", label = TRUE)
DimPlot(no_correction, reduction = "umap", group.by = "cell_type", label = TRUE)

batch_correction@meta.data <- batch_correction@meta.data %>% mutate(cell_type = ifelse(is.na(cell_type), "NA", cell_type))
na_batch_correction <- subset(batch_correction, subset = cell_type == "NA")

FeatureScatter(na_batch_correction, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "seurat_clusters")

DimPlot(na_batch_correction, reduction = "umap", group.by = "seurat_clusters")

