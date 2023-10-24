#Author: Dzhansu Hasanova 10232023
#Cell annotation, manual and SingleR

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(DoubletFinder)
library(celldex)
library(SingleR)

cell_annotation <- read_excel("Documents/ETH/HS23/data/cell_annotation_R.xlsx")

norm <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_normalised_TCR_BCR_filtered.rds")
norm_no_rbc <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_no_rbc_normalised_TCR_BCR_filtered.rds")

batchcorr <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_normalised_TCR_BCR_filtered_batchcorr.rds")
batchcorr_no_rbc <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_subsampled_no_rbc_normalised_TCR_BCR_filtered_batchcorr.rds")

# Manual annotation
DimPlot(norm, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)

DimPlot(no_correction, reduction = 'umap', group.by = "seurat_clusters", label = TRUE) +
  DimPlot(no_correction, reduction = 'umap', group.by = "orig.ident")

cluster1.markers <- FindMarkers(norm, ident.1 = 1, min.pct = 0.25) #Monocytes
cluster7.markers <- FindMarkers(norm, ident.1 = 7, min.pct = 0.25)
cluster9.markers <- FindMarkers(norm, ident.1 = 9, min.pct = 0.25) #tbd
cluster3.markers <- FindMarkers(norm, ident.1 = 3, min.pct = 0.25) #Monocytes
cluster6.markers <- FindMarkers(norm, ident.1 = 6, min.pct = 0.25) #Monoctes
cluster10.markers <- FindMarkers(norm, ident.1 = 10, min.pct = 0.25) #Dendritic cells

cluster14.markers <- FindMarkers(norm, ident.1 = 14, min.pct = 0.25) #?

cluster13.markers <- FindMarkers(norm, ident.1 = 13, min.pct = 0.25) #? 

cluster12.markers <- FindMarkers(norm, ident.1 = 12, min.pct = 0.25) #B cells

cluster5.markers <- FindMarkers(norm, ident.1 = 5, min.pct = 0.25) #T cells
cluster2.markers <- FindMarkers(norm, ident.1 = 2, min.pct = 0.25)
cluster8.markers <- FindMarkers(norm, ident.1 = 8, min.pct = 0.25) #Cytotoxic T cells
cluster0.markers <- FindMarkers(norm, ident.1 = 0, min.pct = 0.25) #Cytotoxic T cells
cluster4.markers <- FindMarkers(norm, ident.1 = 4, min.pct = 0.25) #NK cells
cluster11.markers <- FindMarkers(norm, ident.1 = 11, min.pct = 0.25) #NK cells

table(norm@meta.data$seurat_clusters)

no_correction@meta.data$manual_annotation <- no_correction@meta.data$seurat_clusters

no_correction@meta.data$manual_annotation[no_correction@meta.data$manual_annotation == "0"] <- "CD8+"

no_correction@meta.data$manual_annotation <-replace(no_correction@meta.data$manual_annotation, 
                                                    no_correction@meta.data$manual_annotation==0, "CD8+")

names(new.cluster.ids) <- levels(no_correction)
no_correction <- RenameIdents(no_correction, new.cluster.ids)
DimPlot(no_correction, reduction = "umap", label = TRUE, pt.size = 0.5)

if (!dir.exists("/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/")) dir.create("/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/")
saveRDS(norm_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/norm_manann.rds")
saveRDS(norm_no_rbc_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/norm_norbc_manann.rds")

saveRDS(batchcorr_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/batchcorr_manann.rds")
saveRDS(batchcorr_no_rbc_annot, "/Users/dhasanova/Documents/ETH/HS23/data/output/annotated/batchcorr_norbc_manann.rds")







