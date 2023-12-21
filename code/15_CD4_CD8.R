library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
set.seed(123)
ls <- read.csv("/Users/dhasanova/Documents/ETH/HS23/figures/stator_run1/CellList-2023-11-23.csv")
ls_0.5 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/run1/CellList-2023-12-08_0.5.csv")
md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_subtype.csv")
md_corr <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_corr.csv")
seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/stator_input/rds_new/subset_1.rds")

DE_findmarkers <- function(seu_obj, cluster){
  DE <- FindMarkers(seu_obj, ident.1 = cluster)
  DE <-DE[DE$p_val_adj<0.05,]
  DE$gene<-rownames(DE)
  
  DE$cluster[DE$avg_log2FC>0]<-paste0("up_",cluster)
  DE$cluster[DE$avg_log2FC<0]<-paste0("up_","others")
  
  return(DE)
}

c <- ls[ls$State == paste0('cluster_C:',"1"), ]
selected_cluster_list_optimal <- strsplit(c$CellID, ",\\s*")[[1]]

md_c <- md_corr[md_corr$X%in%selected_cluster_list_optimal,]
seu_obj_1 <- subset(seu_obj, subset = barcodes %in% md_c$X)
md_c <- md_c[match(seu_obj_1@meta.data$barcodes, md_c$X),]
seu_obj_1 <- AddMetaData(seu_obj_1, md_c$Cell.Types , col.name = "cell_corrected")


DefaultAssay(seu_obj_1) <- "RNA"
seu_obj_1 <- NormalizeData(seu_obj_1, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj_1 <- FindVariableFeatures(seu_obj_1, selection.method = "vst", nfeatures = 2000)
seu_obj_1 <- ScaleData(seu_obj_1)
seu_obj_1 <- RunPCA(seu_obj_1, features = VariableFeatures(object = seu_obj_1))
seu_obj_1 <- FindNeighbors(object = seu_obj_1, dims = 1:10)
seu_obj_1 <- FindClusters(object = seu_obj_1, resolution = 0.2)
seu_obj_1 <- RunUMAP(object = seu_obj_1, dims = 1:20)

#seu_obj_1 <- SCTransform(seu_obj_1)
#seu_obj_1 <- RunPCA(object = seu_obj_1)

DimPlot(seu_obj_1, reduction = "pca", dims = c(1, 2), group.by = "cell_corrected")

DimPlot(seu_obj_1, reduction = "umap", dims = c(1, 2), group.by = "cell_type")


DE_c0 <- DE_findmarkers(seu_obj_1, 0)
DE_c1 <- DE_findmarkers(seu_obj_1, 1)
DE_c2 <- DE_findmarkers(seu_obj_1, 2)
DE_c3 <- DE_findmarkers(seu_obj_1, 3)

# --------------- Dice distance 0.5 ------------------

c1 <- ls_0.5[ls_0.5$State == paste0('cluster_C:',"1"), ]
selected_cluster_list1 <- strsplit(c1$CellID, ",\\s*")[[1]]

c2 <- ls_0.5[ls_0.5$State == paste0('cluster_C:',"2"), ]
selected_cluster_list2 <- strsplit(c2$CellID, ",\\s*")[[1]]

c3 <- ls_0.5[ls_0.5$State == paste0('cluster_C:',"3"), ]
selected_cluster_list3 <- strsplit(c3$CellID, ",\\s*")[[1]]

selected_cluster_list <- unique(c(selected_cluster_list1, selected_cluster_list2, selected_cluster_list3))

md_c <- md_corr[md_corr$X %in% selected_cluster_list ,]

seu_obj_1 <- subset(seu_obj, subset = barcodes %in% md_c$X)
md_c <- md_c[match(seu_obj_1@meta.data$barcodes, md_c$X),]
seu_obj_1 <- AddMetaData(seu_obj_1, md_c$Cell.Types , col.name = "cell_corrected")


table(md_c$Cell.Types)

table(md_c$X %in% md_c_optimal$X)



