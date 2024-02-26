library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)

ls <- read.csv("/Users/dhasanova/Documents/ETH/HS23/figures/stator_run1/CellList-2023-11-23.csv")
md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md.csv")
seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_doublet_filtered_md_annot_norm.rds")

DefaultAssay(seu_obj) <- "RNA"

DimPlot(seu_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE)

Idents(object = seu_obj) <- seu_obj@meta.data$cell_type



CD4.markers <- FindMarkers(seu_obj, ident.1 = "CD4 T", ident.2 = "CD8 T", min.pct = 0.25)

CD4.markers_annot <- FindMarkers(seu_obj, ident.1 = "CD4 T", min.pct = 0.25)

CD4.markers$genes <- rownames(CD4.markers)
CD4.markers_annot$genes <- rownames(CD4.markers_annot)

setdiff(rownames(CD4.markers), rownames(CD4.markers_annot))

c1 <- ls[ls$State == 'cluster_C:1', ]
selected_cluster_list1<-strsplit(c1$CellID, ",\\s*")[[1]]



md_c1 <- md_CD8[md_CD8$X%in%selected_cluster_list1,]

md_c1_NR <- md_c1[md_c1$Cell.State == "NR",]
md_c1_R <- md_c1[md_c1$Cell.State == "R",]

DE<-FindMarkers(object = seu_obj,ident.1 = md_c1_NR$X,ident.2 = md_c1_R$X,logfc.threshold = 0.25)
DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_","NR")
DE$cluster[DE$avg_log2FC<0]<-paste0("up_","R")

print(table(DE$cluster))
