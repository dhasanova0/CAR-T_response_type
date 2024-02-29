#Author: Dzhansu Hasanova
#Initial analysis of baseline samples

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library("ggpubr")
library(enrichR)
library(EnhancedVolcano)

set.seed(12342)
library("ggstatsplot")


cols1 = c("#1F78C8", "#33a02c", "#ff0000", "#6A33C2","#ff7f00", "#36648B", "#FFD700", "#a6cee3",
          "#FB6496", "#999999", "#b2df8a", "#CAB2D6", "#FDBF6F", "#EEE685", "#565656",
          "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#00E2E5", "#00FF00")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Load Seurat object
seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state.rds"))

#Change metadata to get cell distribution plot
seu_obj@meta.data$orig.ident_label_product <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$label, "_", seu_obj@meta.data$product)
seu_obj@meta.data$cell_label <- paste0(seu_obj@meta.data$cell_type, "_", seu_obj@meta.data$label)
seu_obj@meta.data$orig.ident_label <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$label)

#Create Cell distribution plot
cell_id <- data.frame(table(seu_obj@meta.data[, c("cell_label", "orig.ident_label")]))
cell_id$orig.ident_label <- as.character(cell_id$orig.ident_label)

cell_id$label_from_cell <- sub(".*_", "", cell_id$cell_label)
cell_id$label_from_id <- sub(".*_", "", cell_id$orig.ident_label)
cell_id <- cell_id[cell_id$label_from_id==cell_id$label_from_cell,]

cell <- data.frame(table(seu_obj@meta.data[, c("cell_type", "orig.ident_label")]))
cell$orig.ident_label <- as.character(cell$orig.ident_label)
cell$label <- sub('.*_', '', cell$orig.ident_label)


ggplot(cell, aes(fill=label, y=Freq, x=cell_type)) + 
  geom_boxplot(position="dodge", alpha=0.5, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.8), 
             aes(color=factor(label)), size = 1.5, alpha = 1) +
  scale_fill_manual(name = "Response Type", values=c("deepskyblue4", "#E69F00")) +
  scale_color_manual(name = "Response Type", values=c("deepskyblue4", "#E69F00")) +
  theme_bw()  +
  theme(axis.text=element_text(size=12), legend.text=element_text(size=12))+
  xlab("") +
  ylab("Counts")
ggsave("boxplot_cell.png", path = paste0(wd, "figures/baseline/analysis/", width = 10, height = 5))


#Create Dimplots

DimPlot(seu_obj, group.by = "cell_type", reduction = "umap", cols = cols1) +ggtitle("Cell Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_celltype.png", path = paste0(wd,"figures/baseline/analysis/"), width = 6, height = 5)

DimPlot(seu_obj, group.by = "label", reduction = "umap", cols = cols1) +ggtitle("Response Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_label.png", path = paste0(wd,"figures/baseline/analysis/"), width = 6, height = 5)

DimPlot(seu_obj, group.by = "seurat_clusters", reduction = "umap", label = "TRUE", cols = cols1) + ggtitle("Seurat Cluster") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("cluster.png", path = paste0(wd,"figures/baseline/analysis/"), width = 6, height = 5)

DimPlot(seu_obj, group.by = "orig.ident", reduction = "umap", cols = cols1) + ggtitle("Patient ID") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_patient.png", path = paste0(wd,"figures/baseline/analysis/"), width = 6, height = 5)


#Perform DE for R and NR
Idents(seu_obj) <- "label"
DE <- FindMarkers(seu_obj, ident.1 = "R")
DE <- DE[DE$p_val_adj<0.05,]

write.csv(DE, paste0(wd, "figures/baseline/analysis/DE_RvsNR_log.csv"))

DE <- DE[order(DE$avg_log2FC),]

EnhancedVolcano(DE , 
                rownames(DE ),
                x ="avg_log2FC",
                y ="p_val_adj",
                selectLab = c(rownames(DE)[1:5], rev(rownames(DE)[(nrow(DE)-3):nrow(DE)])),
                title = "R vs NR",
                subtitle = 'Baseline samples',
                FCcutoff = 0.5,
                labSize = 4.0)
ggsave("DE_volcano.png", path = paste0(wd, "figures/baseline/analysis/", width = 8, height = 7))



Idents(seu_obj_filtered) <- "label"
DE_filtered <- FindMarkers(seu_obj_filtered, ident.1 = "R")
DE_filtered <- DE_filtered[DE_filtered$p_val_adj<0.05,]


EnhancedVolcano(DE_filtered , 
                rownames(DE_filtered ),
                x ="avg_log2FC", 
                y ="p_val_adj",
                FCcutoff = 0.5,
                labSize = 4.0)
ggsave("DE_volcano_filtered.png", path = paste0(wd, "figures/baseline/analysis/"), width = 6, height = 7)

#Perform DE for between R and NR in T cells

T_seu <- subset(seu_obj, subset = cell_type == c("CD4 T", "CD8 T", "T cell"))

T_seu <- NormalizeData(T_seu, normalization.method = "LogNormalize", scale.factor = 10000)
T_seu <- FindVariableFeatures(T_seu, selection.method = "vst", nfeatures = 2000)
T_seu <- ScaleData(T_seu)

Idents(T_seu) <- "label"
DE_1_T <- FindMarkers(object = T_seu, ident.1 = "R", logfc.threshold = 0.25)
DE_1_T <- DE_1_T[DE_1_T$p_val_adj<0.05,]
DE_1_T$genes <- rownames(DE_1_T)
write.csv(DE_1_T, paste0(wd, "figures/baseline/analysis/DE_T_R_NR.csv"))

#Perform DE for between R and NR in Monocytes

mono <- subset(seu_obj, subset = cell_type == c("Monocyte"))

mono <- NormalizeData(mono, normalization.method = "LogNormalize", scale.factor = 10000)
mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = 2000)
mono <- ScaleData(mono)

Idents(mono) <- "label"
DE_mono <- FindMarkers(object = mono, ident.1 = "R", logfc.threshold = 0.25)
DE_mono <- DE_mono[DE_mono$p_val_adj<0.05,]
DE_mono$genes <- rownames(DE_mono)
write.csv(DE_mono, paste0(wd, "figures/baseline/analysis/DE_mono_R_NR.csv"))

#Perform DE for between R samples who have and didn't recieve ASCT
seu_R <- subset(seu_obj, subset = label == "R")

Idents(seu_R) <- "autologous_transplant"
DE_asct <- FindMarkers(seu_R, ident.1 = "yes")
DE_asct <- DE_asct[DE_asct$p_val_adj<0.05,]

write.csv(DE_asct, paste0(wd, "figures/baseline/analysis/DE_asct_log.csv"))

DE_asct <- DE_asct[order(DE_asct$avg_log2FC),]

EnhancedVolcano(DE_asct , 
                rownames(DE_asct ),
                x ="avg_log2FC", 
                y ="p_val_adj",
                selectLab = c(rownames(DE)[1:5], rev(rownames(DE_asct)[(nrow(DE_asct)-6):nrow(DE_asct)])),
                title = "ASCT vs no ASCT",
                subtitle = "Baseline samples from R ",
                FCcutoff = 0.5,
                labSize = 4.0)
ggsave("DE_volcano_asct.png", path = paste0(wd, "figures/baseline/analysis/"), width = 8, height = 7)


#Check clusters which seem to have more cells from a given response type in UMAP

#Check cluster 5

Idents(seu_obj) <- "seurat_clusters"

C5 <- subset(seu_obj, idents = 5)
table(C5@meta.data$label)
table(C5@meta.data$cell_type)
table(C5@meta.data$orig.ident_label_product)

#Check cluster 13

C13 <- subset(seu_obj, idents = 13)
table(C13@meta.data$label)
table(C13@meta.data$cell_type)
table(C13@meta.data$orig.ident_label_product)

#Check cluster 15

C15 <- subset(seu_obj, idents = 15)
table(C15@meta.data$label)
table(C15@meta.data$cell_type)
table(C15@meta.data$orig.ident_label_product)#Check cluster 18



