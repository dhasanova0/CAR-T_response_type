set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(EnhancedVolcano)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


seu_obj <- readRDS(paste0(wd, "data/output_posttreatment/states/post_state_axiD7.rds"))
states <- read.csv(paste0(wd, "data/stator_results/axi_D7_1/Shiny_results/State_Table-2024-01-03.csv"))
md <- read.csv(paste0(wd, "data/stator_results/axi_D7_1/md_cell_label.csv"))

seu_obj <- get_state(seu_obj, states)

seu_obj@meta.data$timepoint <- sapply(strsplit(seu_obj@meta.data$orig.ident, "-"), function(x) x[2])

Idents(seu_obj) <- "timepoint"

seu_obj_axi_D7 <- subset(seu_obj, ident = "D7", subset = Product == "Axi-cel")

#These steps are performed on eddie (output: seu_obj_axi_D7_umap.rds)
#DefaultAssay(seu_obj_axi_D7) <- "RNA"
#seu_obj_n[['SCT']] <- NULL

#seu_obj_axi_D7 <- NormalizeData(seu_obj_axi_D7, normalization.method = "LogNormalize", scale.factor = 10000)
#seu_obj_axi_D7 <- FindVariableFeatures(seu_obj_axi_D7, selection.method = "vst", nfeatures = 2000)
#seu_obj_axi_D7 <- ScaleData(seu_obj_axi_D7)

#seu_obj_axi_D7 <- RunPCA(seu_obj_axi_D7, features = VariableFeatures(object = seu_obj), assay = 'RNA')
#seu_obj_axi_D7 <- FindNeighbors(seu_obj_axi_D7, dims = 1:20)
#seu_obj_axi_D7 <- FindClusters(seu_obj_axi_D7, resolution = 0.5)
#seu_obj_axi_D7 <- RunUMAP(seu_obj_axi_D7, dims = 1:20, assay = DefaultAssay(object = seu_obj))

seu_obj_axi_D7 <- readRDS(paste0(wd, "data/output_posttreatment/states/seu_obj_axi_D7_umap.rds"))

#Number of cells in at least one state
length(rownames(seu_obj_axi_D7@meta.data[seu_obj_axi_D7@meta.data[, "Cluster:2"] == "yes" |
                                           seu_obj_axi_D7@meta.data[, "Cluster:3"] == "yes" |
                                           seu_obj_axi_D7@meta.data[, "Cluster:4"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:1"] == "yes" |
                                           seu_obj_axi_D7@meta.data[, "Cluster:5"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:6"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:7"] == "yes" |
                                           seu_obj_axi_D7@meta.data[, "Cluster:8"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:9"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:10"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:11"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:12"] == "yes"|
                                           seu_obj_axi_D7@meta.data[, "Cluster:13"] == "yes",]))


#Remove cells from which states were inferred and re-normalize

Idents(seu_obj_axi_D7) <- "barcode"
seu_obj_axi_D7 <- subset(seu_obj_axi_D7, ident = md$X, invert = TRUE)

seu_obj_axi_D7 <- NormalizeData(seu_obj_axi_D7, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj_axi_D7 <- FindVariableFeatures(seu_obj_axi_D7, selection.method = "vst", nfeatures = 2000)
seu_obj_axi_D7 <- ScaleData(seu_obj_axi_D7)

saveRDS(seu_obj_axi_D7, paste0(wd, "data/output_posttreatment/states/seu_obj_axi_D7_valid.rds"))


cols1 = c("#ff0000", "#1F78C8", "#6A33C2","#33a02c", "#36648B", "#FFD700", "#ff7f00", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#00E2E5", "#C814FA", "#C8308C",  "#0000FF", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#00FF00")


DimPlot(seu_obj_axi_D7, reduction = "umap", group.by = "Cluster:1", cols = c("gray", "#1F78C8")) + ggtitle("State 1") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state1.png", path =paste0(wd, "data/output_posttreatment/states/axi_D7/umap/"), width = 6, height = 5)

DimPlot(seu_obj_axi_D7, reduction = "umap", group.by = "Cluster:10", cols = c("gray", "#1F78C8")) + ggtitle("State 10") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state10.png", path =paste0(wd, "data/output_posttreatment/states/axi_D7/umap/"), width = 6, height = 5)

DimPlot(seu_obj_axi_D7, reduction = "umap", group.by = "Cluster:9", cols = c("gray", "#1F78C8")) + ggtitle("State 9") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state9.png", path =paste0(wd, "data/output_posttreatment/states/axi_D7/umap/"), width = 6, height = 5)

DimPlot(seu_obj_axi_D7, reduction = "umap", group.by = "cell_type", cols = cols1) + ggtitle("Cell Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_cell_type.png", path =paste0(wd, "data/output_posttreatment/states/axi_D7/umap/"), width = 6, height = 5)

Idents(seu_obj_axi_D7) <- "Cluster:9"
DE_c9<- FindMarkers(seu_obj_axi_D7, ident.1 = "yes")
DE_c9 <- DE_c9[DE_c9$p_val_adj<0.05,]
DE_c9$gene<-rownames(DE_c9)
write.csv(DE_c9, paste0(wd, "data/output_posttreatment/states/axi_D7/DE_c9.csv"))

EnhancedVolcano(DE_c9 , 
                rownames(DE_c9),
                x ="avg_log2FC",
                y ="p_val_adj",
                title = "Cells in state 9 vs all other cells",
                subtitle = '',
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_C9_vs_all.png", path = paste0(wd, "data/output_posttreatment/states/axi_D7/"), width = 10, height = 10)

Idents(seu_obj_axi_D7) <- "Cluster:10"
DE_c10<- FindMarkers(seu_obj_axi_D7, ident.1 = "yes")
DE_c10 <- DE_c10[DE_c10$p_val_adj<0.05,]
DE_c10$gene<-rownames(DE_c10)
write.csv(DE_c10, paste0(wd, "data/output_posttreatment/states/axi_D7/DE_c10.csv"))

table(seu_obj_axi_D7@meta.data$`Cluster:10`, seu_obj_axi_D7@meta.data$cell_type)

EnhancedVolcano(DE_c10 , 
                rownames(DE_c10),
                x ="avg_log2FC",
                y ="p_val_adj",
                title = "Cells in state 10 vs all other cells",
                subtitle = '',
                FCcutoff = 0.5,
                labSize = 6.0)

ggsave("DE_C10_vs_all.png", path = paste0(wd, "data/output_posttreatment/states/axi_D7/"), width = 10, height = 10)

axi_D7@meta.data$label_ID <- paste0(axi_D7@meta.data$Response_classification,"_", axi_D7@meta.data$orig.ident)

Idents(axi_D7) <- "Response_classification"
levels(axi_D7) <- unique(axi_D7@meta.data$label_ID)[order(unique(axi_D7@meta.data$label_ID))]

VlnPlot(axi_D7, features = c("CD69"), pt.size = 0) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) 

