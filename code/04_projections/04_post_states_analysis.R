#Author: Dzhansu Hasanova
#Projection of states from baseline subset 1 onto post-treatmetn samples
#Analysis of states in remaining cells

set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbreak) 
library(patchwork)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/post_state_baseline.rds"))
states <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))


#Get cell composition of post-treatment samples and create boxplot

seu_obj@meta.data$orig.ident_label_product <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$Response_classification, "_", seu_obj@meta.data$Product)
seu_obj@meta.data$cell_label <- paste0(seu_obj@meta.data$cell_type, "_", seu_obj@meta.data$Response_classification)
seu_obj@meta.data$orig.ident_label <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$Response_classification)

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
  ylab("Counts") + scale_y_break(c(10000, 25000))
ggsave("boxplot_cell_post.png", path = paste0(wd, "figures/posttreatment/"), width = 10, height = 6)



#Specify in what state the cells are
seu_obj <- get_state(seu_obj, states)

length(rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:2"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:3"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:4"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:1"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:5"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:6"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:7"] == "yes",]))


#Create UMAPs
cols1 = c("#ff0000", "#1F78C8", "#6A33C2","#33a02c", "#36648B", "#FFD700", "#ff7f00", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#00E2E5", "#C814FA", "#C8308C",  "#0000FF", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#00FF00")



DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:1", cols = c("gray", "#1F78C8")) + ggtitle("State 1") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state1.png", path =paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:2", cols = c("gray", "#1F78C8")) + ggtitle("State 2") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state2.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:3", cols = c("gray", "#1F78C8")) + ggtitle("State 3") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state3.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:4", cols = c("gray", "#1F78C8")) + ggtitle("State 4") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state4.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:5", cols = c("gray", "#1F78C8")) + ggtitle("State 5") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state5.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:6", cols = c("gray", "#1F78C8")) + ggtitle("State 6") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state6.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:7", cols = c("gray", "#1F78C8")) + ggtitle("State 7") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state7.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "cell_type", cols = cols1) + ggtitle("Cell Types") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_cell.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Response_classification", cols = cols1) + ggtitle("Response Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_label.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Product", cols = c("#ff7f00", "#1F78C8")) + ggtitle("Product") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_product.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "orig.ident") + ggtitle("Sample ID") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_sampleID.png", path = paste0(wd, "figures/posttreatment/states/umap/"), width = 7, height = 5)

#Perform DE for states which are different between baseline and post-treatment

Idents(seu_obj) <- "Cluster:1"
DE_1 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_1 <- DE_1[DE_1$p_val_adj<0.05,]
write.csv(DE_1, paste0(wd, "figures/posttreatment/states/DE/DE_s1.csv"))

Idents(seu_obj) <- "Cluster:3"
DE_3 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_3 <- DE_3[DE_3$p_val_adj<0.05,]
write.csv(DE_3, paste0(wd, "figures/posttreatment/states/DE/DE_s3.csv"))

Idents(seu_obj) <- "Cluster:4"
DE_4 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_4 <- DE_4[DE_4$p_val_adj<0.05,]
write.csv(DE_4, paste0(wd, "figures/posttreatment/states/DE/DE_s4.csv"))

Idents(seu_obj) <- "Cluster:5"
DE_5 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_5 <- DE_5[DE_5$p_val_adj<0.05,]
write.csv(DE_5, paste0(wd, "figures/posttreatment/states/DE/DE_s5.csv"))

Idents(seu_obj) <- "Cluster:2"
DE_2 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_2 <- DE_2[DE_2$p_val_adj<0.05,]
write.csv(DE_2, paste0(wd, "figures/posttreatment/states/DE/DE_s2.csv"))

Idents(seu_obj) <- "Response_classification"
DE_label <- FindMarkers(seu_obj, ident.1 = "R")
DE_label <- DE_label[DE_label$p_val_adj<0.05,]
DE_label$gene<-rownames(DE_label)
write.csv(DE_label, paste0(wd, "figures/posttreatment/states/DE/DE_RvsNR.csv"))




#Subset T cells
Idents(seu_obj) <- "cell_type"
T_seu <- subset(x = seu_obj, idents = c("CD4 T", "CD8 T", "Infusion T"))

T_seu <- NormalizeData(T_seu, normalization.method = "LogNormalize", scale.factor = 10000)
T_seu <- FindVariableFeatures(T_seu, selection.method = "vst", nfeatures = 2000)
T_seu <- ScaleData(T_seu)


Idents(T_seu) <- "Response_classification"
DE_T_label <- FindMarkers(T_seu, ident.1 = "R")
DE_T_label <- DE_T_label[DE_T_label$p_val_adj<0.05,]
DE_T_label$gene<-rownames(DE_T_label)
write.csv(DE_T_label, paste0(wd, "figures/posttreatment/states/DE/DE_RvsNR_T.csv"))

Idents(T_seu) <- "Cluster:2"
DE_T_s2<- FindMarkers(T_seu, ident.1 = "yes")
DE_T_s2 <- DE_T_s2[DE_T_s2$p_val_adj<0.05,]
DE_T_s2$gene<-rownames(DE_T_s2)
write.csv(DE_T_s2, paste0(wd, "figures/posttreatment/states/DE/DE_s2_T.csv"))



Idents(T_seu) <- "orig.ident" 
DE_NR_69 <- FindMarkers(object = T_seu, ident.1 = c("Patient14-D7", "Patient20-D7", "Patient25-D7", "Patient31-D7" ), logfc.threshold = 0.25)
DE_NR_69<-DE_NR_69[DE_NR_69$p_val_adj<0.05,]
DE_NR_69$gene<-rownames(DE_NR_69)

Idents(T_seu) <- "orig.ident" 
DE_R_69 <- FindMarkers(object = T_seu, ident.1 = c("Patient12-D7", "Patient21-D7", "Patient30-D7"), logfc.threshold = 0.25)
DE_R_69 <- DE_R_69[DE_R_69$p_val_adj<0.05,]
DE_R_69$gene<-rownames(DE_R_69)

Idents(T_seu) <- "orig.ident" 
DE_69 <- FindMarkers(object = T_seu, ident.1 = c("Patient14-D7", "Patient20-D7", "Patient25-D7", "Patient31-D7" ,"Patient12-D7", "Patient21-D7", "Patient30-D7"), logfc.threshold = 0.25)
DE_69 <- DE_69[DE_69$p_val_adj<0.05,]
DE_69$gene<-rownames(DE_69)

Idents(T_seu) <- "orig.ident" 
DE_69_R_NR <- FindMarkers(object = T_seu, ident.1 = c("Patient14-D7", "Patient20-D7", "Patient25-D7", "Patient31-D7" ), ident.2 = c("Patient12-D7", "Patient21-D7", "Patient30-D7"), logfc.threshold = 0.25)
DE_69_R_NR <- DE_69_R_NR[DE_69_R_NR$p_val_adj<0.05,]
DE_69_R_NR$gene<-rownames(DE_69_R_NR)


DoHeatmap(T_seu, features = c("FOXP3", 'PFCD1', "SLAMF6", "LAG3", "TGFB1")) + NoLegend()

Idents(seu_obj) <- "Response_classification"
VlnPlot(seu_obj, features = "CD69", pt.size = 0, cols = c("dodgerblue3", "darkorange2")) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means()
ggsave("violin_CD69_label.png", path = paste0(wd, "figures/posttreatment/"), width = 7, height = 5)

seu_obj@meta.data$label_ID <- paste0(seu_obj@meta.data$Response_classification,"_", seu_obj@meta.data$orig.ident)

Idents(seu_obj) <- "label_ID"
levels(seu_obj) <- unique(seu_obj@meta.data$label_ID)[order(unique(seu_obj@meta.data$label_ID))]

VlnPlot(seu_obj, features = c("IL10"), pt.size = 0) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means()

VlnPlot(seu_obj, features = c("TGFB1", "LAG3", "PDCD1", "CTLA4", "FOXP3", "IL10", "SLAMF6", "CD244", "CD160", "HAVCR2"), pt.size = 0) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means()
ggsave("violin_CD69_sample.png", path = paste0(wd, "figures/posttreatment/"), width = 7, height = 5)

Idents(seu_obj) <- "Response_classification"
DoHeatmap(seu_obj, features = c("TGFB1", "LAG3", "PDCD1", "CTLA4", "FOXP3", "IL10"))
