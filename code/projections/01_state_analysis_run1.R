#Author: Dzhansu Hasanova
#Projection of states from baseline subset 1 onto the remaining cell
#Analysis of states in remaining cells
set.seed(12342)

library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)
library(ggsignif)
library(ggpubr)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Get data with interactions
seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state.rds"))
states <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))
states_0.5 <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-0.5.csv"))
states_0.5 <- states_0.5[match(states$X, states_0.5$X),]
md <- read.csv(paste0(wd, "data/stator_results/run1/md/subset1_md_cell.csv"))

Idents(seu_obj) <- "barcodes"
seu_obj <- subset(seu_obj, idents = md$X, invert = TRUE )

#Log normalize
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

#saveRDS(seu_obj, paste0(wd, "data/output_baseline/counts/baseline_state_valid.rds"))
#seu_obj <- readRDS(paste0(wd, "data/output_baseline/counts/baseline_state_valid.rds"))

seu_obj <- get_state(seu_obj, states)

seu_obj_0.5 <- seu_obj
seu_obj_0.5@meta.data <- seu_obj_0.5@meta.data %>% select(-contains('Cluster'))
seu_obj_0.5 <- get_state(seu_obj_0.5, states_0.5)

#Number of cells in at least one state
length(rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:2"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:3"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:4"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:1"] == "yes" |
                                    seu_obj@meta.data[, "Cluster:5"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:6"] == "yes"|
                                    seu_obj@meta.data[, "Cluster:7"] == "yes",]))

cells_state_R <- rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:2"] == "yes" |
                                            seu_obj@meta.data[, "Cluster:3"] == "yes" |
                                            seu_obj@meta.data[, "Cluster:4"] == "yes" ,])


cells_state_NR <- rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:1"] == "yes" |
                                              seu_obj@meta.data[, "Cluster:5"] == "yes" |
                                              seu_obj@meta.data[, "Cluster:6"] == "yes" |
                                              seu_obj@meta.data[, "Cluster:7"] == "yes",])

cells_state_1 <- rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:1"] == "yes",])

#Remove cells with both R and NR state
overlap<-intersect(cells_state_R,cells_state_NR)

cells_state_R<-cells_state_R[!cells_state_R%in%overlap]
cells_state_NR<-cells_state_NR[!cells_state_NR%in%overlap]

seu_obj@meta.data$R_NR_state <- NA
seu_obj@meta.data$R_NR_state<- ifelse(rownames(seu_obj@meta.data) %in% cells_state_R, "R_state", NA)
seu_obj@meta.data[is.na(seu_obj@meta.data$R_NR_state),]$R_NR_state <- ifelse(rownames(seu_obj@meta.data[is.na(seu_obj@meta.data$R_NR_state),]) %in% cells_state_NR, "NR_state", NA)

#Perform DE between states enriched in R and NR and create heatmap of top 20 DE gene in R and NR states
Idents(seu_obj) <- "barcodes"

DE<-FindMarkers(object = seu_obj,ident.1 = cells_state_R,ident.2 = cells_state_NR,logfc.threshold = 0.25)

DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_",paste0("R",collapse = "_"))
DE$cluster[DE$avg_log2FC<0]<-paste0("up_",paste0("NR",collapse = "_"))

write.csv(DE, paste0(wd, "figures/stator_run1/thesis/valid/DE_RvsNR_entire_valid.csv"))

EnhancedVolcano(DE , 
                rownames(DE ),
                x ="avg_log2FC",
                y ="p_val_adj",
                title = "R state vs NR state",
                subtitle = 'Baseline samples',
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_volcano_entire.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 12, height = 10)

DE <- DE[order(DE$avg_log2FC),]

top20_NR <- DE[1:20,]
top20_R <- DE[(nrow(DE)-19):nrow(DE),]

Idents(seu_obj) <- "R_NR_state"
levels(seu_obj) <- c("R_state", "NR_state")

DoHeatmap(seu_obj, features = c(rownames(top20_NR), rownames(top20_R)), cells = c(cells_state_R, cells_state_NR), group.colors = c("darkorange2", "dodgerblue3"), slot = "counts",
          disp.min = 0,disp.max = 4) +scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
ggsave("DE_R_NR_entire.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 7, height = 7)


#Create UMAPs of states
cols1 = c("#ff0000", "#1F78C8", "#6A33C2","#33a02c", "#36648B", "#FFD700", "#ff7f00", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#00E2E5", "#C814FA", "#C8308C",  "#0000FF", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#00FF00")

DimPlot(seu_obj, reduction = "umap", group.by = "R_NR_state") + ggtitle("Cells in R or NR states") + ylab("UMAP 2") + xlab("UMAP 1") + 
  scale_color_manual(labels=c('NR State', 'R State'), values = c("#1F78C8", "#ff0000"))
ggsave("umap_R_NR_state.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:1", cols = c("gray", "#1F78C8")) + ggtitle("State 1") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state1.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:2", cols = c("gray", "#1F78C8")) + ggtitle("State 2") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state2.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:3", cols = c("gray", "#1F78C8")) + ggtitle("State 3") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state3.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:4", cols = c("gray", "#1F78C8")) + ggtitle("State 4") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state4.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:5", cols = c("gray", "#1F78C8")) + ggtitle("State 5") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state5.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:6", cols = c("gray", "#1F78C8")) + ggtitle("State 6") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state6.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:7", cols = c("gray", "#1F78C8")) + ggtitle("State 7") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state7.png", path = paste0(wd, "figures/stator_run1/thesis/valid/states/"), width = 6, height = 5)

#Perform DE on clsuter 1
Idents(seu_obj) <- 'barcodes'
Idents(seu_obj) <- "Cluster:1"

DE<-FindMarkers(object = seu_obj,ident.1 = "yes",ident.2 = "no",logfc.threshold = 0.25)

DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_",paste0("C1",collapse = "_"))
DE$cluster[DE$avg_log2FC<0]<-paste0("up_",paste0("other",collapse = "_"))

write.csv(DE, paste0(wd, "figures/stator_run1/thesis/valid/DE_C1_RvsNR.csv"))

EnhancedVolcano(DE , 
                rownames(DE ),
                x ="avg_log2FC",
                y ="p_val_adj",
                drawConnectors = TRUE,
                selectLab = c("CD69", "JUN", "LTB", "IL7R", "TNFAIP3", "FOS", "CXCR4"),
                title = "Cells in state 1 vs all other cells",
                subtitle = 'States inferred from first subset',
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_C1_vs_all.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 15, height = 10)


#Check CD69 expression ViolinPlot between R and NR
Idents(seu_obj) <- "label"
VlnPlot(seu_obj, features = "CD69", pt.size = 0, cols = c("dodgerblue3", "darkorange2")) +
  stat_compare_means()
ggsave("CD69.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 4, height = 5)


#Check CD69 expression ViolinPlot in patient samples
seu_obj@meta.data$label_orig.ident <- paste0(seu_obj@meta.data$label, "_" ,seu_obj@meta.data$orig.ident)
Idents(seu_obj) <- "label_orig.ident"
levels(seu_obj) <- unique(seu_obj@meta.data$label_orig.ident)[order(unique(seu_obj@meta.data$label_orig.ident))]

VlnPlot(seu_obj, features = "CD69", pt.size = 0, cols = cols1) + NoLegend() + theme(axis.text.x = element_text(angle = 90))
ggsave("violin_CD69_sample.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 7, height = 5)

#Check expression of Cksuter 1 genes in R and NR
MFI_1 <- c("TSPYL2", "TNFAIP3", "JUN", "CD69", "CD247", "IL7R", "KLRG1", "CST3" )
DotPlot(seu_obj, features = MFI_1)
ggsave("DotPlot_C1.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 4, height = 5)

#---------------------- DE analyses ---------------------------

#DE analyses for states

Idents(seu_obj) <- "Cluster:1"
DE_1 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_1 <- DE_1[DE_1$p_val_adj<0.05,]
write.csv(DE_1, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s1.csv"))

Idents(seu_obj) <- "Cluster:3"
DE_3 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_3 <- DE_3[DE_3$p_val_adj<0.05,]
write.csv(DE_3, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s3.csv"))

Idents(seu_obj) <- "Cluster:4"
DE_4 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_4 <- DE_4[DE_4$p_val_adj<0.05,]
write.csv(DE_4, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s4.csv"))

Idents(seu_obj) <- "Cluster:5"
DE_5 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_5 <- DE_5[DE_5$p_val_adj<0.05,]
write.csv(DE_5, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s5.csv"))

Idents(seu_obj) <- "Cluster:2"
DE_2 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_2 <- DE_2[DE_2$p_val_adj<0.05,]
write.csv(DE_2, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s2.csv"))

#DE analyses T cells

Idents(seu_obj) <- "cell_type"
T_seu <- subset(x = seu_obj, idents = c("CD4 T", "CD8 T", "T cell"))

T_seu <- NormalizeData(T_seu, normalization.method = "LogNormalize", scale.factor = 10000)
T_seu <- FindVariableFeatures(T_seu, selection.method = "vst", nfeatures = 2000)
T_seu <- ScaleData(T_seu)

Idents(T_seu) <- "Cluster:1"
DE_1_T <- FindMarkers(object = T_seu, ident.1 = "yes", logfc.threshold = 0.25)
DE_1_T <- DE_1_T[DE_1_T$p_val_adj<0.05,]
DE_1_T$genes <- rownames(DE_1_T)
write.csv(DE_1_T, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s1_T.csv"))

EnhancedVolcano(DE_1_T , 
                rownames(DE_1_T ),
                x ="avg_log2FC",
                y ="p_val_adj",
                drawConnectors = TRUE,
                selectLab = c("FOS","CD69", "JUN", "LTB", "IL7R", "DUSP1", "NKG7", "GZMH",
                              "FGFBP2", "GNLY", "GZMB"),
                title = "T Cells in state 1 vs all other cells",
                subtitle = '',
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_C1_vs_all_T.png", path = paste0(wd, "figures/stator_run1/thesis/valid/"), width = 15, height = 10)


Idents(T_seu) <- "Cluster:2"
DE_2_T <- FindMarkers(object = T_seu, ident.1 = "yes", logfc.threshold = 0.25)
DE_2_T <- DE_2_T[DE_2_T$p_val_adj<0.05,]
write.csv(DE_2_T, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s2_T.csv"))

Idents(T_seu) <- "Cluster:3"
DE_3_T <- FindMarkers(object = T_seu, ident.1 = "yes", logfc.threshold = 0.25)
DE_3_T <- DE_3_T[DE_3_T$p_val_adj<0.05,]
write.csv(DE_3_T, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s3_T.csv"))

Idents(T_seu) <- "Cluster:4"
DE_4_T <- FindMarkers(object = T_seu, ident.1 = "yes", logfc.threshold = 0.25)
DE_4_T <- DE_4_T[DE_4_T$p_val_adj<0.05,]
write.csv(DE_4_T, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_s4_T.csv"))

Idents(T_seu) <- "label"
DE_label_T <- FindMarkers(object = T_seu, ident.1 = "R", logfc.threshold = 0.25)
DE_label_T <- DE_label_T[DE_label_T$p_val_adj<0.05,]
DE_label_T$genes <- rownames(DE_label_T)
write.csv(DE_label_T, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_RvsNR_T.csv"))

#DE for Monocytes

Idents(seu_obj) <- "cell_type"
mono <- subset(x = seu_obj, idents = c("Monocyte"))

mono <- NormalizeData(mono, normalization.method = "LogNormalize", scale.factor = 10000)
mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = 2000)
mono <- ScaleData(mono)

Idents(mono) <- "Cluster:6"
DE_mono_6 <- FindMarkers(object = mono, ident.1 = "yes", logfc.threshold = 0.25)
DE_mono_6 <- DE_mono_6[DE_mono_6$p_val_adj<0.05,]
DE_mono_6$genes <- rownames(DE_mono_6)
write.csv(DE_mono_6, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_mono_6.csv"))


Idents(mono) <- "Cluster:7"
DE_mono_7 <- FindMarkers(object = mono, ident.1 = "yes", logfc.threshold = 0.25)
DE_mono_7 <- DE_mono_7[DE_mono_7$p_val_adj<0.05,]
DE_mono_7$genes <- rownames(DE_mono_7)
write.csv(DE_mono_7, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_mono_7.csv"))



#---------------------- Dice Distance 0.5 Monocyte state 17 and CD4T cell state 5 ---------------------------


mono <- subset(x = seu_obj_0.5, subset = cell_type == "Monocyte")

mono <- NormalizeData(mono, normalization.method = "LogNormalize", scale.factor = 10000)
mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = 2000)
mono <- ScaleData(mono)

Idents(mono) <- "Cluster:17"
DE_mono_17_all <- FindMarkers(object = mono, ident.1 = "yes", logfc.threshold = 0.25)
DE_mono_17_all <- DE_mono_17_all[DE_mono_17_all$p_val_adj<0.05,]
DE_mono_17_all$genes <- rownames(DE_mono_17_all)
write.csv(DE_mono_17_all, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_mono_17_vs_all.csv"))


CD4 <- subset(x = seu_obj_0.5, subset = cell_type == "CD4 T")

CD4 <- NormalizeData(CD4, normalization.method = "LogNormalize", scale.factor = 10000)
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", nfeatures = 2000)
CD4 <- ScaleData(CD4)

Idents(CD4) <- "Cluster:5"
DE_CD4_c5_all <- FindMarkers(object = CD4, ident.1 = "yes", logfc.threshold = 0.25)
DE_CD4_c5_all <- DE_CD4_c5_all[DE_CD4_c5_all$p_val_adj<0.05,]
DE_CD4_c5_all$genes <- rownames(DE_CD4_c5_all)
write.csv(DE_CD4_c5_all, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_CD4_c5_vs_all.csv"))



