#Author: Dzhansu Hasanova
#Projection of states from baseline subset 2 onto the remaining cell
#Analysis of states in remaining cells
set.seed(12342)

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Get data with interactions
seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_run2.rds"))
states <- read.csv(paste0(wd, "data/stator_results/run3/results/State_Table-2024-02-01.csv"))
md <- read.csv(paste0(wd, "data/stator_results/run3/Shiny/md_cell_label.csv"))

seu_obj <- get_state(seu_obj, states)

#Subset cells not in 
Idents(seu_obj) <- "barcodes"
seu_obj <- subset(seu_obj, idents = md$X, invert = TRUE)

#Log normalize
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

saveRDS(seu_obj, paste0(wd, "data/output_baseline/states/baseline_state_run2_valid.rds"))

cells_state_NR <- rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:2"] == "yes" |
                                               seu_obj@meta.data[, "Cluster:8"] == "yes" |
                                               seu_obj@meta.data[, "Cluster:3"] == "yes" |
                                               seu_obj@meta.data[, "Cluster:1"] == "yes"|
                                               seu_obj@meta.data[, "Cluster:6"] == "yes",])



cells_state_45791011 <- rownames(seu_obj@meta.data[seu_obj@meta.data[, "Cluster:4"] == "yes" | 
                                                seu_obj@meta.data[, "Cluster:5"] == "yes" |
                                                seu_obj@meta.data[, "Cluster:7"] == "yes" |
                                                seu_obj@meta.data[, "Cluster:9"] == "yes" |
                                                seu_obj@meta.data[, "Cluster:10"] == "yes" |
                                                seu_obj@meta.data[, "Cluster:11"] == "yes",])


seu_obj@meta.data$NR_state <- NA
seu_obj@meta.data$NR_state<- ifelse(rownames(seu_obj@meta.data) %in% cells_state_NR, "NR_state", "no")

overlap<-intersect(cells_state_NR,cells_state_45791011)

cells_state_NR<-cells_state_NR[!cells_state_NR%in%overlap]


#Perform DE between states enriched in NR and all other states create heatmap of top 20 DE gene NR states and all other cells
Idents(seu_obj) <- "barcodes"

DE<-FindMarkers(object = seu_obj,ident.1 = cells_state_NR,logfc.threshold = 0.25)

DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_",paste0("NR",collapse = "_"))
DE$cluster[DE$avg_log2FC<0]<-paste0("up_",paste0("other",collapse = "_"))

write.csv(DE, paste0(wd, "data/stator_results/run3/results/valid/DE_NRvsother_entire.csv"))

EnhancedVolcano(DE , 
                rownames(DE ),
                x ="avg_log2FC",
                y ="p_val_adj",
                title = "Cells in NR state vs all other cells",
                subtitle = 'Baseline samples',
                FCcutoff = 0.5,
                labSize = 6.0)
ggsave("DE_volcano_entire.png", path = paste0(wd, "data/stator_results/run3/results/valid/"), width = 12, height = 10)

DE <- DE[order(DE$avg_log2FC),]

top20_NR <- DE[1:20,]
top20_R <- DE[(nrow(DE)-19):nrow(DE),]

Idents(seu_obj) <- "NR_state"
levels(seu_obj) <- c("NR_state", "no")

DoHeatmap(seu_obj, features = c(rownames(top20_NR), rownames(top20_R)), cells = c(cells_state_NR), group.colors = c("darkorange2", "dodgerblue3"), slot = "counts",
          disp.min = 0,disp.max = 4) +scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
ggsave("DE_R_NR_entire.png", path = paste0(wd, "data/stator_results/run3/results/valid/"), width = 7, height = 7)


#Create UMAPs of states
cols1 = c("#1F78C8", "#ff0000", "#33a02c", "#6A33C2","#36648B", "#ff7f00", "#FFD700", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#00E2E5", "#00FF00")

DimPlot(seu_obj, reduction = "umap", group.by = "NR_state") + ggtitle("Cells in R or NR states") + ylab("UMAP 2") + xlab("UMAP 1") + 
  scale_color_manual(labels=c('other', 'NR State'), values = c("#1F78C8", "#ff0000"))
ggsave("umap_NR_state.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:1", cols = c("gray", "#1F78C8")) + ggtitle("State 1") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state1.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:2", cols = c("gray", "#1F78C8")) + ggtitle("State 2") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state2.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:3", cols = c("gray", "#1F78C8")) + ggtitle("State 3") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state3.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:4", cols = c("gray", "#1F78C8")) + ggtitle("State 4") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state4.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:5", cols = c("gray", "#1F78C8")) + ggtitle("State 5") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state5.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:6", cols = c("gray", "#1F78C8")) + ggtitle("State 6") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state6.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:7", cols = c("gray", "#1F78C8")) + ggtitle("State 7") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state7.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:8", cols = c("gray", "#1F78C8")) + ggtitle("State 8") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state8.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:9", cols = c("gray", "#1F78C8")) + ggtitle("State 9") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state9.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:10", cols = c("gray", "#1F78C8")) + ggtitle("State 10") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state10.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)

DimPlot(seu_obj, reduction = "umap", group.by = "Cluster:11", cols = c("gray", "#1F78C8")) + ggtitle("State 11") + ylab("UMAP 2") + xlab("UMAP 1") + NoLegend()
ggsave("umap_state11.png", path = paste0(wd, "data/stator_results/run3/results/umap/"), width = 6, height = 5)


#Cluster 8 DE (activated cell state)
Idents(seu_obj) <- "cell_type"
T_seu <- subset(x = seu_obj, idents = c("CD4 T", "CD8 T", "T cell"))

Idents(seu_obj) <- "Cluster:8"

DE_c8 <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
DE_c8 <- DE_c8[DE_c8$p_val_adj<0.05,]
DE_c8$genes <- rownames(DE_c8)
write.csv(DE_c8, paste0(wd, "figures/stator_run1/thesis/valid/DE/DE_c8_T_vs_all.csv"))


