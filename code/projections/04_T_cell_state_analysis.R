#Author: Dzhansu Hasanova
#Projection of states from baseline T cell subset onto the remaining T cells cell
#Analysis of states in remaining cells

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

seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/base_T_cell.rds"))
states <- read.csv(paste0(wd,"data/stator_results/T_cell/shiny/State_Table-2024-02-22.csv"))
md <- read.csv(paste0(wd,"data/stator_results/T_cell/md/md_cell_label.csv"))

Idents(seu_obj) <- "barcodes"

seu_obj <- get_state(seu_obj, states)

seu_obj <- subset(seu_obj, idents = md$X, invert = TRUE)

#Log normalize
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

saveRDS(seu_obj, paste0(wd,"data/output_baseline/states/baseline_T_state_valid.rds"))

R_states <- seu_obj@meta.data[seu_obj@meta.data[, "Cluster:1"] == "yes" |
                               seu_obj@meta.data[, "Cluster:2"] == "yes" |
                                 seu_obj@meta.data[, "Cluster:3"] == "yes" |
                               seu_obj@meta.data[, "Cluster:4"] == "yes" ,]

#DE for cluster 3
Idents(seu_obj) <- "Cluster:3"

DE<-FindMarkers(object = seu_obj,ident.1 = "yes",ident.2 = "no",logfc.threshold = 0.25)

DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_",paste0("C3",collapse = "_"))
DE$cluster[DE$avg_log2FC<0]<-paste0("up_",paste0("other",collapse = "_"))

#DE for states enriched in R
Idents(seu_obj) <- "barcodes"

DE<-FindMarkers(object = seu_obj,ident.1 = R_states,logfc.threshold = 0.25)

DE<-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_",paste0("R_states",collapse = "_"))
DE$cluster[DE$avg_log2FC<0]<-paste0("up_",paste0("other",collapse = "_"))

#DE for CD8 T EM cells
Idents(seu_obj) <- "subtype"
seu_obj_EM <- subset(seu_obj, idents = "CD8 T_EM")

table(seu_obj_EM@meta.data$`Cluster:1`, seu_obj_EM@meta.data$label)
table(seu_obj_EM@meta.data$`Cluster:2`, seu_obj_EM@meta.data$label)
table(seu_obj_EM@meta.data$`Cluster:3`, seu_obj_EM@meta.data$label)
table(seu_obj_EM@meta.data$`Cluster:4`, seu_obj_EM@meta.data$label)

#DE for Cluster 2
Idents(seu_obj_EM) <- "Cluster:2"
seu_obj_EM_2 <- subset(seu_obj_EM, idents = "yes")

Idents(seu_obj_EM_2) <- "label"
DE_2_EM<-FindMarkers(object = seu_obj_EM_2,ident.1 = "R",logfc.threshold = 0.25)

DE_2_EM<-DE_2_EM[DE_2_EM$p_val_adj<0.05,]
DE_2_EM$gene<-rownames(DE_2_EM)

DE_2_EM$cluster[DE_2_EM$avg_log2FC>0]<-paste0("up_",paste0("R",collapse = "_"))
DE_2_EM$cluster[DE_2_EM$avg_log2FC<0]<-paste0("up_",paste0("other",collapse = "_"))


