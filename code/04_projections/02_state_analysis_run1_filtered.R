#Author: Dzhansu Hasanova
#Projection of states from baseline subset 1 filtered for patient specific genes onto the remaining cell
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

#Import Data
seu_obj <- readRDS(paste0(wd, "data/output_baseline/states/base_filtered_states.rds"))


md_corr <- read.csv(paste0(wd, "data/stator_results/run1/md/subset1_md_cell.csv"))

states_0.83 <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-filtered.csv"))
states_0.5 <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-2024-0.5.csv"))
states_0.5 <- states_0.5[match(states_0.83$X, states_0.5$X),]

#Filter

#BCR filtering
seu_obj <- seu_obj[!grepl("^IG[HKL]V", rownames(seu_obj)), ]
seu_obj <- seu_obj[!grepl("^IG[HKL]J", rownames(seu_obj)), ]
seu_obj <- seu_obj[!grepl("^IG[HKL]C", rownames(seu_obj)), ]
seu_obj <- seu_obj[!grepl("^IGH[ADEGM]", rownames(seu_obj)), ]

#TCR filtering
seu_obj <- seu_obj[!grepl("^TR[ABDG][VJC]", rownames(seu_obj)), ]

#MHC filetring
seu_obj <- seu_obj[!grepl("^HLA", rownames(seu_obj)), ]


X_linked <- read.delim(paste0(wd, "data/X_linked.txt"))
Y_linked <- read.delim(paste0(wd, "data/Y_linked.txt"))

X_linked$Approved.symbol <- gsub("-", ".", X_linked$Approved.symbol)
Y_linked$Approved.symbol <- gsub("-", ".", Y_linked$Approved.symbol)

#Filter X and Y-linked genes
seu_obj <- seu_obj[is.na(!match(rownames(seu_obj), X_linked$Approved.symbol)), ]
seu_obj <- seu_obj[is.na(!match(rownames(seu_obj), Y_linked$Approved.symbol)), ]


#Subset cells on which Stator was not run
Idents(seu_obj) <- "barcodes"
seu_obj <- subset(seu_obj, idents = md_corr$X, invert = TRUE)

#Log normalize
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

#saveRDS(seu_obj, paste0(wd, "data/output_baseline/counts/base_filtered_states_valid.rds"))
#seu_obj <- readRDS(paste0(wd, "data/output_baseline/counts/base_filtered_states_valid.rds"))

seu_obj_0.83 <- get_state(seu_obj, states_0.83)
seu_obj_0.5 <- get_state(seu_obj, states_0.5)

#Perform DE based on different enrichment between R and NR

Idents(seu_obj_0.5) <- "Cluster:27"
CD8_27 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "CD8 T")

Idents(CD8_27) <- "label"
DE_CD8_7 <- FindMarkers(object = CD8_27, ident.1 = "R", logfc.threshold = 0.25)
DE_CD8_7 <- DE_CD8_7[DE_CD8_7$p_val_adj<0.05,]
DE_CD8_7$genes <- rownames(DE_CD8_7)
write.csv(DE_CD8_7, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_CD8_27.csv"))

Idents(seu_obj_0.5) <- "Cluster:22"
CD8_22 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "CD8 T")

Idents(CD8_22) <- "label"
DE_CD8_22 <- FindMarkers(object = CD8_22, ident.1 = "R", logfc.threshold = 0.25)
DE_CD8_22 <- DE_CD8_22[DE_CD8_22$p_val_adj<0.05,]
DE_CD8_22$genes <- rownames(DE_CD8_22)
write.csv(DE_CD8_22, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_CD8_22.csv"))

Idents(seu_obj_0.5) <- "Cluster:23"
CD8_23 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "CD8 T")

Idents(CD8_23) <- "label"
DE_CD8_23 <- FindMarkers(object = CD8_23, ident.1 = "R", logfc.threshold = 0.25)
DE_CD8_23 <- DE_CD8_23[DE_CD8_23$p_val_adj<0.05,]
DE_CD8_23$genes <- rownames(DE_CD8_23)
write.csv(DE_CD8_23, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_CD8_23.csv"))

Idents(seu_obj_0.5) <- "Cluster:26"
mono_26 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "Monocyte")

Idents(mono_26) <- "label"
DE_mono_26 <- FindMarkers(object = mono_26, ident.1 = "R", logfc.threshold = 0.25)
DE_mono_26 <- DE_mono_26[DE_mono_26$p_val_adj<0.05,]
DE_mono_26$genes <- rownames(DE_mono_26)
write.csv(DE_mono_26, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_mono_26.csv"))

Idents(seu_obj_0.5) <- "Cluster:24"
NK_24 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "NK")

Idents(NK_24) <- "label"
DE_NK_24 <- FindMarkers(object = NK_24, ident.1 = "R", logfc.threshold = 0.25)
DE_NK_24 <- DE_NK_24[DE_NK_24$p_val_adj<0.05,]
DE_NK_24$genes <- rownames(DE_NK_24)
write.csv(DE_NK_24, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_NK_24.csv"))

Idents(seu_obj_0.5) <- "Cluster:21"
NK_21 <- subset(x = seu_obj_0.5, idents = "yes", subset = cell_type == "NK")

Idents(NK_21) <- "label"
DE_NK_21 <- FindMarkers(object = NK_21, ident.1 = "R", logfc.threshold = 0.25)
DE_NK_21 <- DE_NK_21[DE_NK_21$p_val_adj<0.05,]
DE_NK_21$genes <- rownames(DE_NK_21)
write.csv(DE_NK_21, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_R_vsNR_NK_21.csv"))


#Cluster vs all other cells from the cell type

CD8 <- subset(x = seu_obj_0.5, subset = cell_type == "CD8 T")

Idents(CD8) <- "Cluster:27"

DE_CD8_7 <- FindMarkers(object = CD8, ident.1 = "yes", logfc.threshold = 0.25)
DE_CD8_7 <- DE_CD8_7[DE_CD8_7$p_val_adj<0.05,]
DE_CD8_7$genes <- rownames(DE_CD8_7)
write.csv(DE_CD8_7, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_CD8_27_all.csv"))

Idents(CD8) <- "Cluster:22"
DE_CD8_22 <- FindMarkers(object = CD8, ident.1 = "yes", logfc.threshold = 0.25)
DE_CD8_22 <- DE_CD8_22[DE_CD8_22$p_val_adj<0.05,]
DE_CD8_22$genes <- rownames(DE_CD8_22)
write.csv(DE_CD8_22, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_CD8_22_all.csv"))

Idents(CD8) <- "Cluster:23"

DE_CD8_23 <- FindMarkers(object = CD8, ident.1 = "yes", logfc.threshold = 0.25)
DE_CD8_23 <- DE_CD8_23[DE_CD8_23$p_val_adj<0.05,]
DE_CD8_23$genes <- rownames(DE_CD8_23)
write.csv(DE_CD8_23, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_CD8_23_all.csv"))

mono <- subset(x = seu_obj_0.5, subset = cell_type == "Monocyte")
Idents(mono) <- "Cluster:26"

DE_mono_26 <- FindMarkers(object = mono, ident.1 = "yes", logfc.threshold = 0.25)
DE_mono_26 <- DE_mono_26[DE_mono_26$p_val_adj<0.05,]
DE_mono_26$genes <- rownames(DE_mono_26)
write.csv(DE_mono_26, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_mono_26_all.csv"))

NK <- subset(x = seu_obj_0.5, subset = cell_type == "NK")

Idents(NK) <- "Cluster:24"

DE_NK_24 <- FindMarkers(object = NK, ident.1 = "yes", logfc.threshold = 0.25)
DE_NK_24 <- DE_NK_24[DE_NK_24$p_val_adj<0.05,]
DE_NK_24$genes <- rownames(DE_NK_24)
write.csv(DE_NK_24, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_NK_24_all.csv"))

Idents(NK) <- "Cluster:21"

DE_NK_21 <- FindMarkers(object = NK, ident.1 = "yes", logfc.threshold = 0.25)
DE_NK_21 <- DE_NK_21[DE_NK_21$p_val_adj<0.05,]
DE_NK_21$genes <- rownames(DE_NK_21)
write.csv(DE_NK_21, paste0(wd, "data/stator_results/run1_filtered/projected/DE/DE_NK_21_all.csv"))









