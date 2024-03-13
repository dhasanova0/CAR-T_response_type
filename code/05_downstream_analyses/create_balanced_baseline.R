#Author: Dzhansu Hasanova
#This script creates seurat objects with balanced cells across cell types between R and NR
#and balanced cells across patients

set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(DESeq2)
library(groupdata2)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seu_obj_base <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_valid.rds"))

# --------------------- Balance across patients ---------------------

min(table(seu_obj_base@meta.data$orig.ident))

md <- balance(seu_obj_base@meta.data,
              'min',
              "orig.ident")

Idents(seu_obj_base) <- 'barcodes'

seu_obj_balanced_id <- subset(seu_obj_base, idents = md$barcodes)

saveRDS(seu_obj_balanced_id, paste0(wd, "data/output_baseline/states/baseline_state_balanced_id.rds"))

# --------------------- Balance R and NR across cell types ---------------------
monocytes <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "Monocyte" & seu_obj_base@meta.data$label == "NR",], 16931 - 9591)$barcodes
CD8 <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "CD8 T" & seu_obj_base@meta.data$label == "R",], 9175 - 6212)$barcodes
CD4 <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "CD4 T" & seu_obj_base@meta.data$label == "NR",], 6206 - 3486)$barcodes
NK <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "NK" & seu_obj_base@meta.data$label == "NR",], 4689 - 4064)$barcodes
pDC <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "pDC" & seu_obj_base@meta.data$label == "NR",], 43 - 21)$barcodes
T_cell <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "T cell" & seu_obj_base@meta.data$label == "R",], 2319 - 2237)$barcodes
Platelet <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "Platelet" & seu_obj_base@meta.data$label == "NR",], 530 - 237)$barcodes
Unknown <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "Unknown" & seu_obj_base@meta.data$label == "NR",], 1512 - 990)$barcodes
mDC <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "mDC" & seu_obj_base@meta.data$label == "NR",], 634 - 462)$barcodes
B <- sample_n(seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "B" & seu_obj_base@meta.data$label == "R",], 712 - 17)$barcodes
Infusion <- seu_obj_base@meta.data[seu_obj_base@meta.data$cell_type == "Infusion T",]$barcodes

Idents(seu_obj_base) <- 'barcodes'

seu_obj_balanced <- subset(seu_obj_base, idents = c(monocytes, CD8, CD4, NK, pDC, Platelet, Unknown, mDC, B, T_cell, Infusion), invert = TRUE)

saveRDS(seu_obj_balanced, paste0(wd, "data/output_baseline/states/baseline_state_balanced_cell.rds"))


