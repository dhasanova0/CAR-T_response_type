library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
source("~/CAR-T_response_type/code/DE_Stator.R")

ls <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/axi_D7_1/Shiny results/CellList-2024-01-04.csv")
ls_0.5 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/axi_D7_1/Shiny results/CellList-0.5.csv")
md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/axi_D7_1/md_cell_label.csv")
seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_posttreatment/stator_input/axi_cel/rds_D7/subset_1.rds")

DefaultAssay(seu_obj) <- "RNA"
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

DE_Infusion_9 <- DE_stator(md, seu_obj,ls, "Infusion T", "9")
DE_Infusion_12 <- DE_stator(md, seu_obj,ls, "Infusion T", "12")

DE_NK_12 <- DE_stator(md, seu_obj,ls, "NK", "12")

DE_B_10 <- DE_stator(md, seu_obj,ls, "B", "10")

DE_Monocyte_3 <- DE_stator(md, seu_obj,ls, "Monocyte", "3")
DE_Monocyte_4 <- DE_stator(md, seu_obj,ls, "Monocyte", "4")
DE_Monocyte_6 <- DE_stator(md, seu_obj,ls, "Monocyte", "6")
DE_Monocyte_7 <- DE_stator(md, seu_obj,ls, "Monocyte", "7")

DE_DC_1 <- DE_stator(md, seu_obj,ls, "mDC", "1")
DE_DC_2 <- DE_stator(md, seu_obj,ls, "mDC", "2")
DE_DC_5 <- DE_stator(md, seu_obj,ls, "mDC", "5")
DE_DC_7 <- DE_stator(md, seu_obj,ls, "mDC", "7")

# Dice distance: 0.5

DE_CD8_54 <- DE_stator(md, seu_obj,ls_0.5, "CD8 T", "54")

DE_CD4_38 <- DE_stator(md, seu_obj,ls_0.5, "CD4 T", "38")
DE_CD4_47 <- DE_stator(md, seu_obj,ls_0.5, "CD4 T", "47")
DE_CD4_43 <- DE_stator(md, seu_obj,ls_0.5, "CD4 T", "43")

write.csv(DE_CD4_43, paste0("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/axi_D7_1/Shiny results/DE0.5_CD4_c43_RvsNR.csv"), quote=FALSE)





colnames(md) <- c("X", "Cell.State", "Cell.Types")
md <- md[md$Cell.Types == "CD8",]

cluster = "54"

c <- ls_0.5[ls_0.5$State == paste0('cluster_C:',cluster), ]
selected_cluster_list <- strsplit(c$CellID, ",\\s*")[[1]]

md_c <- md[md$X%in%selected_cluster_list,]

md_c_NR <- md_c[md_c$Cell.State == "NR",]
md_c_R <- md_c[md_c$Cell.State == "R",]

Idents(seu_obj) <- seu_obj@meta.data$barcode

DE<-FindMarkers(object = seu_obj, ident.1 = md_c_NR$X, ident.2 = md_c_R$X, logfc.threshold = 0.25)
DE <-DE[DE$p_val_adj<0.05,]
DE$gene<-rownames(DE)

DE$cluster[DE$avg_log2FC>0]<-paste0("up_","NR")
DE$cluster[DE$avg_log2FC<0]<-paste0("up_","R")

write.csv(DE, paste0("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/axi_D7_1/Shiny results/DE_CD4_c38_RvsNR.csv"), quote=FALSE)
