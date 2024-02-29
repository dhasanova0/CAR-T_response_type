#Author: Dzhansu Hasanova
#Comparison states inferred from baseline samples projected onto post-treatment
#and baseline cells

set.seed(123)

library(ggstatsplot)
library(dplyr)
library(tidyr)
library(Seurat)
library(DESeq2)
library(ggsignif)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seu_obj_post <- readRDS(paste0(wd, "data/output_baseline/states/post_state_baseline.rds"))
seu_obj_base <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state.rds"))

states <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))

seu_obj_post <- get_state(seu_obj_post, states)
seu_obj_base <- get_state(seu_obj_base, states)

get_state_DE <- function(seu_obj, ident){
  
  #DE analyses to compare between pos-treatment samples
  
  Idents(seu_obj) <- ident
  DE <- FindMarkers(object = seu_obj, ident.1 = "yes", logfc.threshold = 0.25)
  DE <- DE[DE$p_val_adj<0.05,]

  return(DE)
  
}

DE_post_1 <- get_state_DE(seu_obj_post, "Cluster:1")
DE_post_3 <- get_state_DE(seu_obj_post, "Cluster:3")
DE_post_4 <- get_state_DE(seu_obj_post, "Cluster:4")
DE_post_5 <- get_state_DE(seu_obj_post, "Cluster:5")

DE_base_1 <- get_state_DE(seu_obj_base, "Cluster:1")
DE_base_3 <- get_state_DE(seu_obj_base, "Cluster:3")
DE_base_4 <- get_state_DE(seu_obj_base, "Cluster:4")
DE_base_5 <- get_state_DE(seu_obj_base, "Cluster:5")

intersect(rownames(DE_post_1[DE_post_1$avg_log2FC>0,]), rownames(DE_base_1[DE_base_1$avg_log2FC>0,]))
intersect(rownames(DE_post_3[DE_post_3$avg_log2FC>0,]), rownames(DE_base_1[DE_base_3$avg_log2FC>0,]))
intersect(rownames(DE_post_4[DE_post_4$avg_log2FC>0,]), rownames(DE_base_1[DE_base_4$avg_log2FC>0,]))
intersect(rownames(DE_post_5[DE_post_5$avg_log2FC>0,]), rownames(DE_base_1[DE_base_5$avg_log2FC>0,]))

#Merge the seu objects to see differences in the cells in states enriched  differently in R and NR

colnames(seu_obj_post@meta.data)[colnames(seu_obj_post@meta.data) == "Response_classification"] <- "label"

seu_obj_base@meta.data$timepoint <- "baseline"
seu_obj_post@meta.data$timepoint <- "post"

column_names <- intersect(colnames(seu_obj_base@meta.data), colnames(seu_obj_post@meta.data))

seu_obj_base@meta.data <- seu_obj_base@meta.data[,column_names ]
seu_obj_post@meta.data <- seu_obj_post@meta.data[,column_names ]

colnames(seu_obj_base@meta.data)
colnames(seu_obj_post@meta.data)

merged <- merge(seu_obj_base, y = seu_obj_post,  project = "timepoints")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)

merged@meta.data$barcodes <- rownames(merged@meta.data)

get_DE_timepoint <- function(seu_obj, cluster){
  
  cells_state_c1_base <- rownames(seu_obj@meta.data[seu_obj@meta.data[, cluster] == "yes" & seu_obj@meta.data[, "timepoint"] == "baseline" ,])
  cells_state_c1_post <- rownames(seu_obj@meta.data[seu_obj@meta.data[, cluster] == "yes" & seu_obj@meta.data[, "timepoint"] == "post" ,])
  
  Idents(seu_obj) <- "barcodes"
  DE <- FindMarkers(object = seu_obj, ident.1 = cells_state_c1_base, ident.2 = cells_state_c1_post, logfc.threshold = 0.25)
  DE <- DE[DE$p_val_adj<0.05,]
  return(DE)

}

DE_1 <- get_DE_timepoint(merged, "Cluster.1")
DE_3 <- get_DE_timepoint(merged, "Cluster.3")
DE_4 <- get_DE_timepoint(merged, "Cluster.4")
DE_5 <- get_DE_timepoint(merged, "Cluster.5")

#Check for exhaustion markers in T cells in state 1 NR cells
Idents(merged) <- "cell_type"

unique(merged@meta.data$cell_type)

T_c1_NR_merged <- subset(merged, idents = c("CD8 T", "T cell", "CD4 T", "Infusion T"), subset = Cluster.1 == "yes")


Idents(T_c1_NR_merged) <- "timepoint"
DE <- FindMarkers(object = T_c1_NR_merged, ident.1 = "baseline", logfc.threshold = 0.25)
DE <- DE[DE$p_val_adj<0.05,]

#Compare CD69 expression between the two timepoints
Idents(merged) <- "timepoint"
VlnPlot(merged, features = "CD69", pt.size = 0, cols = c("dodgerblue3", "darkorange2")) + NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(paired = FALSE)
ggsave("violin_CD69_timepoints.png", path = paste0(wd, "figures/posttreatment/", width = 7, height = 5))


