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


seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/counts/baseline_state.rds")



seu_obj@meta.data$orig.ident_label_product <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$label, "_", seu_obj@meta.data$product)
seu_obj@meta.data$cell_label <- paste0(seu_obj@meta.data$cell_type, "_", seu_obj@meta.data$label)
seu_obj@meta.data$orig.ident_label <- paste0(seu_obj@meta.data$orig.ident, "_", seu_obj@meta.data$label)

cell_id <- data.frame(table(seu_obj@meta.data[, c("cell_label", "orig.ident_label")]))
cell_id$orig.ident_label <- as.character(cell_id$orig.ident_label)

cell_id$label_from_cell <- sub(".*_", "", cell_id$cell_label)
cell_id$label_from_id <- sub(".*_", "", cell_id$orig.ident_label)
cell_id <- cell_id[cell_id$label_from_id==cell_id$label_from_cell,]

cell <- data.frame(table(seu_obj@meta.data[, c("cell_type", "orig.ident_label")]))
cell$orig.ident_label <- as.character(cell$orig.ident_label)
cell$label <- sub('.*_', '', cell$orig.ident_label)

find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

cell_id <- cell_id %>%
  group_by(cell_label) %>%
  mutate(outlier = ifelse(find_outlier(Freq), orig.ident_label, NA))

cell_id$outlier_2 <- NA

cell_id$outlier_2[cell_id$outlier == 'Patient24-Baseline_NR' & cell_id$cell_label == 'Monocyte_NR' ] <- 'Patient24-Baseline_NR'
cell_id$outlier_2[cell_id$outlier == 'Patient13-Baseline_R' & cell_id$cell_label == 'CD8 T_R'] <- 'Patient13-Baseline_R'


ggplot(cell, aes(x=cell_type, y=Freq)) + 
  geom_boxplot(outlier.colour="darkred", outlier.shape=8, outlier.size = 2)+
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(label)), size = 1.5, alpha = 0.5)+
  scale_color_manual(name = "Response Type", values=c("deepskyblue4", "#E69F00"))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )+
  ylab("Counts")+xlab("")
ggsave("boxplot_cell.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 8, height = 5)

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
ggsave("boxplot_cell.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 10, height = 5)


# remove outlier Patient24-Baseline_NR
Idents(seu_obj) <- "orig.ident"
seu_obj_filtered <- subset(seu_obj, idents = "Patient24-Baseline", invert = TRUE)



#saveRDS(seu_obj_filtered, "/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/baseline_doublet_filtered_md_annot_norm_rmoutlier.rds")

#Perform normalization with eddie
#seu_obj_filtered <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/baseline_norm_no_outlier.rds")

cell_filtered <- data.frame(table(seu_obj_filtered@meta.data[, c("cell_type", "orig.ident_label")]))
cell_filtered$orig.ident_label <- as.character(cell_filtered$orig.ident_label)
cell_filtered$label <- sub('.*_', '', cell_filtered$orig.ident_label)

ggplot(cell_filtered, aes(x=cell_type, y=Freq)) + 
  geom_boxplot(outlier.colour="darkred", outlier.shape=8, outlier.size = 2)+
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(color=factor(label)), size = 1.5, alpha = 0.5)+
  scale_color_manual(name = "Response Type", values=c("deepskyblue4", "#E69F00"))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )+
  ylab("Counts")+xlab("")
ggsave("boxplot_cell_filtered.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 8, height = 5)



DimPlot(seu_obj, group.by = "cell_type", reduction = "umap", cols = cols1) +ggtitle("Cell Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_celltype.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 6, height = 5)

DimPlot(seu_obj, group.by = "label", reduction = "umap", cols = cols1) +ggtitle("Response Type") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_label.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 6, height = 5)

DimPlot(seu_obj, group.by = "seurat_clusters", reduction = "umap", label = "TRUE", cols = cols1) + ggtitle("Seurat Cluster") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("cluster.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 6, height = 5)

DimPlot(seu_obj, group.by = "orig.ident", reduction = "umap", cols = cols1) + ggtitle("Patient ID") + ylab("UMAP 2") + xlab("UMAP 1")
ggsave("umap_patient.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 6, height = 5)


#Perform DE for R and NR
Idents(seu_obj) <- "label"
DE <- FindMarkers(seu_obj, ident.1 = "R")
DE <- DE[DE$p_val_adj<0.05,]

write.csv(DE, "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/DE_RvsNR_log.csv")

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
ggsave("DE_volcano.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 8, height = 7)



Idents(seu_obj_filtered) <- "label"
DE_filtered <- FindMarkers(seu_obj_filtered, ident.1 = "R")
DE_filtered <- DE_filtered[DE_filtered$p_val_adj<0.05,]


EnhancedVolcano(DE_filtered , 
                rownames(DE_filtered ),
                x ="avg_log2FC", 
                y ="p_val_adj",
                FCcutoff = 0.5,
                labSize = 4.0)
ggsave("DE_volcano_filtered.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 6, height = 7)


seu_R <- subset(seu_obj, subset = label == "R")

Idents(seu_R) <- "autologous_transplant"
DE_asct <- FindMarkers(seu_R, ident.1 = "yes")
DE_asct <- DE_asct[DE_asct$p_val_adj<0.05,]

write.csv(DE_asct, "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/DE_asct_log.csv")

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
ggsave("DE_volcano_asct.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/baseline/analysis/", width = 8, height = 7)


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



