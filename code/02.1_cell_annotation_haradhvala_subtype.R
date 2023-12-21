#Author: Dzhansu Hasanova 10232023
#Cell annotation, manual and SingleR

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readr)
#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

# --------------------Baseline---------------------
#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output_baseline/baseline_raw_md.rds"))

#Load cell annotation for T cell subset
annotations_haradhvala_sub <- read_delim("Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)
annotations_baseline_sub <- annotations_haradhvala_sub[grep("Baseline", annotations_haradhvala_sub$timepoint), ]
class <- annotations_baseline_sub[c("...1", "subtype", "class")]
class$subtype[class$subtype == 'Unknown'] <- 'T cell'
class$class[class$class == 'Unknown'] <- 'T cell'
class$class <- paste0(class$subtype, "_", class$class)

#Load cell annotation for entire dataset
annotations_haradhvala <- read.csv("Documents/ETH/HS23/data/metadata/Fig1_CART_global_obs.csv", row.names=1)
annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_haradhvala$barcode <- rownames(annotations_haradhvala)

#Merge df with T cell subtypes with the whole dataset
annotations_baseline$barcodes <- rownames(annotations_baseline)

merged <- left_join(annotations_baseline, class, by = join_by(barcodes == ...1))
merged$class <- ifelse(is.na(merged$class), merged$cell_type, merged$class)
rownames(merged) <- merged$barcodes

#To DO: add cell subtypes to md of seurat object

orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcodes <- rownames(seurat_object@meta.data)

merged <- merged %>%
  filter(merged$barcodes %in% orig_rownames)

md <- merge(merged, seurat_object@meta.data, by ="barcodes", all.y = TRUE)
md <- md[match(orig_rownames, md$barcodes),]


seurat_object <- AddMetaData(seurat_object, md$class , col.name = "cell_subtype")

#saveRDS(seurat_object, paste0(wd, "data/output_baseline/baseline_raw_md_annot.rds"))

# --------------------Post-treatment---------------------

#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output_posttreatment/post_raw_md.rds"))

#Load cell annotation for T cell subset
annotations_haradhvala_sub <- read_delim("Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#annotations_baseline_sub <- annotations_haradhvala_sub[grep("Baseline", annotations_haradhvala_sub$timepoint), ]
class <- annotations_haradhvala_sub[c("...1", "subtype", "class")]
class$subtype[class$subtype == 'Unknown'] <- 'T cell'
class$class[class$class == 'Unknown'] <- 'T cell'
class$class <- paste0(class$subtype, "_", class$class)

#Load cell annotation for entire dataset
annotations_haradhvala <- read.csv("Documents/ETH/HS23/data/metadata/Fig1_CART_global_obs.csv", row.names=1)
#annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_haradhvala$barcode <- rownames(annotations_haradhvala)

#Merge df with T cell subtypes with the whole dataset
#annotations_baseline$barcodes <- rownames(annotations_baseline)

merged <- left_join(annotations_haradhvala, class, by = join_by(barcode == ...1))
merged$class <- ifelse(is.na(merged$class), merged$cell_type, merged$class)
rownames(merged) <- merged$barcode

#To DO: add cell subtypes to md of seurat object

orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcode <- rownames(seurat_object@meta.data)

merged <- merged %>%
  filter(merged$barcode %in% orig_rownames)

md <- merge(merged, seurat_object@meta.data, by ="barcode", all.y = TRUE)
md <- md[match(orig_rownames, md$barcode),]

md$class <- replace(md$class, is.na(md$class), "not_defined")
md$cell_type <- replace(md$cell_type, is.na(md$cell_type), "not_defined")

seurat_object <- AddMetaData(seurat_object, md$class , col.name = "cell_subtype")
seurat_object <- AddMetaData(seurat_object, md$cell_type , col.name = "cell_type")

saveRDS(seurat_object, paste0(wd, "data/output_posttreatment/post_raw_md_annot.rds"))



#Get for subset 1 Stator run 1

md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md.csv")
md.x <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_cell_id.csv")

md <- left_join(md, md.x, by = join_by(X == X))
md <- md[,-3]
colnames(md) <- c("X", "label", "orig.ident", "Cell.Types")

md_new <- left_join(md, merged[merged$barcodes %in% md$X,], by = join_by(X == barcodes))
md_new$subtype <- ifelse(is.na(md_new$subtype), md_new$cell_type, md_new$subtype)

md_new[is.na(md_new$class),]$class <- "not_defined"
md_new[is.na(md_new$subtype),]$subtype <- "not_defined"

rownames(md_new) <- md_new$X
md_new <- md_new[,-1]

md_sub <- md_new[c("label", "class")]
colnames(md_sub) <- c("Cell.State", "Cell.Types")

write.csv(md_sub, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_subtype.csv", quote=FALSE)

md_joint <- md_sub
md_joint$Cell.Types <- paste0(md_joint$Cell.State, "_", md_joint$Cell.Types)

write.csv(md_joint, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_subtype_joint.csv", quote=FALSE)


md_corr <- md_new[c("label", "subtype")]
colnames(md_corr) <- c("Cell.State", "Cell.Types")
write.csv(md_corr, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_corr.csv", quote=FALSE)

md_corr_joint <- md_corr
md_corr_joint$Cell.Types <- paste0(md_corr_joint$Cell.State, "_", md_corr_joint$Cell.Types)
write.csv(md_corr_joint, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_corr_joint.csv", quote=FALSE)

md_corr_id_cell <- md_new[c("orig.ident", "subtype")]
colnames(md_corr_id_cell) <- c("Cell.State", "Cell.Types")
write.csv(md_corr_id_cell, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_corr_cell_id.csv", quote=FALSE)







