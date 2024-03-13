#Author: Dzhansu Hasanova 10232023
#Cell annotation with Haradhvala et.al.
#Cells were annotated wit cell type, T cells were additionally annotated with subtypes

library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readr)
#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


#Import annotations from Haradhvala et al.
annotations_haradhvala <- read.csv(paste0(wd, "data/metadata/Fig1_CART_global_obs.csv"), row.names=1)

#Import subtype annotations from Haradhvala et al.
annotations_haradhvala_sub <- read_delim("Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# --------------------Baseline---------------------
#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output_baseline/baseline_raw_md.rds"))

#Create column for subtype annotation in baseline samples
annotations_baseline_sub <- annotations_haradhvala_sub[grep("Baseline", annotations_haradhvala_sub$timepoint), ]
annotations_baseline_sub <- annotations_baseline_sub[c("...1", "subtype", "class")]
annotations_baseline_sub$subtype[annotations_baseline_sub$subtype == 'Unknown'] <- 'T cell'
annotations_baseline_sub$class[annotations_baseline_sub$class == 'Unknown'] <- 'T cell'
annotations_baseline_sub$class <- paste0(annotations_baseline_sub$subtype, "_", annotations_baseline_sub$class)


#Subset baseline samples in cell annotation data frame
annotations_baseline <- annotations_haradhvala[grep("Baseline", annotations_haradhvala$timepoint), ]
annotations_baseline$barcodes <- rownames(annotations_baseline)
rownames(annotations_baseline) <- NULL

#Join subtype and cell type annotations, if a cell is annotate in subtype but not in cell type, replace the NA value by the respective annotation in subtype
annotations_baseline <- left_join(annotations_baseline, annotations_baseline_sub, by = join_by(barcodes == ...1))
annotations_baseline$subtype <- ifelse(is.na(annotations_baseline$subtype), annotations_baseline$cell_type, annotations_baseline$subtype)
annotations_baseline$class <- ifelse(is.na(annotations_baseline$class), annotations_baseline$cell_type, annotations_baseline$class)

annotations_baseline <- annotations_baseline[ , c("subtype", "class","barcodes")]

#Crate metadata containing cell annotation
orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcodes <- rownames(seurat_object@meta.data)

annotations_baseline <- annotations_baseline %>%
  filter(annotations_baseline$barcodes %in% orig_rownames)

md <- merge(annotations_baseline, seurat_object@meta.data, by ="barcodes", all.y = TRUE)
md <- md[match(orig_rownames, md$barcodes),]


seurat_object <- AddMetaData(seurat_object, md$subtype , col.name = "cell_type")
seurat_object <- AddMetaData(seurat_object, md$class , col.name = "subtype")

#Add metadata to seurat object
seurat_object@meta.data$cell_type <- replace(seurat_object@meta.data$cell_type, is.na(seurat_object@meta.data$cell_type), "Unknown")
seurat_object@meta.data$subtype <- replace(seurat_object@meta.data$subtype, is.na(seurat_object@meta.data$subtype), "Unknown")

saveRDS(seurat_object, paste0(wd, "data/output_baseline/baseline_raw_md_annot.rds"))

# --------------------Post-treatment---------------------

#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output_posttreatment/post_raw_md.rds"))

#Cell annotation for T cell subset

class <- annotations_haradhvala_sub[c("...1", "subtype", "class")]
class$subtype[class$subtype == 'Unknown'] <- 'T cell'
class$class[class$class == 'Unknown'] <- 'T cell'
class$class <- paste0(class$subtype, "_", class$class)

#rename rows for cell annotation data frame
annotations_haradhvala$barcode <- rownames(annotations_haradhvala)

#Merge df with T cell subtypes with the whole dataset

merged <- left_join(annotations_haradhvala, class, by = join_by(barcode == ...1))
merged$class <- ifelse(is.na(merged$class), merged$cell_type, merged$class)
rownames(merged) <- merged$barcode


orig_rownames <- rownames(seurat_object@meta.data)
seurat_object@meta.data$barcode <- rownames(seurat_object@meta.data)

merged <- merged %>%
  filter(merged$barcode %in% orig_rownames)

md <- merge(merged, seurat_object@meta.data, by ="barcode", all.y = TRUE)
md <- md[match(orig_rownames, md$barcode),]

#Reannotate with subtype information for CD8 and CD4 T cells, keep Infusion T cells in cell_type entry
md$subtype[md$cell_type %in% c("Infusion T")] <- paste0(md$subtype[md$cell_type %in% c("Infusion T")], "_", md$cell_type[md$cell_type %in% c('Infusion T')])
md$cell_type[md$subtype %in% c('CD8 T', 'CD4 T')] <- md$subtype[md$subtype %in% c('CD8 T', 'CD4 T')]

md$class <- replace(md$class, is.na(md$class), "Unknown")
md$cell_type <- replace(md$cell_type, is.na(md$cell_type), "Unknown")

#Add metadata to seurat object
seurat_object <- AddMetaData(seurat_object, md$class , col.name = "cell_subtype")
seurat_object <- AddMetaData(seurat_object, md$cell_type , col.name = "cell_type")

saveRDS(seurat_object, paste0(wd, "data/output_posttreatment/post_raw_md_annot.rds"))


