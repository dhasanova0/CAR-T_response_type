library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readr)

wd = "/Users/dhasanova/Documents/ETH/HS23/"

annotations_haradhvala_sub <- read_delim("/Users/dhasanova/Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

subset_CD8_raw_18755 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/subset_CD8_raw_18755.csv")
trainData <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/trainingData_18755Cells_1000Genes.csv")

md_CD8 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/CD8_md_joint.csv")

subset_CD8_sc <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/md_scanpy.csv")
subset_CD8_md <- read.csv("~/Documents/ETH/HS23/data/output_baseline/stator_input/CD8/subset_CD8_md.csv")


subset_CD8_md$X <- gsub(".*_","",subset_CD8_md$X)
subset_CD8_md <- subset_CD8_md[subset_CD8_md$X %in% subset_CD8_sc$X,]
subset_CD8_md <- subset_CD8_md[match(subset_CD8_sc$X, subset_CD8_md$X),]

subtype <- annotations_haradhvala_sub[annotations_haradhvala_sub$...1 %in% subset_CD8_md$X,]
subtype <- subtype[,c("...1", "class")]

subset_CD8_md_sub <- left_join(subset_CD8_md, subtype, by = join_by(barcodes == ...1))
subset_CD8_md_sub$class <- replace(subset_CD8_md_sub$class, is.na(subset_CD8_md_sub$class), "CD8 T")

md_label_joint <- subset_CD8_md[, c("X","label", "label_id")]
rownames(md_label_joint) <- md_label_joint$X
md_label_joint <- md_label_joint[,-1]

md_subtype <- subset_CD8_md_sub[, c("X","label", "class")]
rownames(md_subtype) <- md_subtype$X
md_subtype <- md_subtype[,-1]
write.csv(md_subtype, "/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/CD8_md_subtype.csv")

md_sex <- subset_CD8_md_sub[, c("X","label", "sex")]
rownames(md_sex) <- md_sex$X
md_sex <- md_sex[,-1]
write.csv(md_sex, "/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/CD8_md_sex.csv")

subset_CD8_md_sub$sex_id <- paste0(subset_CD8_md_sub$orig.ident,"_", subset_CD8_md_sub$sex)
md_sex_joint <- subset_CD8_md_sub[, c("X","label", "sex_id")]
rownames(md_sex_joint) <- md_sex_joint$X
md_sex_joint <- md_sex_joint[,-1]
write.csv(md_sex_joint, "/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/CD8_md_sex_joint.csv")


subset_CD8_raw <- subset_CD8_raw[subset_CD8_raw$X %in% subset_CD8_sc$X,]

subset_CD8_raw_18755 <- subset_CD8_raw_18755[match(subset_CD8_sc$X, subset_CD8_raw_18755$X),]
rownames(subset_CD8_raw_18755) <- subset_CD8_raw_18755$X
subset_CD8_raw_18755 <- subset_CD8_raw_18755[,-1]





write.csv(md_label_joint, "/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/CD8_md_joint.csv")
write.csv(subset_CD8_raw_18755, "/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/subset_CD8_raw_18755.csv", quote=FALSE)
