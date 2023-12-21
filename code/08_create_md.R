library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)


md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_corr.csv")
wd = "/Users/dhasanova/Documents/ETH/HS23/"

input_path1 <- paste0(wd, "data/output/stator_input/rds_new/")
seurat_obj <- readRDS(paste0(input_path1, "subset_1.rds"))

df <- as.data.frame(seurat_obj@meta.data)

#Sex

df_sex <- df[, c("barcodes","sex")]
df_sex <- df_sex[order(match(rownames(df_sex),md$X)),]
md_sex <- md
md_sex$Cell.Types <- paste0(df_sex$sex, "_", md$Cell.State)

rownames(md_sex) <- md_sex$X
md_sex <- md_sex[,-1]


write.csv(md_sex, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_sex.csv", quote=FALSE)

#patient ID

df_id <- df[, c("no_proior_therapy", "bridging_therapy")]
#df_id <- df_id[order(match(rownames(df_id),md$X)),]
df_id$id_label <- paste0(df_id$orig.ident, "_", df_id$label)
df_id <- df_id[, -2]
colnames(df_id) <- c("Cell.State", "Cell.Types")

write.csv(df_id, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_id.csv", quote=FALSE)

unique(df_id$Cell.Types)

#cell and R/NR
df_joint <- md
df_joint$Cell.Types <- paste0(md$Cell.State, "_", md$Cell.Types)
rownames(df_joint) <- df_joint$X
df_joint <- df_joint[,-1]

write.csv(df_joint, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_joint.csv", quote=FALSE)

#cell and patient joint

df_cell_id <- df[, c("label", "orig.ident", "cell_type")]
df_cell_id$Cell.Types <- paste0(df_cell_id$orig.ident, "_", df_cell_id$cell_type)
df_cell_id <- df_cell_id[order(match(rownames(df_cell_id),md$X)),]
df_cell_id <- df_cell_id[, -c(2,3)]
colnames(df_cell_id) <- c("Cell.State", "Cell.Types")

write.csv(df_cell_id, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_cell_id_joint.csv", quote=FALSE)

#cell and patient

df_cell_id <- df[, c("orig.ident", "cell_type")]
df_cell_id <- replace(df_cell_id, is.na(df_cell_id), "not_defined")
df_cell_id <- df_cell_id[order(match(rownames(df_cell_id),md$X)),]
colnames(df_cell_id) <- c("Cell.State", "Cell.Types")

write.csv(df_cell_id, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_cell_id.csv", quote=FALSE)


#Bridging

df_bridging <- df[, c("label", "bridging_therapy")]
colnames(df_bridging) <- c("Cell.State", "Cell.Types")
write.csv(df_bridging, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_bridging.csv", quote=FALSE)

#Number therapies

df_no_therapy <- df[, c("label", "no_proior_therapy")]
colnames(df_no_therapy) <- c("Cell.State", "Cell.Types")
write.csv(df_no_therapy, "/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_no_therapy.csv", quote=FALSE)






