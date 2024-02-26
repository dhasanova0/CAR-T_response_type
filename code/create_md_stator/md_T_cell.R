#Author: Dzhansu Hasanova
#Create meta data as input for Stator
#Metadata created for T cells from baseline

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


seu_obj <- readRDS(paste0(wd, "data/output_baseline/stator_input_T/rds/subset_1.rds"))

md <- seu_obj@meta.data

md_label <- md[, c("label", "cell_type")]
write.csv(md_label, paste0(wd,"data/stator_results/T_cell/md/md_cell_label.csv"))

md_subtype <- md[, c("label", "subtype")]
write.csv(md_subtype, paste0(wd,"data/stator_results/T_cell/md/md_cell_subtype.csv"))

md_joint_subtype <- md_subtype
md_joint_subtype$subtype <- paste0(md_joint_subtype$subtype, "_", md_joint_subtype$label)
write.csv(md_joint_subtype, paste0(wd,"data/stator_results/T_cell/md/md_subtype_label_joint.csv"))

md_joint <- md_label
md_joint$cell_type <- paste0(md_joint$cell_type, "_", md_joint$label)
write.csv(md_joint, paste0(wd,"data/stator_results/T_cell/md/md_cell_label_joint.csv"))


md_ID <- md[, c("label", "orig.ident")]
md_ID$orig.ident <- paste0(md_ID$label, "_", md_ID$orig.ident)
write.csv(md_ID, paste0(wd,"data/stator_results/T_cell/md/md_ID.csv"))


md_asct <- md[, c("label", "autologous_transplant")]
write.csv(md_asct, paste0(wd,"data/stator_results/T_cell/md/md_asct.csv"))


md_product <- md[, c("label", "product")]
write.csv(md_product, paste0(wd,"data/stator_results/T_cell/md/md_product.csv"))


md_sex <- md[, c("label", "sex")]
md_sex$sex <- paste0(md_sex$label, "_", md_sex$sex)
write.csv(md_sex, paste0(wd,"data/stator_results/T_cell/md/md_sex.csv"))


md_bridging <- md[, c("no_proior_therapy", "bridging_therapy")]
write.csv(md_bridging, paste0(wd,"data/stator_results/T_cell/md/md_therapy.csv"))
