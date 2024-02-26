#Author: Dzhansu Hasanova
#Create meta data as input for Shiny App
#Metadata created for cells from baseline (first subset)

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seu_obj <- readRDS(paste0(wd, "data/output_baseline/stator_input/rds_new/subset_1.rds"))

md_seu <- seu_obj@meta.data

md_cell <- md_seu[c("label", "cell_type")]
write.csv(md_cell, paste0(wd, "data/stator_results/run1/md/subset1_md_cell.csv"), quote=FALSE)

md_joint <- md_cell
md_joint$label_cell <- paste0(md_cell$label, "_", md_cell$cell_type)
write.csv(md_joint, paste0(wd,"data/stator_results/run1/md/subset1_md_cell_joint.csv"), quote=FALSE)


md_sub <- md_seu[c("label", "subtype")]
write.csv(md_sub, paste0(wd,"data/stator_results/run1/md/subset1_md_subtype.csv"), quote=FALSE)

md_seu$label_id <- paste0(md_seu$label, "_", md_seu$orig.ident)
md_id <- md_seu[c("label", "label_id")]
write.csv(md_id, paste0(wd,"data/stator_results/run1/md/subset1_md_id.csv"), quote=FALSE)

md_asct <- md_seu[, c("label", "autologous_transplant")]
write.csv(md_asct, paste0(wd,"data/stator_results/run1/md/subset1_md_asct.csv"))


md_product <- md_seu[, c("label", "product")]
write.csv(md_product, paste0(wd,"data/stator_results/run1/md/subset1_md_product.csv"))


md_sex <- md_seu[, c("label", "sex")]
md_sex$sex <- paste0(md_sex$label, "_", md_sex$sex)
write.csv(md_sex, paste0(wd,"stator_results/run1/md/subset1_md_sex.csv"))


md_bridging <- md_seu[, c("no_proior_therapy", "bridging_therapy")]
write.csv(md_bridging,paste0(wd,"data/stator_results/run1/md/subset1_md_bridging.csv"))







