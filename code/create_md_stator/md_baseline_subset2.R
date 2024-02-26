#Author: Dzhansu Hasanova
#Create meta data as input for Shiny App
#Metadata created for cells from baseline (second subset)

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


seu_obj <- readRDS(paste0(wd, "data/output_baseline/stator_input/rds_new/subset_2.rds"))

md_s2 <- seu_obj@meta.data

md_s2_label <- md_s2[, c("label", "cell_type")]
write.csv(md_s2_label, paste0(wd, "data/stator_results/run3/Shiny/md_cell_label.csv"))

md_s2_joint <- md_s2_label
md_s2_joint$cell_type <- paste0(md_s2_joint$cell_type, "_", md_s2_joint$label)
write.csv(md_s2_joint, paste0(wd, "data/stator_results/run3/Shiny/md_cell_label_joint.csv"))


md_s2_ID <- md_s2[, c("label", "orig.ident")]
md_s2_ID$orig.ident <- paste0(md_s2_ID$label, "_", md_s2_ID$orig.ident)
write.csv(md_s2_ID, paste0(wd, "data/stator_results/run3/Shiny/md_ID.csv"))


md_s2_asct <- md_s2[, c("label", "autologous_transplant")]
write.csv(md_s2_asct, paste0(wd, "data/stator_results/run3/Shiny/md_asct.csv"))


md_s2_product <- md_s2[, c("label", "product")]
write.csv(md_s2_product, paste0(wd, "data/stator_results/run3/Shiny/md_product.csv"))


md_s2_sex <- md_s2[, c("label", "sex")]
md_s2_sex$sex <- paste0(md_s2_sex$label, "_", md_s2_sex$sex)
write.csv(md_s2_sex, paste0(wd, "data/stator_results/run3/Shiny/md_sex.csv"))


md_s2_bridging <- md_s2[, c("no_proior_therapy", "bridging_therapy")]
write.csv(md_s2_bridging, paste0(wd, "data/stator_results/run3/Shiny/md_therapy.csv"))
