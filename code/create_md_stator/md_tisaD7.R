#Author: Dzhansu Hasanova
#Create meta data as input for Shiny App
#Metadata created for cells from post-treatment tisa-cel day 7

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


seu_obj <- readRDS(paste0(wd, "data/output_baseline/counts/post_state_baseline.rds"))
md <- read.csv(paste0(wd,"data/output_posttreatment/stator_input/tisa_cel/raw_md_D7/subset1_md.csv"))

md <- seu_obj@meta.data[rownames(seu_obj@meta.data) %in% md$X,]

md_cell <- md[, c("Response_classification", "cell_type")]

write.csv(md_cell, paste0(wd, "data/stator_results/tisa_D7/md_cell.csv"))


md$ID_label <- paste0(md$Response_classification, "_", md$orig.ident)

md_ID <- md[, c("Response_classification", "ID_label")]
write.csv(md_ID, paste0(wd, "data/stator_results/tisa_D7/md_ID.csv"))

md$cell_label <- paste0(md$Response_classification, "_", md$cell_type)
md_joint <- md[, c("Response_classification", "ID_label")]
write.csv(md_joint, paste0(wd, "data/stator_results/tisa_D7/md_celljoint.csv"))



