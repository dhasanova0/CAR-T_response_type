#Author: Dzhansu Hasanova
#Create meta data as input for Shiny App
#Metadata created for cells from post-treatment axi-cel day 7

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


subset1_md <- read.csv(paste0(wd, "data/stator_results/axi_D7_1/subset1_md.csv"), header=TRUE)

rownames(subset1_md) <- subset1_md$X

md_1 <- subset1_md[,c("label", "cell_type")]
write.csv(md_1, paste0(wd, "data/stator_results/axi_D7_1/md_cell_label.csv"))

md_2 <- md_1
md_2$cell_label <- paste0(md_1$cell_type, "_",md_1$label)
md_2 <- md_2[,-2]
write.csv(md_2, paste0(wd,"data/stator_results/axi_D7_1/md_cell_label_joint.csv"))

md_3 <- subset1_md[,c("label", "cell_subtype")]
write.csv(md_3, paste0(wd,"data/stator_results/axi_D7_1/md_subtype_label.csv"))

md_4 <- md_3
md_4$cell_label <- paste0(md_3$cell_subtype, "_",md_3$label)
md_4 <- md_4[,-2]
write.csv(md_4, paste0(wd,"data/stator_results/axi_D7_1/md_subtype_label_joint.csv"))

md_5 <- subset1_md[,c("label", "Sex")]
write.csv(md_5, paste0(wd,"data/stator_results/axi_D7_1/md_sex_label.csv"))

md_6 <- subset1_md[,c("label", "orig.ident")]
write.csv(md_6, paste0(wd,"data/stator_results/axi_D7_1/md_id_label.csv"))

md_7 <- md_6
md_7$id_label <- paste0(md_6$orig.ident, "_",md_6$label)
md_7 <- md_7[,-2]
write.csv(md_7, paste0(wd,"data/stator_results/axi_D7_1/md_id_label_joint.csv"))

