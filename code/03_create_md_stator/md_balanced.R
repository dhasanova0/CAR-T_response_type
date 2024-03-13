#Author: Dzhansu Hasanova
#Create meta data as input for Shiny App
#Metadata created for cells from balanced baseline (first subset)

library(Seurat)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


md <- read.csv(paste0(wd, "data/stator_results/run1/md/subset1_md_cell.csv"))
balanced <- read.csv(paste0(wd, "data/stator_results/balanced/subset1_cell_balanced.csv"))$X

md <- md[md$X %in% balanced,]
rownames(md) <- md$X
md <- md[,-1]
write.csv(md, paste0(wd, "data/stator_results/balanced/md_cell.csv"))

md$Cell.Types <- paste0(md$Cell.State, "_", md$Cell.Types)
write.csv(md, paste0(wd, "data/stator_results/balanced/md_cell_joint.csv"))


