library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)

#Define paths
wd = "/Users/dhasanova/Documents/ETH/HS23/"
path_fig <- paste0(wd, "figures/")
output_path <- (paste0(wd, "data/output/stator_input/CD8/"))
if (file.exists(output_path)) {cat("The folder already exists")} else {dir.create(output_path)}

# Import Seurat object
input_path <- paste0(wd, "data/output/")
seurat_obj <- readRDS(paste0(input_path, "baseline_doublet_filtered_md_annot.rds"))
seurat_obj_corr <- readRDS(paste0(input_path, "baseline_raw_md_annot_corr.rds"))

#Correct annotation in filtered baseline_doublet_filtered_md_annot.rds
md1 <- data.frame(seurat_obj@meta.data)
md2 <- data.frame(seurat_obj_corr@meta.data)

md_adj <- md2[md2$barcodes %in% md1$barcodes,]
md_adj <- md_adj[match(md1$barcodes, md_adj$barcodes),]

md1$cell_type <- md_adj$cell_type

#Add correct annotation to seurat object

seurat_obj <- AddMetaData(seurat_obj, md1$cell_type , col.name = "cell_type")


# Subset CD8 T cells
seurat_CD8 <- subset(seurat_obj, subset = cell_type == 'CD8 T')

# Get md as metadata and plot number of T cells per patient ID
md <- data.frame(seurat_CD8@meta.data)
md$label_id <- paste0(md$orig.ident, "_", md$label)

counts <- data.frame(table(md$label_id))

ggplot(data=counts, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("CD8_id_dist.png", path = path_fig, height = 7, width = 9)

# Save .csv files (raw counts and md) for Stator input
df <- t(as.data.frame(seurat_CD8@assays[["RNA"]]@counts))
rownames(df) <- sub(".*_", "", rownames(df))

write.csv(df, file=paste0(output_path, "subset_CD8","_raw.csv"), quote=FALSE)
write.csv(md, file=paste0(output_path, "subset_CD8","_md.csv"))


DefaultAssay(seurat_CD8) <- "RNA"
