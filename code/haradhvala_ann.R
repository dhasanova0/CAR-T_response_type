library(readr)

md_scanpy <- read.csv("~/Documents/ETH/HS23/data/haradhvala_scanpy/md_scanpy.csv")
md_fig <- read_delim("Documents/ETH/HS23/data/metadata/Fig3_ALLT_knnsubtypes_annotated_obs.txt", 
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_corr.csv")



md_scanpy <- md_scanpy[grep("Baseline", md_scanpy$timepoint), ]

table(md_scanpy$subtype)
table(md_scanpy$cell_type)


md_CD8 <- md_scanpy[md_scanpy$X %in% md[grep("CD8 T", md$Cell.Types), ]$X,]

table(md[grep("CD8 T", md$Cell.Types), ]$X %in% md_scanpy$X)

table(md_scanpy$X %in% md[grep("CD8 T", md$Cell.Types), ]$X)
