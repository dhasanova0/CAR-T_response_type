set.seed(1234)

library(dplyr)

X_linked <- read.delim("~/Documents/ETH/HS23/data/X_linked.txt")
Y_linked <- read.delim("~/Documents/ETH/HS23/data/Y_linked.txt")

X_linked$Approved.symbol <- gsub("-", ".", X_linked$Approved.symbol)
Y_linked$Approved.symbol <- gsub("-", ".", Y_linked$Approved.symbol)

subset1 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/stator_input/all_cells/subset1_raw.csv", header=TRUE)
genes <- colnames(subset1)


subset1_filtered <- subset1[,!(colnames(subset1) %in% Y_linked$Approved.symbol)]
subset1_filtered <- subset1_filtered[,!(colnames(subset1_filtered) %in% X_linked$Approved.symbol)]


subset1_filtered <- subset1_filtered[, !grepl("^IG[HKL]V", colnames(subset1_filtered))]
subset1_filtered <- subset1_filtered[, !grepl("^IG[HKL]J", colnames(subset1_filtered))]
subset1_filtered <- subset1_filtered[, !grepl("^IG[HKL]C", colnames(subset1_filtered))]
subset1_filtered <- subset1_filtered[, !grepl("^IGH[ADEGM]", colnames(subset1_filtered))]

#TCR filtering
subset1_filtered <- subset1_filtered[,!grepl("^TR[ABDG][VJC]", colnames(subset1_filtered))]

#MHC filetring
subset1_filtered <- subset1_filtered[, !grepl("^HLA", colnames(subset1_filtered))]

rownames(subset1_filtered) <- subset1_filtered$X

subset1_filtered <- subset1_filtered[, -1]

write.csv(subset1_filtered, "/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/stator_input/all_cells/subset1_filtered_raw_2.csv", quote=FALSE)

md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_corr.csv")

table(md$Cell.Types)

monocytes <- sample_n(md[md$Cell.Types == "Monocyte",], 5075)$X
CD8 <- sample_n(md[md$Cell.Types == "CD8 T",], 2100)$X
CD4 <- sample_n(md[md$Cell.Types == "CD4 T",], 161)$X

md_filtered <- md[!md$X %in% c(monocytes, CD8, CD4),]

balanced_cells <- subset1[subset1$X %in% md_filtered$X,]

rownames(balanced_cells) <- balanced_cells$X

balanced_cells <- balanced_cells[, -1]

write.csv(balanced_cells, "/Users/dhasanova/Documents/ETH/HS23/data/output_baseline/stator_input/all_cells/subset1_cell_balanced.csv", quote=FALSE)


