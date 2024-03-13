set.seed(1234)

library(dplyr)

wd = "/Users/dhasanova/Documents/ETH/HS23/"

X_linked <- read.delim(paste0(wd, "data/X_linked.txt"))
Y_linked <- read.delim(paste0(wd, "data/Y_linked.txt"))

X_linked$Approved.symbol <- gsub("-", ".", X_linked$Approved.symbol)
Y_linked$Approved.symbol <- gsub("-", ".", Y_linked$Approved.symbol)

subset1 <- read.csv(paste0(wd, "data/output_baseline/stator_input/all_cells/subset1_raw.csv"), header=TRUE)
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

write.csv(subset1_filtered, paste0(wd, "data/output_baseline/stator_input/all_cells/subset1_filtered_raw.csv"), quote=FALSE)

