library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(pheatmap)


md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/subset1_md_corr_cell_id.csv")

g1 <- as.data.frame(table(md[md$Cell.State %in% c("Patient6-Baseline", "Patient13-Baseline", "Patient17-Baseline" ),]$Cell.Types))
g2 <- as.data.frame(table(md[md$Cell.State %in% c("Patient24-Baseline", "Patient18-Baseline", 
                                                  "Patient19-Baseline", "Patient9-Baseline", "Patient8-Baseline",
                                                  "Patient14-Baseline", "Patient20-Baseline"),]$Cell.Types))
g3 <- as.data.frame(table(md[md$Cell.State %in% c("Patient14-Baseline", "Patient20-Baseline", "Patient12-Baseline",
                                                  "Patient25-Baseline", "Patient21-Baseline", "Patient30-Baseline", 
                                                  "Patient31-Baseline"),]$Cell.Types))
g4 <- as.data.frame(table(md[md$Cell.State %in% c("Patient30-Baseline", "Patient31-Baseline", "Patient15-Baseline",
                                                  "Patient23-Baseline", "Patient22-Baseline", "Patient10-Baseline", "Patient11-Baseline" ),]$Cell.Types))
g5 <- as.data.frame(table(md[md$Cell.State %in% c("Patient10-Baseline", "Patient11-Baseline"),]$Cell.Types))

groups <- list(g1, g2, g3, g4, g5)


r=5
col=length(unique(md$Cell.Types))
Data_mtrix<-array(data=NA,dim = c(r,col))
colnames(Data_mtrix) <- unique(md$Cell.Types)
rownames(Data_mtrix) <- c("group1", "group2", "group3", "group4", "group5")

P_set=rep(NA,length(unique(md$Cell.Types)))
names(P_set)<-unique(md$Cell.Types)

fold_set=rep(NA,length(unique(md$Cell.Types)))
names(fold_set)<-unique(md$Cell.Types)


for (i in 1:length(groups)){
  
  for (cell in groups[[i]]$Var1){
    
    print(cell)
    q = groups[[i]][groups[[i]]$Var1 == cell,]$Freq
    m = sum(md$Cell.Types == cell)
    n = nrow(md) - m
    k = sum(groups[[i]]$Freq)
    p_value <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
    fold_value <- (nrow(md)*(m-q))/(m*(nrow(md)-k))
    P_set[cell] <- p_value
    fold_set[cell] <- fold_value
    
    for (diff in setdiff(unique(md$Cell.Types), groups[[i]]$Var1)){
      P_set[diff] <- NA
    }
    
  }
  Data_mtrix[i,] <- P_set
}


#q = g3[g3$Var1 == "Monocyte",]$Freq
#m = sum(md$Cell.Types == "Monocyte")
#n = nrow(md) - m
#k = sum(g3$Freq)

#p_value <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)


Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
log10FDR<--log10(Mydata_raw_m)


log10FDR[log10FDR==Inf]<-400
colnames(log10FDR) <- unique(md$Cell.Types)
rownames(log10FDR) <- c("group1", "group2", "group3", "group4", "group5")

p_label <- Mydata_raw_m
p_label[p_label <= 0.05] <- "*"
p_label[p_label > 0.05] <- ""

pdf("/Users/dhasanova/Documents/ETH/HS23/figures/stator_run1/heatmap_cell_id_corr.pdf", width=6, height=5)
ComplexHeatmap::pheatmap(log10FDR, display_numbers = p_label, fontsize_number=15, cellheight=20, 
                         cluster_rows = FALSE,
                         cluster_cols = FALSE, heatmap_legend_param=list(title_position = "lefttop-rot"),
                         name = "-Log10FDR", color = (colorRampPalette(c("white", "red", "darkred"))(50)), 
                         main = "Cell type enrichment in patient groups")
dev.off()








