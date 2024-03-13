#Author: Dzhansu Hasanova
#script adapted from https://github.com/YuelinYao/MFIs/blob/main/app/tabs/heatmap/heatmap.R
#Create heatmap of enrichment analysis
set.seed(100)

library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(pheatmap)


#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

enrich_analysis <- function(ls, col_names, md, row_names, title, path ){
  #Input: ls: list of df for which enrichment is performed
  #       col_names: Column names for heatmap
  #       md: data for which enrichment is performed
  #       row_name: Row names of heatmap
  #       title: Title of heatmap
  #Output: Heatmap of enrichment analysis
  
  #Prepare matrices for enrichment analysis
  r=length(ls)
  col=length(unique(col_names))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  colnames(Data_mtrix) <- c(unique(col_names))
  rownames(Data_mtrix) <- row_names
  
  Fold<-array(data=NA,dim = c(r,col))
  colnames(Fold) <- c(unique(col_names))
  rownames(Fold) <- row_names
  
  P_set=rep(NA,length(unique(col_names)))
  names(P_set)<-unique(col_names)
  
  fold_set=rep(NA,length(unique(col_names)))
  names(fold_set)<-unique(col_names)
  
  #Perform enrichment analysis
  for (i in 1:length(ls)){
    
    for (cell in ls[[i]]$Var1){
      
      print(cell)
      q = as.numeric(ls[[i]][ls[[i]]$Var1 == cell,]$Freq) #Number of cell of interest in group
      m = as.numeric(sum(col_names == cell)) #Total number of cell of interest in entire data
      n = as.numeric(nrow(md) - m) #Number of cell without cell of interest
      k = as.numeric(sum(ls[[i]]$Freq)) #Number of cells in group
      p_value <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
      fold_value <- (q*(m+n))/(m*k)
      P_set[cell] <- p_value
      fold_set[cell] <- fold_value
      
      for (diff in setdiff(unique(col_names), ls[[i]]$Var1)){
        P_set[diff] <- NA
        fold_set[diff] <- NA
      }
      
    }
    Data_mtrix[i,] <- P_set
    Fold[i,] <- fold_set
  }
  
  #Adjust p-values and get -log10
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)
  
  
  colnames(log10FDR) <- unique(col_names)
  rownames(log10FDR) <- row_names
  
  #Significance
  p_label <- Mydata_raw_m
  p_label[p_label <= 0.05] <- "*"
  p_label[p_label > 0.05] <- ""
  
  
  ht <- ComplexHeatmap::pheatmap(log10FDR, display_numbers = p_label, fontsize_number=10, cellwidth = 15,cellheight = 15, 
                                 cluster_rows = FALSE,show_row_dend = FALSE, border_color = NA,fontsize = 12,
                                 cluster_cols = FALSE, heatmap_legend_param=list(title_position = "lefttop-rot"),
                                 name = "-Log10FDR", color = colorRampPalette(c("white","firebrick3"))(100), 
                                 main = title)
  write.csv(Fold, path)
  return(ht)
  
}

#Prepare data

#Baseline
md <- read.csv(paste0(wd, "code/MFIs/data/md/subset1_md_corr_cell_id.csv"))
md$Cell.Types[md$Cell.Types == "not_defined"] <- "Unknown"
md_label <- read.csv(paste0(wd,"data/stator_results/run1/md/subset1_md_label_ID_joint.csv"))
md$Cell.State <- md_label$label_ID

seu_obj <- readRDS(paste0(wd,"data/output_baseline/states/baseline_state.rds"))
md_seurat <- data.frame(seu_obj@meta.data)

R <- as.data.frame(table(md_seurat[md_seurat$label %in% c("R"),]$cell_type))
NR <- as.data.frame(table(md_seurat[md_seurat$label %in% c("NR"),]$cell_type))

groups2 <- list(R, NR)

R_sub <- as.data.frame(table(md_seurat[md_seurat$label %in% c("R"),]$subtype))
NR_sub <- as.data.frame(table(md_seurat[md_seurat$label %in% c("NR"),]$subtype))

groups_subtype <- list(R_sub, NR_sub)

#Patient subset 1 baseline
patient <- c("NR_Patient6-Baseline", "R_Patient13-Baseline", "R_Patient17-Baseline", "NR_Patient24-Baseline",
             "R_Patient19-Baseline", "NR_Patient18-Baseline", "R_Patient8-Baseline", "NR_Patient9-Baseline",
             "NR_Patient14-Baseline", "NR_Patient20-Baseline", "R_Patient12-Baseline","R_Patient21-Baseline", 
             "NR_Patient25-Baseline","R_Patient30-Baseline", "NR_Patient31-Baseline", "R_Patient15-Baseline",
             "NR_Patient23-Baseline", "R_Patient22-Baseline", "R_Patient11-Baseline","NR_Patient10-Baseline")
patient_l <- c()

for (i in 1:20){
  print(i)
  df <- as.data.frame(table(md[md$Cell.State %in% patient[i],]$Cell.Types))
  patient_l[[i]] <- df
}


#Posttreatment
seu_obj_post <- readRDS(paste0(wd,"data/output_baseline/states/post_state_baseline.rds"))
md_seurat_post <- data.frame(seu_obj_post@meta.data)

R_post <- as.data.frame(table(md_seurat_post[md_seurat_post$Response_classification %in% c("R"),]$cell_type))
NR_post <- as.data.frame(table(md_seurat_post[md_seurat_post$Response_classification %in% c("NR"),]$cell_type))

groups2_post <- list(R_post, NR_post)


#Create heatmaps
pdf(paste0(wd,"figures/baseline/analysis/heatmap_cell_enrichment.pdf"), width=6, height=5)
ht <- enrich_analysis(groups2, md_seurat$cell_type, md_seurat, c("R", "NR"), 
                "Cell type enrichment in R and NR", paste0(wd,"figures/baseline/analysis/Fold_cell.csv"))
draw(ht,heatmap_legend_side = "left")
dev.off()

pdf(paste0(wd,"figures/baseline/analysis/heatmap_subtype_enrichment.pdf"), width=10, height=5)
ht <- enrich_analysis(groups_subtype, md_seurat$subtype, md_seurat, c("R", "NR"), 
                      "Cell subtype enrichment in R and NR", paste0(wd,"figures/baseline/analysis/Fold_subtype.csv"))
draw(ht,heatmap_legend_side = "left")
dev.off()

pdf(paste0(wd,"figures/stator_run1/thesis/heatmap_cell_enrichment_patient_label.pdf"), width=7, height=10)
ht <- enrich_analysis(patient_l, md$Cell.Types, md, patient, 
                      "Cell type enrichment in sample", paste0(wd,"figures/stator_run1/thesis/Fold_cell.csv"))
draw(ht,heatmap_legend_side = "left")
dev.off()

pdf(paste0(wd,"figures/posttreatment/cell_enrichment_post.pdf"), width=7, height=10)
ht <- enrich_analysis(groups2_post, md_seurat_post$cell_type, md_seurat_post, c("R", "NR"), 
                      "Cell type enrichment in post-treatment",paste0(wd,"figures/posttreatment/Fold_cell.csv" ))
draw(ht,heatmap_legend_side = "left")
dev.off()
