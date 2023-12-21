library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)

ls <- read.csv("/Users/dhasanova/Documents/ETH/HS23/figures/stator_run1/CellList-2023-11-23.csv")
ls_0.5 <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/run1/CellList-2023-12-08_0.5.csv")
md <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_subtype.csv")
md_corr <- read.csv("/Users/dhasanova/Documents/ETH/HS23/code/MFIs/data/md/subset1_md_corr.csv")
seu_obj <- readRDS("/Users/dhasanova/Documents/ETH/HS23/data/output/stator_input/rds_new/subset_1.rds")

DefaultAssay(seu_obj) <- "RNA"
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)

DE_stator <- function(md, seurat_object, ls, cell_type, cluster){
  # DE for Monocytes from cluster 7
  md <- md[md$Cell.Types == cell_type,]
  
  
  c <- ls[ls$State == paste0('cluster_C:',cluster), ]
  selected_cluster_list <- strsplit(c$CellID, ",\\s*")[[1]]
  
  md_c <- md[md$X%in%selected_cluster_list,]
  
  md_c_NR <- md_c[md_c$Cell.State == "NR",]
  md_c_R <- md_c[md_c$Cell.State == "R",]

  
  DE<-FindMarkers(object = seurat_object,ident.1 = md_c_NR$X,ident.2 = md_c_R$X,logfc.threshold = 0.25)
  DE <-DE[DE$p_val_adj<0.05,]
  DE$gene<-rownames(DE)
  
  DE$cluster[DE$avg_log2FC>0]<-paste0("up_","NR")
  DE$cluster[DE$avg_log2FC<0]<-paste0("up_","R")
  
  #write.csv(DE, paste0("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/run1/DE/DE_",cell_type, cluster, "_R_vs_NR_raw.csv"), quote=FALSE)
  
  return(md_c)
}

DE_Monocytes <- DE_stator(md, seu_obj,ls, "Monocyte", "7") #R:11, NR:75

DE_mDC <- DE_stator(md, seu_obj,ls, "mDC", "6") #R:44, NR:86, no DE genes

DE_CD8_EMRA <- DE_stator(md, seu_obj,ls, "CD8 T_EMRA", "1") #R:190, NR:119

DE_CD8_EM_3 <- DE_stator(md, seu_obj, ls, "CD8 T_EM", "3") #R:338, NR:87

DE_CD8_EM_4 <- DE_stator(md, seu_obj, ls, "CD8 T_EM", "4") #R:148, NR:21

DE_T_EMRA_2 <- DE_stator(md, seu_obj, ls, "T cell_EMRA", "2") #R:77, NR:3
#DE_T_EMRA_3 <- DE_stator(md, seu_obj, "T cell_EMRA", "3") #R:36, NR:2

DE_CD4_CTL_1 <- DE_stator(md, seu_obj, ls, "CD4 T_CD4+ CTL", "1") #R:25 ,NR:95
DE_CD4_CTL_2 <- DE_stator(md, seu_obj, ls, "CD4 T_CD4+ CTL", "2") #R:29 ,NR:108
DE_CD4_CTL_3 <- DE_stator(md, seu_obj, ls, "CD4 T_CD4+ CTL", "3") #R:5 ,NR:54

DE_T_CTL_3 <- DE_stator(md, seu_obj,ls, "T cell_CD4+ CTL", "3") #R:11, NR:34
#DE_T_CTL_1 <- DE_stator(md, seu_obj, "T cell_CD4+ CTL", "4") #R:6, NR:2 #fewer than 3 cells doesn't work


DE_CD4_Naive_1 <- DE_stator(md, seu_obj, ls,"CD4 T_Naive", "1") #R:6, NR:19

DE_CD4_EMRA_1 <- DE_stator(md, seu_obj,ls, "CD4 T_EMRA", "1") #R:3, NR:4

DE_CD4 <- DE_stator(md_corr, seu_obj, ls_0.5 ,"CD4 T", "5")

md_seurat <- seu_obj@meta.data[seu_obj@meta.data$barcodes %in% DE_CD8_EM_3$X,]
md_seurat$label_id <- paste0(md_seurat$orig.ident, "_" ,md_seurat$sex,"_", md_seurat$label)

counts <- data.frame(table(md_seurat$label_id))

counts[c('Var1', 'sex', "label")] <- str_split_fixed(counts$Var1, '_', 3)
counts$Var1 <- paste0(counts$Var1, "_", counts$label)

ggplot(data=counts, aes(x=reorder(Var1, -Freq), y=Freq, fill=sex)) +
  geom_bar(stat="identity")+
  geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  ggtitle("CD8 EM Cluster 3")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
