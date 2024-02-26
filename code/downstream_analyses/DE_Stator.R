DE_stator <- function(md, seurat_object, ls, cell_type, cluster){
  # DE for Monocytes from cluster 7
  colnames(md) <- c("X", "Cell.State", "Cell.Types")
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
  
  return(DE)
}