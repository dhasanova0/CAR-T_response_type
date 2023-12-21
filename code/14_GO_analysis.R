library("org.Hs.eg.db",character.only = T)

DE <- read.csv("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/run1/DE_CD8c1_R_vs_NR_raw.csv")

gene_list_NR <- DE$gene[DE$avg_log2FC>0]
gene_list_R <- DE$gene[DE$avg_log2FC<=0]

GO <- function(gene_list, GO_aspect, responder){
  
  entrez <- mapIds(org.Hs.eg.db, gene_list, "ENTREZID", "SYMBOL")
  cluster_GOBP <- enrichGO(gene = entrez,
                           OrgDb= "org.Hs.eg.db",
                           ont = GO_aspect,
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
  
  
  pdf(paste0("/Users/dhasanova/Documents/ETH/HS23/figures/stator_run1/GO_CD8c1_",GO_aspect,"_",responder,".pdf"), width=11, height=9)
  print(barplot(cluster_GOBP, showCategory=10))
  dev.off()
  
}

GO(gene_list_NR, "BP", "NR")
GO(gene_list_R, "BP", "R")

GO(gene_list_NR, "CC", "NR")
GO(gene_list_R, "CC", "R")

GO(gene_list_NR, "MF", "NR")
GO(gene_list_R, "MF", "R")






