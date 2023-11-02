#Author: Dzhansu Hasanova 10162023
#Remove Doublets with DoupletFinder for each sample

#Import libraries
library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(DoubletFinder)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

seurat_obj <- readRDS(paste0(wd, "/data/output/baseline_filtered_md_annot.rds"))
seurat_obj_no_rbc <- readRDS(paste0(wd, "data/output/baseline_filtered_no_rbc_md_annot.rds"))

seurat_obj.split <- SplitObject(seurat_obj, split.by = "orig.ident")
seurat_obj_no_rbc.split <- SplitObject(seurat_obj_no_rbc, split.by = "orig.ident")

identify_doublets <- function(seurat_object, create_plots, plot_label){
  
  for (i in 1:length(seurat_object)) {
    # print the sample we are on
    print(paste0("Sample ",i))
    
    pbmc_1 <- SCTransform(seurat_object[[i]])
    pbmc_1 <- RunPCA(object = pbmc_1)
    pbmc_1 <- FindNeighbors(object = pbmc_1, dims = 1:20)
    pbmc_1 <- FindClusters(object = pbmc_1)
    pbmc_1 <- RunUMAP(object = pbmc_1, dims = 1:20)
    
    ## pK Identification (no ground-truth)
    sweep.res.list_pbmc1 <- paramSweep_v3(pbmc_1, PCs = 1:20, sct = TRUE)
    sweep.stats_pbmc1 <- summarizeSweep(sweep.res.list_pbmc1, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats_pbmc1)
    
    pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
      filter(BCmetric == max(BCmetric)) %>%
      select(pK) 
    pK <- as.numeric(as.character(pK[[1]]))
    
    ## Homotypic doublet proportion estimate
    annotations <- pbmc_1@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations) 
    nExp.poi <- round(0.046 * nrow(pbmc_1@meta.data))
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    
    # run doubletFinder 
    seurat_object[[i]] <- doubletFinder_v3(pbmc_1, 
                                              PCs = 1:20, 
                                              pN = 0.25, 
                                              pK = pK, 
                                              nExp = nExp.poi.adj,
                                              reuse.pANN = FALSE, sct = TRUE)
    
    colnames(seurat_object[[i]]@meta.data)[ncol(seurat_object[[i]]@meta.data)-1] <- "pANN"
    colnames(seurat_object[[i]]@meta.data)[ncol(seurat_object[[i]]@meta.data)] <- "doublets"
    
  }
  
  # visualize doublets
  if (create_plots == TRUE){
    for (i in 1:length(seurat_object)){
      
      dimplot <- DimPlot(seurat_object[[i]], reduction = 'umap', group.by = "doublets")
      ggsave(paste0("dimplot_",i,plot_label,".png"), plot = dimplot,  path = "/Users/dhasanova/Documents/ETH/HS23/figures/QC/doublets")
      
      scatter <- FeatureScatter(seurat_object[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",  group.by = "doublets")
      ggsave(paste0("scatter_",i,plot_label,".png"), plot = scatter,  path = "/Users/dhasanova/Documents/ETH/HS23/figures/QC/doublets")
      
    }
    
  }
  
  gc()
  baseline_filtered_doublets <- merge(seurat_object[[1]], y = seurat_object[-c(1)], add.cell.ids = ls(seurat_object), project = "baseline")
  baseline_filtered <- subset(baseline_filtered_doublets, subset = doublets == "Singlet")
  
  return(baseline_filtered)
  
}




baseline_filtered <- identify_doublets(seurat_obj.split, FALSE, "")
baseline_filtered_no_rbc <- identify_doublets(seurat_obj_no_rbc.split,FALSE,  "_norbc" )

if (file.exists(paste0(wd, "data/output/"))) {cat("The folder already exists")} else {dir.create(output_path)}

saveRDS(baseline_filtered, paste0(wd, "data/output/baseline_doublet_filtered_md_annot.rds"))
saveRDS(baseline_filtered_no_rbc, paste0(wd, "data/output/baseline_doublet_no_rbc_filtered_md_annot.rds"))


# Create visualizations
#plot1 <- FeatureScatter(baseline_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  #xlab("nGenes")+
  #ylab("MT%")+
  #ggtitle("")
#plot2 <- FeatureScatter(baseline_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  #xlab("nUMIs")+
  #ylab("nGenes")
#plot1 + theme(legend.position="none") + 
  #plot2 + theme(legend.position="none")
#ggsave("QC_after_doublet.png", path = "/Users/dhasanova/Documents/ETH/HS23/figures/QC/")

#FeatureScatter(baseline_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  #xlab("nUMIs")+
  #ylab("nGenes")+ theme(legend.position="none")
#ggsave(paste0("dimplot_",i,".png"), path = "/Users/dhasanova/Documents/ETH/HS23/figures/QC/")

#table(baseline_filtered@meta.data$doublets)

#counts <- data.frame(table(baseline_filtered@meta.data$orig.ident))
#cat("sum", sum(counts$Freq))
#cat("\nmean", mean(counts$Freq))
#cat("\nmax", max(counts$Freq))
#cat("\nmin", min(counts$Freq))



