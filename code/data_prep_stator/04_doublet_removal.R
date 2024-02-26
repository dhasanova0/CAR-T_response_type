#Author: Dzhansu Hasanova 10162023
#Remove Doublets with DoupletFinder for each sample

#Import libraries
library(Seurat)
library(SeuratObject)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(DoubletFinder)
library(Matrix)

source("~/CAR-T_response_type/code/qc_figures.R")

set.seed(12342)

#Specify path to save figures
path_qc_fig_post <- paste0(wd, "figures/posttreatment/QC/")
if (file.exists(path_qc_fig_post)) {cat("The folder already exists")} else {dir.create(path_qc_fig_post)}

path_qc_fig_baseline <- paste0(wd, "figures/baseline/QC/")
if (file.exists(path_qc_fig_baseline)) {cat("The folder already exists")} else {dir.create(path_qc_fig_baseline)}

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import seurat objects
seurat_obj <- readRDS(paste0(wd, "/data/output_baseline/baseline_filtered_md_annot.rds"))
seurat_obj_post <- readRDS(paste0(wd, "data/output_posttreatment/post_filtered_md_annot.rds"))


#Splot Seurat objects for doublet removal
seurat_obj.split <- SplitObject(seurat_obj, split.by = "orig.ident")
seurat_obj_post.split <- SplitObject(seurat_obj_post, split.by = "orig.ident")



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
    sweep.res.list_pbmc1 <- paramSweep(pbmc_1, PCs = 1:20, sct = TRUE)
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
    seurat_object[[i]] <- doubletFinder(pbmc_1, 
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
      ggsave(paste0("dimplot_",i,plot_label,".png"), plot = dimplot,  path = "/Users/dhasanova/Documents/ETH/HS23/figures/posttreatment/QC/doublets")
      
      scatter <- FeatureScatter(seurat_object[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",  group.by = "doublets")
      ggsave(paste0("scatter_",i,plot_label,".png"), plot = scatter,  path = "/Users/dhasanova/Documents/ETH/HS23/figures/posttreatment/QC/doublets")
      
    }
    
  }
  
  gc()
  baseline_filtered_doublets <- merge(seurat_object[[1]], y = seurat_object[-c(1)], add.cell.ids = ls(seurat_object), project = "baseline")
  baseline_filtered <- subset(baseline_filtered_doublets, subset = doublets == "Singlet")
  
  DefaultAssay(baseline_filtered) <- "RNA"
  baseline_filtered[['SCT']] <- NULL
  
  return(baseline_filtered)
  
}



#------------------- Baseline -------------------

baseline_filtered <- identify_doublets(seurat_obj.split, FALSE, "")
simple_qc_plot(baseline_filtered, "QC_after_doublet.png", path_qc_fig_baseline)
get_counts(baseline_filtered@meta.data$orig.ident,"Number of cells after filtering", "after_filtering_doublet.png", path_qc_fig_baseline)
QC_plots(baseline_filtered@meta.data, "after_doublet", path_qc_fig_baseline)

saveRDS(baseline_filtered, paste0(wd, "data/output/baseline_doublet_filtered_md_annot.rds"))

#------------------- Post-treatment -------------------

seurat_obj_post <- identify_doublets(seurat_obj_post.split, FALSE,  "" )
simple_qc_plot(seurat_obj_post, "QC_after_doublet.png", path_qc_fig_post )
get_counts(seurat_obj_post@meta.data$orig.ident,"Number of cells after filtering", "after_filtering_doublet.png", path_qc_fig_post)
QC_plots(seurat_obj_post@meta.data, "after_doublet", path_qc_fig_post)

saveRDS(baseline_filtered_post, "/exports/igmm/eddie/ponting-lab/dzhansu/data/output/post_doublet_filtered_md_annot.rds")




