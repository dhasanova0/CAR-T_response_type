#Author: Dzhansu Hasanova 10042023
#First inspection of baseline samples

#Import libraries
library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)

#Get Data
data_dir <- '/Users/dhasanova/Documents/ETH/HS23/data/GSE197268_RAW/baseline'
baseline_samples <- list.dirs(data_dir)
baseline_samples <- baseline_samples[-c(1)]

#Create lists to save names of samples and seurat objects
names <- c()
object_list <- c()

#Create a seurat object for each sample
for (i in baseline_samples){
  data <- Read10X(data.dir = i)
  object_list <- c(object_list, CreateSeuratObject(counts = data, project = str_sub(i, 65)))
  names <- c(names, str_sub(i, 65))
}

# Merge seurat objects of baseline samples
baseline_seurat <- merge(object_list[[1]], y = object_list[-c(1)], add.cell.ids = names, project = "baseline")
rm(object_list)

gc()

#--------------------QC--------------------

# Create functions for figures

get_counts <- function(orig_ident, plot_title, plot_name, path_name){
  counts <- data.frame(table(orig_ident))
  cat("sum", sum(counts$Freq))
  cat("\nmean", mean(counts$Freq))
  cat("\nmax", max(counts$Freq))
  cat("\nmin", min(counts$Freq))

  ggplot(data=counts, aes(x=orig_ident, y=Freq)) +
    geom_bar(stat="identity", fill = "deepskyblue4")+
    ylab("Number of cells")+xlab("")+
    theme_bw()+
    ggtitle(plot_title)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  ggsave(plot_name, path = path_name)
}

qc_plot <- function(seurat_object, plot_name, plot_path){
  
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")+
    xlab("nGenes")+
    ylab("MT%")+
    ggtitle("")
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    xlab("nUMIs")+
    ylab("nGenes")
  plot1 + theme(legend.position="none") + 
    plot2 + theme(legend.position="none")
  ggsave(plot_name, path = plot_path)
}


#Create plots for QC
QC_plots <- function(df, name, plot_path){
  df %>% 
    ggplot(aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density")+
    xlab("nGenes")
  ggsave(paste0("Gene_",name,".png"), path = plot_path, width=10, height=7)
  
  df %>% 
    ggplot(aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density")+
    xlab("nUMIs")
  ggsave(paste0("UMI_",name,".png"), path = plot_path, width=10, height=7)
  
  
  df %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
    geom_point()+
    stat_smooth(method=lm) +
    scale_colour_gradient(low = "gray90", high = "black") +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_hline(yintercept = 200) +
    facet_wrap(~orig.ident)+
    ylab("nGenes")+
    xlab("nUMIs")
  ggsave(paste0("joint_mt_",name,".png"), path = plot_path, width=13, height=10)
  
  df %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=rbc)) + 
    geom_point()+ 
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_hline(yintercept = 200) +
    facet_wrap(~orig.ident)+
    ylab("nGenes")+
    xlab("nUMIs")
  ggsave(paste0("joint_rb_",name,".png"), path = plot_path, width=13, height=10)
  
}


#Define path for figures
path_qc_fig <- "/Users/dhasanova/Documents/ETH/HS23/figures/QC"

get_counts(baseline_seurat@meta.data$orig.ident,"Number of cells before filtering", "before_filtering.png", path_qc_fig)

#Get mitochondiral genes
baseline_seurat[["percent.mt"]] <- PercentageFeatureSet(baseline_seurat, pattern = "^MT-")
qc_plot(baseline_seurat, "QC_before.png", path_qc_fig)
QC_plots(baseline_seurat@meta.data, "before", path_qc_fig)


VlnPlot(baseline_seurat, features = "nFeature_RNA", pt.size = 0) + NoLegend() + 
  geom_hline(yintercept = 200)+
  ylab("Number of Genes")+ggtitle("Genes")
ggsave("nGenes.png", path = path_qc_fig, width=13, height=9)

#Get RBC
baseline_seurat[["rbc"]] <- colnames(baseline_seurat) %in% 
  colnames(subset(x = baseline_seurat, subset = HBB > 0 | HBA1 > 0 | HBA2 > 0))


#Check samples with high RBC number
VlnPlot(baseline_seurat, features = c("HBB"), pt.size=0, log = TRUE)
rbc <- subset(baseline_seurat, subset = HBB > 0 | HBA1 > 0 | HBA2 > 0)
rbc_fraction <- data.frame(table(rbc@meta.data$orig.ident) / table(baseline_seurat@meta.data$orig.ident))

ggplot(data=rbc_fraction, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  ylab("Fraction")+xlab("")+
  theme_bw()+
  ggtitle("RBC fraction")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("rbc_fraction_patient.png", path = path_qc_fig)

#Create figures for rbc filtering only with the given thresholds obtained from
#median expressed of haemoglobin genes per patient
baseline_seurat_rbc <- subset(baseline_seurat, HBB < 5 & HBA1 < 1 & HBA2 < 1)
QC_plots(baseline_seurat_rbc@meta.data, "after_bc", path_qc_fig)

#Perform filtering only MT and RBC to determine the threshold for gene filtering
baseline_seurat_filtered <- subset(baseline_seurat, subset = percent.mt < 15 & HBB < 5 & HBA1 < 1 & HBA2 < 1)
VlnPlot(baseline_seurat_filtered, features = "nFeature_RNA", pt.size = 0) + NoLegend() + 
  geom_hline(yintercept = 200)+
  ylab("Number of Genes")+ggtitle("Genes")
ggsave("nGenes_filtered_only_RBC_MT.png", path = path_qc_fig, width=13, height=9)

#Perform complete filtering
baseline_seurat_filtered <- subset(baseline_seurat, subset = nFeature_RNA > 200  &  percent.mt < 15 & HBB < 5 & HBA1 < 1 & HBA2 < 1)
baseline_seurat_filtered_remove_rbc <- subset(baseline_seurat, subset = nFeature_RNA > 200  &  percent.mt < 15 & HBB < 1 & HBA1 < 1 & HBA2 < 1)

#Create plots after filtering
qc_plot(baseline_seurat_filtered, "QC_after.png", path_qc_fig )
get_counts(baseline_seurat_filtered@meta.data$orig.ident,"Number of cells after filtering", "after_filtering.png", path_qc_fig)
QC_plots(baseline_seurat_filtered@meta.data, "after", path_qc_fig)

#Remove all cells with hemoglobin genes
baseline_seurat_filtered_remove_rbc <- subset(baseline_seurat_filtered, subset = HBB < 1 & HBA1 < 1 & HBA2 < 1)

#Save Seuarat objects with filtering: 1. filter for certain threshold of hemoglobin genes, 2. Filter out all cells expressing hemoglobin genes

saveRDS(baseline_seurat_filtered, "/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_filtered.rds")
saveRDS(baseline_seurat_filtered_remove_rbc, "/Users/dhasanova/Documents/ETH/HS23/data/output/baseline_filtered_no_rbc.rds")


#Additional plots
#baseline_seurat_filtered@meta.data %>% 
  #ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=rbc)) + 
  #geom_point(alpha = 0.2)+ 
  #scale_x_log10() + 
  #scale_y_log10() + 
  #theme_classic() +
  #ylab("nGenes")+
  #xlab("nUMIs")+
  #geom_hline(yintercept = 200)
#ggsave(paste0("joint_rb_","after",".png"), path = path_qc_fig, width=7, height=7)


#baseline_seurat_filtered@meta.data %>% 
  #ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  #geom_point()+
  #scale_colour_gradient(low = "gray90", high = "black") +
  #scale_x_log10() + 
  #scale_y_log10() + 
  #theme_classic() +
  #ylab("nGenes")+
  #xlab("nUMIs")+
  #geom_hline(yintercept = 200)
#ggsave(paste0("joint_mt_","after",".png"), path = path_qc_fig, width=7, height=7)















