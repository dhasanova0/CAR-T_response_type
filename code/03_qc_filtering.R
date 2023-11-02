#Author: Dzhansu Hasanova
#Perform QC, filtering of cells

#Import libraries
library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import Seurat Object

seurat_object <- readRDS(paste0(wd, "data/output/baseline_raw_md_annot.rds"))


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

simple_qc_plot <- function(seurat_object, plot_name, plot_path){
  
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

#Specify path to save figures
path_qc_fig <- paste0(wd, "figures/QC_updated/")
if (file.exists(path_qc_fig)) {cat("The folder already exists")} else {dir.create(path_qc_fig)}

#--------------------Perform QC--------------------

get_counts(seurat_object@meta.data$orig.ident,"Number of cells before filtering", "before_filtering.png", path_qc_fig)

VlnPlot(seurat_object, features = "nFeature_RNA", pt.size = 0) + NoLegend() + 
  geom_hline(yintercept = 200)+
  ylab("Number of Genes")+ggtitle("Genes")
ggsave("nGenes.png", path = path_qc_fig, width=13, height=9)

#Get mitochondiral genes
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
simple_qc_plot(seurat_object, "QC_before.png", path_qc_fig)
QC_plots(seurat_object@meta.data, "before", path_qc_fig)

#Get RBC
seurat_object[["rbc"]] <- colnames(seurat_object) %in% 
  colnames(subset(x = seurat_object, subset = HBB > 0 | HBA1 > 0 | HBA2 > 0))

#Check samples with high RBC number
VlnPlot(seurat_object, features = c("HBB"), pt.size=0, log = TRUE)
rbc <- subset(seurat_object, subset = HBB > 0 | HBA1 > 0 | HBA2 > 0)
rbc_fraction <- data.frame(table(rbc@meta.data$orig.ident) / table(seurat_object@meta.data$orig.ident))

ggplot(data=rbc_fraction, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  ylab("Fraction")+xlab("")+
  theme_bw()+
  ggtitle("RBC fraction")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("rbc_fraction_patient.png", path = path_qc_fig)

#Create figures for rbc filtering only with the given thresholds obtained from
#median expressed of haemoglobin genes per patient

# 1. Find thresholds
hbb_exp <- FetchData(seurat_object, vars = c("HBB", "orig.ident"))
HBB_thresh <- round(mean(aggregate(HBB ~ orig.ident, data = hbb_exp, FUN = median)$HBB)) 

hba1_exp <- FetchData(seurat_object, vars = c("HBA1", "orig.ident"))
HBA1_thresh <- round(mean(aggregate(HBA1 ~ orig.ident, data = hba1_exp, FUN = median)$HBA1))

hba2_exp <- FetchData(seurat_object, vars = c("HBA2", "orig.ident"))
HBA2_thresh <- round(mean(aggregate(HBA2 ~ orig.ident, data = hba2_exp, FUN = median)$HBA2))

baseline_seurat_rbc <- subset(seurat_object, HBB <= HBB_thresh & HBA1 <= HBA1_thresh & HBA2 <= HBA2_thresh)
QC_plots(baseline_seurat_rbc@meta.data, "after_bc", path_qc_fig)

#Perform filtering only MT and RBC to determine the threshold for gene filtering
baseline_seurat_filtered <- subset(seurat_object, subset = percent.mt < 15 & 
                                     HBB <= HBB_thresh & HBA1 <= HBA1_thresh & HBA2 <= HBA2_thresh)
VlnPlot(baseline_seurat_filtered, features = "nFeature_RNA", pt.size = 0) + NoLegend() + 
  geom_hline(yintercept = 200)+
  ylab("Number of Genes")+ggtitle("Genes")
ggsave("nGenes_filtered_only_RBC_MT.png", path = path_qc_fig, width=13, height=9)

#Perform complete filtering 
#Create two Seurat objects with rbc filtered with the thresholds 
#and with rbc completely removed according to the three hem genes

baseline_seurat_filtered <- subset(seurat_object, subset = nFeature_RNA > 200  &  percent.mt < 15
                                   & HBB <= HBB_thresh & HBA1 <= HBA1_thresh & HBA2 <= HBA2_thresh)
baseline_seurat_filtered_remove_rbc <- subset(seurat_object, subset = nFeature_RNA > 200  &  percent.mt < 15 & HBB < 1 & HBA1 < 1 & HBA2 < 1)

#Create plots after filtering
simple_qc_plot(baseline_seurat_filtered, "QC_after.png", path_qc_fig )
get_counts(baseline_seurat_filtered@meta.data$orig.ident,"Number of cells after filtering", "after_filtering.png", path_qc_fig)
QC_plots(baseline_seurat_filtered@meta.data, "after", path_qc_fig)


#Save Seuarat objects with filtering: 1. filter for certain threshold of hemoglobin genes 
# 2. Filter out all cells expressing hemoglobin genes
output_path <- paste0(wd, "data/output/")
if (file.exists(output_path)) {cat("The folder already exists")} else {dir.create(output_path)}

saveRDS(baseline_seurat_filtered, paste0(output_path, "/baseline_filtered_md_annot.rds"))
saveRDS(baseline_seurat_filtered_remove_rbc, paste0(output_path, "/baseline_filtered_no_rbc_md_annot.rds"))



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







