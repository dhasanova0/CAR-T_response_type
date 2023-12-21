library(Seurat)
library("stringr")
library(ggplot2)
library(SeuratDisk)
library(dplyr)


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