#Author: Dzhansu Hasanova 10042023
#Sequencing data cell count, umi etc.

#Import libraries

library(ggplot2)
library(gridExtra)
library(readxl)


#Import data
samples_counts <- read_excel("Documents/ETH/HS23/data/samples_counts.xlsx")

colnames(samples_counts) <- gsub(" ", "_", colnames(samples_counts))
colnames(samples_counts) <- gsub("/", "_per_", colnames(samples_counts))

#Function to plot various distributions
baseline <- samples_counts[samples_counts$timepoint == 'Baseline',]
post <- samples_counts[samples_counts$timepoint == 'D7' | samples_counts$timepoint == 'D14',]
infusion <- samples_counts[samples_counts$timepoint == 'Infusion',]


path_qc_fig <- "/Users/dhasanova/Documents/ETH/HS23/figures/QC"

grid_plot <- function(df, title, path_plot, w, h){
  cell_count <- ggplot(data=df, aes(x=patient, y=Number_of_cells)) +
    geom_bar(stat="identity", fill = "deepskyblue4")+
    ylab("Number of cells")+xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  cat("sum", sum(df$Number_of_cells), "\n")
  
  umi <- ggplot(data=df, aes(x=patient, y=Median_UMIs_per_cell)) +
    geom_bar(stat="identity", fill = "deepskyblue4")+
    ylab("Median UMIs/cell")+xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  
  genes <- ggplot(data=df, aes(x=patient, y=Median_genes_detected_per_cell)) +
    geom_bar(stat="identity", fill = "deepskyblue4")+
    ylab("Median genes detected/cell")+xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) 
  
  #jpeg(paste0(path_plot, "/", title, ".png"))
  grid.arrange(cell_count,umi,genes, top = title, ncol = 3, nrow  = 1)
  #dev.off()
  
  
 
}

grid_plot(baseline, "Baseline samples", path_qc_fig, 20, 10)
grid_plot(post, "Post-treatment samples", path_qc_fig, 30,10)
grid_plot(infusion, "Infusion samples", path_qc_fig, 40, 10)


