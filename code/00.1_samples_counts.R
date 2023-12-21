#Author: Dzhansu Hasanova
#Data exploration: distributon of cell, umi, genes per patient

#Import libraries

library(ggplot2)
library(gridExtra)
library(readxl)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Import data
samples_counts <- read_excel(paste0(wd, "data/metadata/samples_counts.xlsx"))
#Change column names
colnames(samples_counts) <- gsub(" ", "_", colnames(samples_counts))
colnames(samples_counts) <- gsub("/", "_per_", colnames(samples_counts))

#Create df for the different time points
baseline <- samples_counts[samples_counts$timepoint == 'Baseline',]
post <- samples_counts[samples_counts$timepoint == 'D7' | samples_counts$timepoint == 'D14',]
post$patient_time <- paste0(post$patient, "_", post$timepoint)
colnames(post)[1] <- "patient_id"
colnames(post)[9] <- "patient"
post <- post[-c(31),]
infusion <- samples_counts[samples_counts$timepoint == 'Infusion',]

#Specify path to save figures
path_qc_fig <- paste0(wd, "figures/QC")
if (file.exists(path_qc_fig)) {cat("The folder already exists")} else {dir.create(path_qc_fig)}


#Function to plot various distributions
grid_plot <- function(df, title, path_plot, w, h){
  #Input: 
  #       df: dataframe with metadata
  #       title: title for plots (timepoint)
  #       path_plot: directory where plot is saved
  #       w,h: width and height of plot
  #
  #Output: plot saved in the specified directory
  
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
  
  pdf(paste0(path_plot, "/", title, ".pdf"), width = w, height = h) # Open a new pdf file
  grid.arrange(cell_count,umi,genes, top = title, ncol = 3, nrow  = 1) 
  dev.off()
}

grid_plot(baseline, "Baseline", path_qc_fig, 15, 7)
grid_plot(post, "Post-treatment", path_qc_fig, 17,7)
grid_plot(infusion, "Infusion", path_qc_fig, 25, 7)


