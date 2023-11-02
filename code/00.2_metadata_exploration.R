#Author: Dzhansu Hasanova
#Metadata exploration for baseline samples

library(ggplot2)
library(gridExtra)

set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Get baseline samples
samples_counts <- read_excel(paste0(wd, "data/metadata/samples_counts.xlsx"))
#Create df for baseline
baseline <- samples_counts[samples_counts$timepoint == 'Baseline',]$patient

#import metadata
metadata <- read_excel(paste0(wd, "data/metadata/metadata.xlsx"))
colnames(metadata) <- gsub(" ", "_", colnames(metadata))

metadata_baseline <- metadata[metadata$Sample_ID %in% baseline,]

#Specify path to save figures
path_qc_fig <- paste0(wd, "figures/QC/metadata")
if (file.exists(path_qc_fig)) {cat("The folder already exists")} else {dir.create(path_qc_fig)}

# Create plots with different metadata columns

ggplot(metadata_baseline, aes(x=age_at_diagnosis)) + 
  geom_density(fill = "deepskyblue4", alpha = 0.6)+
  theme_bw() + xlab("Age at diagnosis")
ggsave("age.png", path = path_qc_fig, width = 5, height = 6)

ggplot(metadata_baseline, aes(x=Sex)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()
ggsave("sex.png", path = path_qc_fig, width = 3, height = 4)

ggplot(metadata_baseline, aes(x=reorder(Histology,Histology,
                                          function(x)-length(x)))) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Histology")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("histology.png", path = path_qc_fig, width = 6, height = 6)


ggplot(metadata_baseline, aes(x=reorder(Molecular_subtype,Molecular_subtype,
                                        function(x)-length(x)))) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Molecular Subtype")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("molecular_subtype.png", path = path_qc_fig, width = 5, height =6)

ggplot(metadata_baseline, aes(x=reorder(Briding_therapy,Briding_therapy,
                                        function(x)-length(x)))) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Bridging Therapy")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("birdging_therapy.png", path = path_qc_fig, width = 4, height = 6)


ggplot(metadata_baseline, aes(x=reorder(Number_of_prior_therapies,Number_of_prior_therapies,
                                          function(x)-length(x)))) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Number of prior therpaies")
ggsave("prior_therapies.png", path = path_qc_fig, width = 4, height = 6)

ggplot(metadata_baseline, aes(x=Prior_autologous_transplant)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Prior autologous transplant")
ggsave("autologous_transplant.png", path = path_qc_fig, width = 3, height = 4)

ggplot(metadata_baseline, aes(x=Ann_Arbor_Stage)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+xlab("Ann Arbor Index")
ggsave("ann_arbor.png", path = path_qc_fig, width = 5, height = 6)


