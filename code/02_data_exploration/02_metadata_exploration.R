#Author: Dzhansu Hasanova
#Metadata exploration for baseline samples

library(ggplot2)
library(gridExtra)
library(readxl)
library(tidyverse)
library(plyr)


set.seed(12342)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

#Get baseline samples
samples_counts <- read_excel(paste0(wd, "data/metadata/samples_counts.xlsx"))
#Create df for baseline
baseline <- samples_counts[samples_counts$timepoint == 'Baseline',]$patient
post <- samples_counts[samples_counts$timepoint == 'D7' | samples_counts$timepoint == 'D14',]$patient
#import metadata
metadata <- read_excel(paste0(wd, "data/metadata/metadata.xlsx"))
colnames(metadata) <- gsub(" ", "_", colnames(metadata))

#Correct values to make them consistnet across rows
metadata$Response_classification[metadata$Response_classification == 'PR'] <- 'NR'
metadata$Ann_Arbor_Stage[metadata$Ann_Arbor_Stage == 'iII'] <- 'III'
metadata$Molecular_subtype[metadata$Molecular_subtype == 'GCB-double hit'] <- 'GCB double hit'
metadata$Molecular_subtype[metadata$Molecular_subtype == 'double hit'] <- 'GCB double hit'
metadata$Molecular_subtype[metadata$Molecular_subtype == 'GCB, triple hit'] <- 'GCB triple hit'
metadata$Molecular_subtype[metadata$Molecular_subtype == 'nonGCB, non double expressor'] <- 'nonGCB'

metadata <- metadata %>%
  mutate(Briding_therapy = case_when(
    str_detect(Briding_therapy, " ") ~ "yes",
    TRUE ~ Briding_therapy
  )
  )

colnames(metadata)[colnames(metadata) == "Baseline_tumor_burden_(SPD,_cm2)"] <- "baseline_tumor_burden"

#Convert to correct type
metadata$age_at_diagnosis <- as.numeric(metadata$age_at_diagnosis)
metadata$age_at_diagnosis <- round_any(metadata$age_at_diagnosis, 10)
metadata$age_at_diagnosis <- as.character(metadata$age_at_diagnosis)

metadata$baseline_tumor_burden <- as.numeric(metadata$baseline_tumor_burden)

#Calculate p-values for association of given clinical, patient feature with metadata

correlation <- data.frame(matrix(ncol = 2, nrow = 11))
rownames(correlation) <- c("Sex", "Briding_therapy",
                           "Prior_autologous_transplant",
                           "Product", "Ann_Arbor_Stage", "Number_of_prior_therapies", "age_at_diagnosis", "baseline_tumor_burden", "Molecular_subtype", "ECOG", 
                           "Prior_imid")
colnames(correlation) <- c("p_value", "adj_p_value")

#Perform Fisher's exact test for categorical values
for(i in rownames(correlation)){
  print(i)
  p_val <- fisher.test(table(metadata[, c("Response_classification", i)]))$p.value
  correlation[i,]$p_value <- p_val
}

#Perform Wilcoxon test for numerical values
correlation["baseline_tumor_burden",]$p_value <- wilcox.test(baseline_tumor_burden ~ Response_classification, data = metadata)$p.value

#Adjust p-values with BH
correlation$adj_p_value <- p.adjust(correlation$p_value,method = "BH")

write.csv(correlation, paste0(wd, "data/p_values_md.csv"))


#Specify path to save figures
path_qc_baseline <- paste0(wd, "figures/baseline/QC/metadata")
if (file.exists(path_qc_baseline)) {cat("The folder already exists")} else {dir.create(path_qc_baseline)}

path_qc_post <- paste0(wd, "figures/posttreatment/QC/metadata")
if (file.exists(path_qc_post)) {cat("The folder already exists")} else {dir.create(path_qc_post)}

# Create function for plots with different metadata columns

md_fig <- function(md, path_fig){
  ggplot(md, aes(x=age_at_diagnosis)) + 
    geom_density(fill = "deepskyblue4", alpha = 0.6)+
    theme_bw() + xlab("Age at diagnosis")
  ggsave("age.png", path = path_fig, width = 5, height = 6)
  
  ggplot(md, aes(x=Sex)) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()
  ggsave("sex.png", path = path_fig, width = 3, height = 4)
  
  ggplot(md, aes(x=reorder(Histology,Histology,
                                          function(x)-length(x)))) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Histology")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  ggsave("histology.png", path = path_fig, width = 6, height = 6)
  
  
  ggplot(md, aes(x=reorder(Molecular_subtype,Molecular_subtype,
                                          function(x)-length(x)))) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Molecular Subtype")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  ggsave("molecular_subtype.png", path = path_fig, width = 5, height =6)
  
  ggplot(md, aes(x=reorder(Briding_therapy,Briding_therapy,
                                          function(x)-length(x)))) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Bridging Therapy")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
  ggsave("birdging_therapy.png", path = path_fig, width = 4, height = 6)
  
  
  ggplot(md, aes(x=reorder(Number_of_prior_therapies,Number_of_prior_therapies,
                                          function(x)-length(x)))) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Number of prior therpaies")
  ggsave("prior_therapies.png", path = path_fig, width = 4, height = 6)
  
  ggplot(md, aes(x=Prior_autologous_transplant)) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Prior autologous transplant")
  ggsave("autologous_transplant.png", path = path_fig, width = 3, height = 4)
  
  ggplot(md, aes(x=Ann_Arbor_Stage)) + 
    geom_bar(stat = "count", fill = "deepskyblue4")+
    theme_bw()+xlab("Ann Arbor Index")
  ggsave("ann_arbor.png", path = path_fig, width = 5, height = 6)
}


#Get baseline and post-treatment metadata
metadata_baseline <- metadata[metadata$Sample_ID %in% baseline,]
metadata_post <- metadata[metadata$Sample_ID %in% post,]


md_fig(metadata_baseline, path_qc_baseline)
md_fig(metadata_post, path_qc_post)

