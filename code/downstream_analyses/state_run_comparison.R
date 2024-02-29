#Author: Dzhansu Hasanova
#Compare states between the different Stator runs

set.seed(124)

library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)
library(ggsignif)
library(ggpubr)

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"


baseline_states_run1 <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))
baseline_states_run1$Genes <- strsplit(baseline_states_run1$Genes, "_")
baseline_states_run1$D.tuple <- gsub("\\D", "", baseline_states_run1$D.tuple)

baseline_states_run2 <- read.csv(paste0(wd, "data/stator_results/run3/run3/results/State_Table-2024-02-01.csv"))
baseline_states_run2$Genes <- strsplit(baseline_states_run2$Genes, "_")
baseline_states_run2$D.tuple <- gsub("\\D", "", baseline_states_run2$D.tuple)

baseline_states_run1_filtered <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-2024-0.5.csv"))
baseline_states_run1_filtered$Genes <- strsplit(baseline_states_run1_filtered$Genes, "_")
baseline_states_run1_filtered$D.tuple <- gsub("\\D", "", baseline_states_run1_filtered$D.tuple)


post_axi_D7 <- read.csv(paste0(wd, "data/stator_results/axi_D7_1/Shiny results/State_Table-2024-01-03.csv"))
post_axi_D7$Genes <- strsplit(post_axi_D7$Genes, "_")
post_axi_D7$D.tuple <- gsub("\\D", "", post_axi_D7$D.tuple)

post_tisa_D7 <- read.csv(paste0(wd, "data/stator_results/tisa_D7/Shiny results/State_Table-2024-02-08.csv"))
post_tisa_D7$Genes <- strsplit(post_tisa_D7$Genes, "_")
post_tisa_D7$D.tuple <- gsub("\\D", "", post_tisa_D7$D.tuple)

#post_axi_D14 <- read.csv(paste0(wd, "data/stator_results/axi_D14")
#post_axi_D14$Genes <- strsplit(post_axi_D14$Genes, "_")
#post_axi_D14$D.tuple <- gsub("\\D", "", post_axi_D14$D.tuple)



get_shared_genes <- function(states_run1, states_run2, name_run1, name_run2){
  #Get genes which appear in states if the two Stator runs
  #Input: states_run1: Table of MFI and States of first run
  #       states_run2: Table of MFI and States of second run
  #       name_run1:   Name of first run
  #       name_run2:   Name of second run
  
  states_run1$Genes[states_run1$Genes %in% states_run2$Genes]
  
  genes_run1 <- c()
  genes_run2 <- c()
  
  for (i in 1:nrow(states_run1)){
    genes_run1 <- c(genes_run1, states_run1$Genes[[i]])
    genes_run2 <- c(genes_run2, states_run2$Genes[[i]])
  }
  
  shared <- intersect(genes_run1, genes_run2)
  
  empty_list <- vector(mode='list', length=length(shared))
  names(empty_list) <- shared
  
  for (g in shared){
    empty_list[[g]] <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(empty_list[[g]]) <- c("Gene", "d_tuple", "stator_run")
    for (i in 1:nrow(states_run1)){
      
      
      
      if (g %in% states_run1$Genes[[i]] == TRUE) {
        
        empty_list[[g]][nrow(empty_list[[g]]) + 1,]$Gene <- paste(states_run1$Genes[[i]], collapse = "_")
        empty_list[[g]][nrow(empty_list[[g]]),]$d_tuple <- paste(states_run1$D.tuple[[i]], collapse = "_")
        empty_list[[g]][nrow(empty_list[[g]]),]$stator_run <- name_run1
      }
      
      if (g %in% states_run2$Genes[[i]] == TRUE){
        empty_list[[g]][nrow(empty_list[[g]]) + 1,]$Gene <- paste(states_run2$Genes[[i]], collapse = "_")
        empty_list[[g]][nrow(empty_list[[g]]),]$d_tuple <- paste(states_run2$D.tuple[[i]], collapse = "_")
        empty_list[[g]][nrow(empty_list[[g]]),]$stator_run <- name_run2
        
      }
      
    }
  }
  
  return(empty_list)
  
}

run1_run2 <- get_shared_genes(baseline_states_run1, baseline_states_run2, "run1", "run2")
run1_axiD7 <- get_shared_genes(baseline_states_run1, post_axi_D7, "run1", "axi-D7")
run1_tisaD7 <- get_shared_genes(baseline_states_run1, post_tisa_D7, "run1", "tisa-D7")
axiD7_tisaD7 <- get_shared_genes(post_tisa_D7, post_axi_D7, "axi-D7", "tisa-D7")
run1_filtered <- get_shared_genes(baseline_states_run1, baseline_states_run1_filtered, "run1", "run1_filtered")

run1_filtered$CST3[duplicated(run1_filtered$CST3$Gene),]
run1_filtered$CAVIN2[duplicated(run1_filtered$CAVIN2$Gene),]


