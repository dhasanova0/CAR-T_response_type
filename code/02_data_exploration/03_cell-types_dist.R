#Author: Dzhansu Hasanova
#Create barplots of cell type distributions respective of sample

library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(pals)

set.seed(100)
wd = "/Users/dhasanova/Documents/ETH/HS23/"

input_path_complete <- paste0(wd, "data/output_baseline/")
seurat_obj_complete <- readRDS(paste0(input_path_complete, "baseline_doublet_filtered_md_annot.rds"))

input_path_subset1 <- paste0(wd, "data/output/stator_input/rds_new/")
seurat_obj_subset1 <- readRDS(paste0(input_path_subset1, "subset_1.rds"))

md_complete <- as.data.frame(seurat_obj_complete@meta.data)
md_subset1 <- as.data.frame(seurat_obj_subset1@meta.data)

md_celltype <- md_complete[md_complete$barcodes %in% md_subset1$barcodes,]
md_celltype <- md_celltype[c("orig.ident", "cell_type", "barcodes")]

md_subset1 <- subset(md_subset1, select = -c(cell_type))

md_merged <- left_join(md_subset1, md_celltype, by = join_by(barcodes == barcodes))
md_merged <- md_merged %>% mutate(cell_type = ifelse(is.na(cell_type), "not_defined", cell_type))
md_merged$cell_type[md_merged$cell_type == 'Unknown'] <- 'T cell'

md_merged$label_id <-  paste0(md_merged$orig.ident.x, "_", md_merged$label)
md_merged$cell_id <- paste0(md_merged$orig.ident.x, "_", md_merged$cell_type)


#Specify path to save figures
path_fig <- paste0(wd, "figures/cell_type_distirbutions/")
if (file.exists(path_fig)) {cat("The folder already exists")} else {dir.create(path_fig)}



cols1 = c("#1F78C8", "#ff0000", "#33a02c", "#ff7f00", "#6A33C2","#36648B", "#FFD700", "#a6cee3","#565656",
         "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685",
         "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#00E2E5", "#00FF00", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C")

cols2 = c("#1F78C8", "#33a02c", "#ff0000", "#00E2E5","#ff7f00",  "#6A33C2","#36648B", "#FFD700", "#a6cee3",
          "#999999", "#565656", "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#EEE685",
          "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#00FF00", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C")


get_md <- function(seurat_obj){
  md <- as.data.frame(seurat_obj@meta.data)
  md$label_id <- paste0(md$orig.ident, "_", md$label)
  md$cell_id <- paste0(md$orig.ident, "_", md$cell_type)
  
  return(md)
}
md_complete <- get_md(seurat_obj_complete)
md_subset1 <- get_md(seurat_obj_subset1)

md_post <- get_md(seurat_obj_post)

create_table_celltye_id <- function(md){
  md_id_cell <- md[, c("label_id", "cell_type")]
  
  
  md_id_cell_count <- md_id_cell %>%
    # count number of rows for each combination of server_id and protocol
    group_by(label_id, cell_type) %>%
    tally() %>%
    # pivot the protocol names over the columns
    pivot_wider(names_from=cell_type, values_from=n) %>%
    # replace NA values in all columns with 0
    mutate(across(everything(), .fns=~replace_na(., 0)))
  
  md_id_cell_count <- t(md_id_cell_count)
  colnames(md_id_cell_count) <- md_id_cell_count[1,]
  md_id_cell_count <- md_id_cell_count[-1,]
  
  md_id_cell_count <- as.data.frame(md_id_cell_count)
  
  tbl <- with(md_id_cell, table(cell_type, label_id))
  tbl <- as.data.frame(tbl)
  
  
  
  sorted_columns <- c(rep("Patient6-Baseline_NR", 10), rep("Patient13-Baseline_R", 10), rep("Patient17-Baseline_R",10),
                      rep("Patient24-Baseline_NR",10), rep("Patient18-Baseline_NR",10), rep("Patient19-Baseline_R",10),
                      rep("Patient9-Baseline_NR",10), rep("Patient8-Baseline_R",10), rep("Patient14-Baseline_NR",10),
                      rep("Patient20-Baseline_NR",10), rep("Patient12-Baseline_R",10), rep("Patient25-Baseline_NR",10),
                      rep("Patient21-Baseline_R",10), rep("Patient30-Baseline_R",10), rep("Patient31-Baseline_NR",10),
                      rep("Patient15-Baseline_R",10), rep("Patient23-Baseline_NR",10), rep("Patient22-Baseline_R",10),
                      rep("Patient10-Baseline_NR",10), rep("Patient11-Baseline_R",10))
  # Extracting indices for sorting
  sorted_indices <- order(gsub(".*-Baseline_(NR|R)", "\\1", sorted_columns), as.numeric(gsub("Patient(\\d+)-.*", "\\1", sorted_columns)))
  
  # Sorting the list
  sorted_columns <- sorted_columns[sorted_indices]
  
  order_indices <- order(match(tbl$label_id, sorted_columns))
  tbl <- tbl[order_indices, ]
  
  
  
  #cell_type <- unique(tbl$cell_type)
  #for (cell in cell_type){
  #  tbl[tbl$cell_type == cell,]$percentage <- round(tbl[tbl$cell_type == cell,]$Freq/(sum(tbl[tbl$cell_type == cell,]$Freq))*100, digits = 2)
  
  #}
  
  tbl$percentage <- NA
  
  ID <- unique(tbl$label_id)
  for (id in ID){
    tbl[tbl$label_id == id,]$percentage <- round(tbl[tbl$label_id == id,]$Freq/(sum(tbl[tbl$label_id == id,]$Freq))*100, digits = 2)
    
  }
  return(tbl)
  
}

tbl_complete <- create_table_celltye_id(md_complete)
tbl_susbet1 <- create_table_celltye_id(md_merged)

ggplot(data=tbl_complete, aes(x=factor(reorder(cell_type, -Freq)), Freq, fill = factor(label_id, levels = unique(tbl_complete$label_id))))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=unname(cols25()), name = NULL)+
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )
ggsave("complete_cell-type_id_dist_corr.png", path = path_fig, height = 6, width = 9)

ggplot(data=tbl_complete, aes(x=factor(label_id, levels = unique(tbl_complete$label_id)), percentage, fill = cell_type))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=cols2, name = NULL)+
  ylab("%")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), 
        axis.text=element_text(size=12), legend.text=element_text(size=12))+
  coord_flip() + scale_x_discrete(limits = rev(unique(tbl_complete$label_id)))
ggsave("complete_cell-type_id_dist_per.png", path = path_fig, height = 7, width = 9)

ggplot(data=tbl_susbet1, aes(x=factor(reorder(cell_type, -Freq)), Freq, fill = factor(label_id, levels = unique(tbl_susbet1$label_id))))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=unname(cols25()), name = NULL)+
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )
ggsave("subset1_cell-type_id_dist_corr.png", path = path_fig, height = 6, width = 9)

ggplot(data=tbl_susbet1, aes(x=factor(label_id, levels = unique(tbl_susbet1$label_id)), percentage, fill = cell_type))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=cols1, name = NULL)+
  ylab("%")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )+
  coord_flip() + scale_x_discrete(limits = rev(unique(tbl_susbet1$label_id)))
ggsave("subset1_cell-type_id_dist_per2_corr.png", path = path_fig, height = 6, width = 7)


# --------------- Simple distribution bar plots --------------------------------


md_R <- md[md$label == 'R',]
md_NR <- md[md$label == 'NR',]

counts_R <- data.frame(table(md_R$cell_type))
counts_NR <- data.frame(table(md_NR$cell_type))

counts_R$percentage <- round(counts_R$Freq/(sum(counts_R$Freq))*100, digits = 2)
counts_NR$percentage <- round(counts_NR$Freq/(sum(counts_NR$Freq))*100, digits = 2)

ggplot(data=counts_R, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  ggtitle("Responders")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("cell-dist_R_number.png", path = path_fig, height = 5, width = 6)

ggplot(data=counts_NR, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  ggtitle("Non-Responders")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("cell-dist_NR_number.png", path = path_fig, height = 5, width = 6)


ggplot(data=counts_NR, aes(x=reorder(Var1, -percentage), y=percentage)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  geom_text(aes(label = paste0(percentage, "%")), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("%")+xlab("")+
  theme_bw()+
  ggtitle("Non-Responders")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("cell-dist_NR_perc.png", path = path_fig, height = 5, width = 6)

ggplot(data=counts_R, aes(x=reorder(Var1, -percentage), y=percentage)) +
  geom_bar(stat="identity", fill = "deepskyblue4")+
  geom_text(aes(label = paste0(percentage, "%")), position = position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette = "Set2") +
  ylab("%")+xlab("")+
  theme_bw()+
  ggtitle("Responders")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
ggsave("cell-dist_R_perc.png", path = path_fig, height = 5, width = 6)

