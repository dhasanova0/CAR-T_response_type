library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(pals)
set.seed(100)

wd = "/Users/dhasanova/Documents/ETH/HS23/"

input_path_post <- paste0(wd, "data/output_posttreatment/")
seurat_obj_post <- readRDS(paste0(input_path_post, "post_doublet_filtered_md_annot.rds"))
colnames(seurat_obj_post@meta.data)[26] <- "label"

md_post <- as.data.frame(seurat_obj_post@meta.data)
md_post$label_id <- paste0(md_post$orig.ident, "_", md_post$label) 

path_fig_post <- paste0(wd, "figures/cell_type_distirbutions_post/")
if (file.exists(path_fig_post)) {cat("The folder already exists")} else {dir.create(path_fig_post)}

cols1 = c("#1F78C8", "#ff0000", "#33a02c", "#6A33C2","#36648B", "#ff7f00", "#FFD700", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#00E2E5", "#00FF00", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C")

cols2 = c("#1F78C8", "#33a02c","#ff0000", "#00E2E5", "#6A33C2","#36648B", "#ff7f00", "#FFD700", "#a6cee3",
          "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#999999", "#EEE685", "#565656",
          "#FF83FA", "#C814FA", "#0000FF", "#C8308C", "#00FF00", "#778B00", "#BEBE00", "#8B3B00", "#A52A3C")

cell_type <- as.data.frame(table(seurat_obj_post@meta.data$cell_type))

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
  
  
  
  #sorted_columns <- c(rep("Patient10-D7_NR", 10), rep("Patient10-D7_R", 10), rep("Patient11-D7_NR",10),
                      #rep("Patient11-D7_R",10), rep("Patient12-D14_NR",10), rep("Patient12-D14_R",10),
                      #rep("Patient12-D7_NR",10), rep("Patient12-D7_R",10), rep("Patient13-D7_NR",10),
                      #rep("Patient13-D7_R",10), rep("Patient12-Baseline_R",10), rep("Patient25-Baseline_NR",10),
                      #rep("Patient21-Baseline_R",10), rep("Patient30-Baseline_R",10), rep("Patient31-Baseline_NR",10),
                      #rep("Patient15-Baseline_R",10), rep("Patient23-Baseline_NR",10), rep("Patient22-Baseline_R",10),
                      #rep("Patient10-Baseline_NR",10), rep("Patient11-Baseline_R",10))
  
  #order_indices <- order(match(tbl$label_id, sorted_columns))
  
  
  
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

tbl <- create_table_celltye_id(md_post)

ggplot(data=tbl, aes(x=factor(reorder(cell_type, -Freq)), Freq, fill = factor(label_id, levels = unique(tbl$label_id))))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=unname(cols25()), name = NULL)+
  ylab("Number of cells")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )
ggsave("complete_cell-type_id_dist_corr.png", path = path_fig, height = 6, width = 9)

ggplot(data=tbl, aes(x=factor(label_id, levels = unique(label_id)), percentage, fill = cell_type))+
  geom_col(position = 'stack')+
  scale_fill_manual(values=cols2, name = NULL)+
  ylab("%")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), )+
  coord_flip() + scale_x_discrete(limits = rev(unique(tbl$label_id)))
ggsave("complete_cell-type_id_dist_per_corr.png", path = path_fig_post, height = 6, width = 7)


