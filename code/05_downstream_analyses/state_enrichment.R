#Author: Dzhansu Hasanova
#This script calculates state enrichment for given metadata entry

library(Seurat)
library("stringr")
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(pheatmap)
source("~/CAR-T_response_type/code/get_state_function.R")

#Define working directory
wd = "/Users/dhasanova/Documents/ETH/HS23/"

get_enrichment <- function(seu_obj, md, states, path_pdf, path_fold, order, cluster_row, cluster_col){
  #Calculate enrichment of cells in states with respect ot md value
  #Input: seu_obj: seurat object
  #       md: metadata entry as string for which enrichment is calculated
  #       states: States table as output from Stator
  #       path_df, path_fold: Path for saving enrichment plot and csv file with folds
  
  
  md_l <- unique(seu_obj@meta.data[, md])
  
  cell_list <- vector(mode = "list", length = length(md_l))
  cell_list_en <- vector(mode = "list", length = length(md_l))
  
  for (i in 1:length(md_l)){
    cell_list[[i]] <- as.data.frame(seu_obj@meta.data[seu_obj@meta.data[, md] %in% c(md_l[i]),])
    cell_list_en[[i]] <- data.frame(matrix(NA, nrow =length(unique(states$Cluster)), ncol = 2))
    colnames(cell_list_en[[i]]) <- c("Var1", "Freq")
    cell_list_en[[i]]$Var1 <- gsub("Cluster:", "", unique(states$Cluster))
    cell_list_en[[i]]$Freq <- NA
  }
  
  
  
  for (i in 1:length(unique(states$Cluster))){
    print(i)
    
    for (c in 1:length(cell_list)){
      cell_list_en[[c]]$Freq[i] <- nrow(cell_list[[c]][cell_list[[c]][, unique(states$Cluster)[i]] == "yes",])
      
      
      cell_list_en[[c]]$Var1[i] <- paste0(i, " (",nrow(seu_obj@meta.data[seu_obj@meta.data[, unique(states$Cluster)[i]] == "yes",]), ", ",
                                          round((nrow(seu_obj@meta.data[seu_obj@meta.data[, unique(states$Cluster)[i]] == "yes",]) / nrow(seu_obj@meta.data))*100, digits = 2), "%)")
      
      
    }
    
    
  }
  
  
  r=length(cell_list)
  col=length(cell_list_en[[1]]$Var1)
  Data_mtrix<-array(data=NA,dim = c(r,col))
  colnames(Data_mtrix) <- c(cell_list_en[[1]]$Var1)
  rownames(Data_mtrix) <- md_l
  
  Fold<-array(data=NA,dim = c(r,col))
  colnames(Fold) <- c(cell_list_en[[1]]$Var1)
  rownames(Fold) <- md_l
  
  P_set=rep(NA,length(cell_list_en[[1]]$Var1))
  names(P_set)<-cell_list_en[[1]]$Var1
  
  fold_set=rep(NA,length(cell_list_en[[1]]$Var1))
  names(fold_set)<-cell_list_en[[1]]$Var1
  
  
  #Perform enrichment analysis
  for (i in 1:length(cell_list_en)){
    
    for (cell in cell_list_en[[i]]$Var1){
      
      print(cell)
      q = as.numeric(cell_list_en[[i]][cell_list_en[[i]]$Var1 == cell,]$Freq) #Number of cell of interest in group
      m = as.numeric(nrow(seu_obj@meta.data[seu_obj@meta.data[, paste0("Cluster:", sub("\\ .*", "", cell))] == "yes",])) #Total number of cell of interest in entire data
      n = as.numeric(nrow(seu_obj@meta.data) - m) #Number of cell without cell of interest
      k = as.numeric(nrow(seu_obj@meta.data[seu_obj@meta.data[, md] == md_l[i],])) #Number of celstate_list in group
      p_value <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
      fold_value <- (q*(m+n))/(m*k)
      P_set[cell] <- p_value
      fold_set[cell] <- fold_value
      
      for (diff in setdiff(cell_list_en[[1]]$Var1, cell_list_en[[i]]$Var1)){
        P_set[diff] <- NA
        fold_set[diff] <- NA
      }
      
    }
    Data_mtrix[i,] <- P_set
    Fold[i,] <- fold_set
  }
  
  
  
  #Adjust p-values and get -log10
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<- -log10(Mydata_raw_m)
  
  
  colnames(log10FDR) <- cell_list_en[[1]]$Var1
  rownames(log10FDR) <- md_l
  
  
  row_name <- c()
  
  for (i in md_l){
    row_name <- c(row_name,paste0(i," (",nrow(seu_obj@meta.data[seu_obj@meta.data[, md] == i,]), ", ",
                                  round((nrow(seu_obj@meta.data[seu_obj@meta.data[, md] == i,])/nrow(seu_obj@meta.data))*100, digits = 2), "%)") )
  }
  
  rownames(log10FDR) <- row_name
  
  #Significance
  p_label <- Mydata_raw_m
  p_label[p_label <= 0.05] <- "*"
  p_label[p_label > 0.05] <- ""
  
  if (order == TRUE){
    log10FDR <- arrange(as.data.frame(log10FDR), rownames(log10FDR))
    heatmap_joint_n <- as.data.frame(sapply(log10FDR,as.numeric))
    p_label <- apply(heatmap_joint_n, MARGIN = 2, function(x) ifelse(x > 1.3, "*", ""))
    
    rownames(heatmap_joint_n) <- rownames(log10FDR)
    
    log10FDR <- data.matrix(heatmap_joint_n)
  }
  
  
  pdf(path_pdf, width=10, height=17)
  

  ht <- ComplexHeatmap::pheatmap(log10FDR, display_numbers = p_label, fontsize_number=10, cellwidth = 15,cellheight = 15, 
                                 cluster_rows = cluster_row,show_row_dend = FALSE, border_color = NA,fontsize = 12,
                                 cluster_cols = cluster_col, show_column_dend = FALSE, heatmap_legend_param=list(title_position = "lefttop-rot"),
                                 name = "-Log10FDR", color = colorRampPalette(c("white","firebrick3"))(100), 
                                 main = "State enrichment")
  
  draw(ht,heatmap_legend_side = "left")
  dev.off()
  
  write.csv(Fold, path_fold)
  
}


#Get state enrichment for "validation set" for states from run1
baseline_run1 <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_valid.rds"))
states_baseline_r1 <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-2023-11-17.csv"))

baseline_run1@meta.data$cell_label <- paste0(baseline_run1@meta.data$cell_type, "_", baseline_run1@meta.data$label)
baseline_run1@meta.data$subtype_label <- paste0(baseline_run1@meta.data$subtype, "_", baseline_run1@meta.data$label)

baseline_run1 <- get_state(baseline_run1, states_baseline_r1)



# ---------------- Baseline -----------------------

get_enrichment(baseline_run1, 'label', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_label_base.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_label_base_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_run1, 'cell_type', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_base.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_base_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_run1, 'cell_label', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_joint_base.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_joint_base_fold.csv"), TRUE, FALSE, FALSE)

# ---------------- Baseline 0.5 -----------------------
states_baseline_r1_0.5 <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-0.5.csv"))
baseline_r1_0.5 <- baseline_run1
baseline_r1_0.5@meta.data <- baseline_r1_0.5@meta.data %>% select(-contains('Cluster'))

states_baseline_r1_0.5 <- states_baseline_r1_0.5[match(states_baseline_r1$X, states_baseline_r1_0.5$X),]
baseline_r1_0.5 <- get_state(baseline_r1_0.5, states_baseline_r1_0.5)

#Order columns
# Get column names
col_names <- names(baseline_r1_0.5@meta.data)
# Extract numeric part from column names
numeric_part <- as.numeric(gsub("Cluster:", "", col_names))
# Identify the index of other columns
other_cols_index <- which(!grepl("Cluster", col_names))
ordered_clusters <- paste0("Cluster:", sort(as.numeric(numeric_part[!is.na(numeric_part)])))
cluster_cols_index <- match(ordered_clusters, col_names)

# Order column names based on numeric part
new_col_order <- c(other_cols_index, cluster_cols_index)

# Rearrange columns in the data frame
baseline_r1_0.5@meta.data <- baseline_r1_0.5@meta.data[, new_col_order]

states_baseline_r1_0.5 <- read.csv(paste0(wd, "data/stator_results/run1/State_Table-0.5.csv"))

get_enrichment(baseline_r1_0.5, 'cell_type', states_baseline_r1_0.5, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_base_0.5.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_fold_0.5.csv"), TRUE, FALSE, FALSE)


get_enrichment(baseline_r1_0.5, 'cell_label', states_baseline_r1_0.5, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_joint_base_0.5.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_cell_joint_base_fold_0.5.csv"), TRUE, FALSE, FALSE)


get_enrichment(baseline_r1_0.5, 'subtype_label', states_baseline_r1_0.5, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_subtype_joint_base_0.5.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_subtype_joint_base_fold_0.5.csv"), TRUE, FALSE, FALSE)

get_enrichment(baseline_r1_0.5, 'label', states_baseline_r1_0.5, 
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_label_base_0.5.pdf"),
               paste0(wd, "figures/stator_run1/thesis/heatmap_state_label_base_fold_0.5.csv"), TRUE, FALSE, FALSE)

# ---------------- Baseline Balanced ID -----------------------

#Get state enrichment for balanced id in baseline samples for states from run1 (validation set)
baseline_run1_balanced_id <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_balanced_id.rds"))

get_enrichment(baseline_run1_balanced_id, 'label', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_label_base_balanced_id.pdf"),
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_label_base_balanced_fold_id.csv"), FALSE, FALSE, FALSE )

get_enrichment(baseline_run1_balanced_id, 'cell_type', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_cell_base_balanced_id.pdf"),
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_cell_base_balanced_fold_id.csv"), FALSE, FALSE, FALSE )

baseline_run1_balanced_id@meta.data$id_label <- paste0(baseline_run1_balanced_id@meta.data$label, "_", baseline_run1_balanced_id@meta.data$orig.ident)
get_enrichment(baseline_run1_balanced_id, 'id_label', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_orig_base_balanced_id.pdf"),
               paste0(wd, "figures/stator_run1/thesis/balanced/ID_state_orig_base_balanced_fold_id.csv"), TRUE, FALSE, FALSE )

# ---------------- Baseline Balanced Label -----------------------

#Get state enrichment for balanced cells in baseline samples for states from run1 (validation set)
baseline_run1_balanced_cell <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_balanced.rds"))

get_enrichment(baseline_run1_balanced_cell, 'label', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/balanced/CELL_state_label_base_balanced_cell.pdf"),
               paste0(wd, "figures/stator_run1/thesis/balanced/CELL_state_label_base_balanced_fold_cell.csv"), FALSE, FALSE, FALSE )

get_enrichment(baseline_run1_balanced_cell, 'cell_type', states_baseline_r1, 
               paste0(wd, "figures/stator_run1/thesis/balanced/CELL_state_cell_base_balanced_cell.pdf"),
               paste0(wd, "figures/stator_run1/thesis/balanced/CELL_state_cell_base_balanced_fold_cell.csv"), FALSE, FALSE, FALSE )


# ---------------- Baseline Filtered  0.83 -----------------------

baseline_filtered <- readRDS(paste0(wd, "data/output_baseline/states/base_filtered_states_valid.rds"))
baseline_filtered@meta.data$label_id <- paste0(baseline_filtered@meta.data$label, "_", baseline_filtered@meta.data$orig.ident)
baseline_filtered@meta.data$cell_label <- paste0(baseline_filtered@meta.data$cell_type, "_", baseline_filtered@meta.data$label)

states_filtered_0.83 <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-filtered.csv"))
baseline_filtered_0.83 <- get_state(baseline_filtered, states_filtered_0.83)


get_enrichment(baseline_filtered_0.83, 'label', states_filtered_0.83, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_label_filtered.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_label_filtered_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_filtered_0.83, 'cell_type', states_filtered_0.83, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_filtered.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_filtered_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_filtered_0.83, 'label_id', states_filtered_0.83, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_id_filtered.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_id_filtered_fold.csv"), TRUE, FALSE, FALSE )

get_enrichment(baseline_filtered_0.83, 'cell_label', states_filtered_0.83, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_joint_filtered.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_joint_filtered_fold.csv"), TRUE, FALSE, FALSE )

# ---------------- Baseline Filtered  0.5 -----------------------

states_filtered_0.5 <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-2024-0.5.csv"))
states_filtered_0.5 <- states_filtered_0.5[match(states_filtered_0.83$X, states_filtered_0.5$X),]
baseline_filtered_0.5 <- get_state(baseline_filtered, states_filtered_0.5)
#Order columns
# Get column names
col_names <- names(baseline_filtered_0.5@meta.data)
# Extract numeric part from column names
numeric_part <- as.numeric(gsub("Cluster:", "", col_names))
# Identify the index of other columns
other_cols_index <- which(!grepl("Cluster", col_names))
ordered_clusters <- paste0("Cluster:", sort(as.numeric(numeric_part[!is.na(numeric_part)])))
cluster_cols_index <- match(ordered_clusters, col_names)

# Order column names based on numeric part
new_col_order <- c(other_cols_index, cluster_cols_index)

# Rearrange columns in the data frame
baseline_filtered_0.5@meta.data <- baseline_filtered_0.5@meta.data[, new_col_order]
states_filtered_0.5 <- read.csv(paste0(wd, "data/stator_results/run1_filtered/shiny/State_Table-2024-0.5.csv"))


get_enrichment(baseline_filtered_0.5, 'label', states_filtered_0.5, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_label_filtered_0.5.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_label_filtered_0.5_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_filtered_0.5, 'cell_type', states_filtered_0.5, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_filtered_0.5.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_filtered_0.5_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_filtered_0.5, 'label_id', states_filtered_0.5, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_id_filtered_0.5.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_id_filtered_0.5_fold.csv"), TRUE, FALSE, FALSE )

get_enrichment(baseline_filtered_0.5, 'cell_label', states_filtered_0.5, 
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_joint_filtered_0.5.pdf"),
               paste0(wd, "data/stator_results/run1_filtered/projected/heatmap_state_cell_joint_filtered_0.5_fold.csv"), TRUE, FALSE, FALSE )

# ---------------- Stator run 2 baseline -----------------------

baseline_run2 <- readRDS(paste0(wd, "data/output_baseline/states/baseline_state_run2_valid.rds"))
states_baseline_r2 <- read.csv(paste0(wd, "data/stator_results/run3/run3/results/State_Table-2024-02-01.csv"))

get_enrichment(baseline_run2, 'label', states_baseline_r2, 
               paste0(wd, "data/stator_results/run3/run3/results/valid/heatmap_state_label.pdf"),
               paste0(wd, "data/stator_results/run3/run3/results/valid/heatmap_state_label_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(baseline_run2, 'cell_type', states_baseline_r2, 
               paste0(wd, "data/stator_results/run3/run3/results/valid/heatmap_state_cell.pdf"),
               paste0(wd, "data/stator_results/run3/run3/results/valid/heatmap_state_cell_fold.csv"), FALSE, TRUE, TRUE )



# ---------------- Post Projected States from baseline -----------------------

#Post-treatment with baseline states from run1

post_run1_baseline <- readRDS(paste0(wd, "data/output_baseline/states/post_state_baseline.rds"))
post_run1_baseline <- get_state(post_run1_baseline, states_baseline_r1)

get_enrichment(post_run1_baseline, 'Response_classification', states_baseline_r1, 
               paste0(wd, "figures/posttreatment/states/heatmap_state_label_post.pdf"),
               paste0(wd, "figures/posttreatment/states/heatmap_state_label_post_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(post_run1_baseline, 'cell_type', states_baseline_r1, 
               paste0(wd, "figures/posttreatment/states/heatmap_state_cell_post.pdf"),
               paste0(wd, "figures/posttreatment/states/heatmap_state_cell_post_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(post_run1_baseline, 'Product', states_baseline_r1, 
               paste0(wd, "figures/posttreatment/states/heatmap_state_product_post.pdf"),
               paste0(wd, "figures/posttreatment/states/heatmap_state_product_post_fold.csv"), FALSE, TRUE, TRUE )


# ---------------- Post axi-cel D7 -----------------------

#Post-treatment axi-cel D7


axi_D7 <- readRDS(paste0(wd, "data/output_posttreatment/states/seu_obj_axi_D7_valid.rds"))
states_axi_D7 <- read.csv(paste0(wd, "data/stator_results/axi_D7_1/Shiny_results/State_Table-2024-01-03.csv"))

get_enrichment(axi_D7, 'Response_classification', states_axi_D7, 
               paste0(wd, "data/output_posttreatment/states/axi_D7/heatmap_state_label_post.pdf"),
               paste0(wd, "data/output_posttreatment/states/axi_D7/heatmap_state_label_post_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(axi_D7, 'cell_type', states_axi_D7, 
               paste0(wd, "data/output_posttreatment/states/axi_D7/heatmap_state_cell_post.pdf"),
               paste0(wd, "data/output_posttreatment/states/axi_D7/heatmap_state_cell_post_fold.csv"), FALSE, TRUE, TRUE )

# ---------------- Post tisa-cel D7 -----------------------
#Post-treatment tisa-cel D7

tisa_D7 <- readRDS(paste0(wd, "data/output_posttreatment/states/seu_obj_tisa_D7_umap.rds"))
states_tisa_D7 <- read.csv(paste0(wd, "data/stator_results/tisa_D7/Shiny_results/State_Table-2024-02-08.csv"))

get_enrichment(tisa_D7, 'Response_classification', states_tisa_D7, 
               paste0(wd, "data/output_posttreatment/states/tisa_D7/heatmap_state_label_post.pdf"),
               paste0(wd, "data/output_posttreatment/states/tisa_D7/heatmap_state_label_post_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(tisa_D7, 'cell_type', states_tisa_D7, 
               paste0(wd, "data/output_posttreatment/states/tisa_D7/heatmap_state_cell_post.pdf"),
               paste0(wd, "data/output_posttreatment/states/tisa_D7/heatmap_state_cell_post_fold.csv"), FALSE, TRUE, TRUE )


# ---------------- Baseline T cell (not part of the thesis) -----------------------
#Baseline T cell

T_cell <- readRDS(paste0(wd, "data/output_baseline/states/base_T_cell.rds"))
states_T <- read.csv(paste0(wd, "data/stator_results/T_cell/shiny/State_Table-2024-02-22.csv"))
T_cell <- get_state(T_cell, states_T)

get_enrichment(T_cell, 'label', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_label_all.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_label_fold_all.csv"), FALSE, TRUE, TRUE )

get_enrichment(T_cell, 'cell_type', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_cell.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_cell_fold.csv"), FALSE, TRUE, TRUE )

get_enrichment(T_cell, 'subtype', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_subtype.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_subtype_fold.csv"), FALSE, TRUE, TRUE )

T_cell@meta.data$label_cell <- paste0(T_cell@meta.data$cell_type, "_" ,T_cell@meta.data$label)
T_cell@meta.data$label_id <- paste0(T_cell@meta.data$label, "_" ,T_cell@meta.data$orig.ident)
T_cell@meta.data$label_subtype <- paste0(T_cell@meta.data$subtype, "_" ,T_cell@meta.data$label)

get_enrichment(T_cell, 'label_cell', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_cell_joint.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_cell_joint_fold.csv"), TRUE, FALSE, FALSE )

get_enrichment(T_cell, 'label_subtype', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_subtype_joint.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_subtype_joint_fold.csv"), TRUE, FALSE, FALSE )

get_enrichment(T_cell, 'label_id', states_T, 
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_id.pdf"),
               paste0(wd, "data/stator_results/T_cell/projected/heatmap_state_id_fold.csv"), TRUE, FALSE, FALSE )

