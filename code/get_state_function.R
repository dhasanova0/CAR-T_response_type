get_state <- function(seu_obj, states){
  #Specify in what state the cells are
  for (i in unique(states$Cluster)){
    cells <- c()
    for (m in which(states$Cluster == i)){
      state <- rownames(seu_obj@meta.data[seu_obj@meta.data[, paste0("MFI_", m)] == "yes",])
      cells <- c(cells,state)
    }
    
    cells <- unique(cells)
    seu_obj@meta.data[ , ncol(seu_obj@meta.data) + 1] <- ifelse(rownames(seu_obj@meta.data) %in% cells, "yes", "no")
    colnames(seu_obj@meta.data)[ncol(seu_obj@meta.data)] <- i
  }
  
  return(seu_obj)
}