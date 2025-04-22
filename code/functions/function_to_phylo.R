
#Function to obtain the newick tree format from the ancestry tables

ToPhylo2 <- function(data){
  data.2 <- data
  data.2$repr <- data$spp
  sisters <- levels(as.factor(data$spp))
  mothers <- levels(as.factor(data$ancestor))
  tips <- setdiff(sisters, mothers)
  root <- setdiff(mothers, sisters)
  
  # the root might be ancestor=0
  if(length(root) == 0) root <- 0
  foc.nodes <- unique(data[which(data$spp %in% tips), "ancestor"])
  n <- length(foc.nodes)
  data.2$repr[data.2$spp %in% tips] <- data.2$repr[data.2$spp %in% tips]
  
  while(n > 1){
    foc.nodes2 <- unique(data.2[which(data.2$spp %in% foc.nodes), "ancestor"])
    for(i in 1:n){
      daughters <- data.2[which(data.2$ancestor == foc.nodes[i]), "repr"]
      daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[i]), "distance"]
      
      # This block handles the case where the ancestor is still extant (coexists with mutants)
      if(foc.nodes[i] %in% data.2$spp){
        daughters <- c(daughters, foc.nodes[i])
        daughters.dist <- c(daughters.dist, 0)  # Set ancestor's distance as 0
      }
      
      data.2$repr[data.2$spp == foc.nodes[i]] <- paste0(sister.group(daughters, daughters.dist), foc.nodes[i])
    }
    tips <- foc.nodes
    foc.nodes <- foc.nodes2
    n <- length(foc.nodes)
  }
  
  daughters <- data.2[which(data.2$ancestor == foc.nodes[1]), "repr"]
  daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[1]), "distance"]
  
  paste0(sister.group(daughters, daughters.dist), root, ";")  # Ensuring the tree ends with a semicolon
}