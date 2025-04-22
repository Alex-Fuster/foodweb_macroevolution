




# Compute phylogenetic distances for each timestep from the anc matrix


compute_phylo_dist_all_steps <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  phylo_dist_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    phylo_dist_list[[i]] <- dist_mat_tree
    
    timestep_anc.matrix <- df_anc[i,1]
    
    names(phylo_dist_list)[[i]] <- paste("mat_timestep", timestep_anc.matrix)
    
    
  }
  
  
  return(phylo_dist_list)
  
}










# Compute phylogenetic trees for each timestep from the anc matrix


compute_phylog_all_timesteps <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  Phylo_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    Phylo_list[[i]] <- as.phylo(hclust(as.dist(dist_mat_tree), method = "average"))
    

    
  }
  
  
  return(Phylo_list)
  
}


# Same but without hclust

compute_phylog_all_timesteps_nohclust <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  Phylo_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    Phylo_list[[i]] <- as.dist(dist_mat_tree)
    
    Phylo_list[[i]]
    
    
  }
  
  
  return(Phylo_list)
  
}




# Plot trees of the trees list obtained in the simulation

plot_trees <- function(trees_list, plot_type) { # cladogram, phylogram, radial, fan
  
  list_plots <- list()
  
  for (i in 3:length(trees_list)) {
    
    list_plots[i] <- plot(unroot(trees_list[[i]]),type= plot_type,cex=0.6,
                         use.edge.length=FALSE,lab4ut="axial",
                         no.margin=TRUE)
    
    
  }
  

  
  return(list_plots)
  
}





# set column and row names as "spp_n", n being their position

set_sppNames_icolrows_spp <-  function(matrix) {
  
  mat <- matrix
  rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "spp_")
  colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "spp_")
  
  return(mat)
  
}




# set colnames and rownames as numbers


set_sppNames_numbers <- function(mat) {
  
  rownames(mat) <- seq(1:nrow(mat))
  colnames(mat) <- seq(1:ncol(mat))
  
  return(mat)  
}





#set_sppNames_icolrows <-  function(matrix) {
  
 # mat <- matrix
 # rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "")
 # colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "")
  
#  return(mat)
  
#}



# convert colnames and rownames from numbers to letters

convert_sppnames_toletters <-  function(mat) {
  
  
  rownames(mat) <- chartr("0123456789", "ABCDEFGHIJ", rownames(mat))
  colnames(mat) <- chartr("0123456789", "ABCDEFGHIJ", colnames(mat))
  
  return(mat)
  
}



# convert species names of ancestry table from numbers to letters

change_sppnames_letters_ancdist.table <- function(mat) {
  
  mat[,"ancestor"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"ancestor"])
  
  mat[,"spp"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"spp"])
  
  return(mat)
  
}




# compute dissimilarity matrix computting jaccard index between interacting species

compute_jaccard <- function(matrix, distance) {
  
  mat <- as.matrix(vegdist(matrix[-which(rowSums(matrix) == 0),], method = "jac", binary = T))
  
  return(mat)
  
}


compute_chisq <- function(matrix, distance) {
  
  mat <- as.matrix(vegdist(matrix[-which(rowSums(matrix) == 0),], method = "chisq", binary = T))
  
  return(mat)
  
}


compute_chisq_pred <- function(matrix, distance) {
  
  #Eliminate those species with no interactions (none as predators & none as preys)
  
  if(length(which(rowSums(matrix)==0 & colSums(matrix)==0)) > 0){
    
    mat <- matrix[-which(rowSums(matrix)==0 & colSums(matrix) ==0),
    -which(rowSums(matrix)==0 & colSums(matrix)==0)]

#compute dista ce between spp(columns) as predators

mat_dist <- as.matrix(vegdist(mat, method = distance, binary = T))
    
  }else{
    mat_dist <- as.matrix(vegdist(matrix, method = distance, binary = T))
  }
  
  return(mat_dist)
}


compute_chisq_prey <- function(matrix, distance) {
  
  #Eliminate those species with no interactions (none as predators & none as preys)
  
  if(length(which(rowSums(matrix)==0 & colSums(matrix)==0)) > 0){
    
    mat <- matrix[-which(rowSums(matrix)==0 & colSums(matrix) ==0),
                  -which(rowSums(matrix)==0 & colSums(matrix)==0)]
    
    #compute dista ce between spp(columns) as predators
    
    mat_dist <- as.matrix(vegdist(t(mat), method = distance, binary = T))
    
  }else{
    mat_dist <- as.matrix(vegdist(t(matrix), method = distance, binary = T))
  }
  
  return(mat_dist)
}



# Obtain tree from anc_table 


############################### v1 - Francois

sister.group<-function(sis.names,sis.dist){
  n<-length(sis.dist)
  distances<-paste0(rep(":",n),sis.dist)
  sis<-paste0(sis.names,distances)
  res<-paste0(sis,collapse=",")
  res<-paste0("(",res,")")
  res
}

ToPhylo<-function(data){
  k<-dim(data)[1]
  data.2<-data
  data.2$repr<-data$spp
  sisters<-levels(as.factor(data$spp))
  mothers<-levels(as.factor(data$ancestor))
  tips<-setdiff(sisters,mothers)
  root<-setdiff(mothers,sisters)
  data.2[k+1,"spp"] <- root
  data.2[k+1,"ancestor"] <- root
  data.2[k+1,"distance"] <- 0
  foc.nodes<-unique(data[which(data$spp%in%tips),"ancestor"])
  n<-length(foc.nodes)
  data.2$repr[data.2$spp%in%tips]<-data.2$repr[data.2$spp%in%tips]
  while(n>1){
    #print(foc.nodes)
    foc.nodes2<-unique(data.2[which(data.2$spp%in%foc.nodes),"ancestor"])
    for(i in 1:n){
      daughters<-data.2[which(data.2$ancestor==foc.nodes[i]),"repr"]
      #print(daughters)
      daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[i]),"distance"]
      data.2$repr[data.2$spp==foc.nodes[i]]<-paste0(sister.group(daughters,daughters.dist),foc.nodes[i])
    }
    tips<-foc.nodes
    foc.nodes<-foc.nodes2
    n<-length(foc.nodes)
  }
  #print(foc.nodes)
  data.2<-data.2[-(k+1),]
  daughters<-data.2[which(data.2$ancestor==foc.nodes[1]),"repr"]
  #print(daughters)
  daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[1]),"distance"]
  paste0(sister.group(daughters,daughters.dist),root)
}

############################## v2

# sister.group<-function(sis.names,sis.dist){
#   n<-length(sis.dist)
#   distances<-paste0(rep(":",n),sis.dist)
#   sis<-paste0(sis.names,distances)
#   res<-paste0(sis,collapse=",")
#   res<-paste0("(",res,")")
#   res
# }
# 
# ToPhylo<-function(data){
#   data.2<-data
#   data.2$repr<-data$spp
#   sisters<-levels(as.factor(data$spp))
#   mothers<-levels(as.factor(data$ancestor))
#   tips<-setdiff(sisters,mothers)
#   root<-setdiff(mothers,sisters)
#   foc.nodes<-unique(data[which(data$spp%in%tips),"ancestor"])
#   n<-length(foc.nodes)
#   data.2$repr[data.2$spp%in%tips]<-data.2$repr[data.2$spp%in%tips]
#   while(n>1){
#     foc.nodes2<-unique(data.2[which(data.2$spp%in%foc.nodes),"ancestor"])
#     for(i in 1:n){
#       daughters<-data.2[which(data.2$ancestor==foc.nodes[i]),"repr"]
#       #print(daughters)
#       daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[i]),"distance"]
#       data.2$repr[data.2$spp==foc.nodes[i]]<-paste0(sister.group(daughters,daughters.dist),foc.nodes[i])
#     }
#     tips<-foc.nodes
#     foc.nodes<-foc.nodes2
#     n<-length(foc.nodes)
#   }
#   daughters<-data.2[which(data.2$ancestor==foc.nodes[1]),"repr"]
#   #print(daughters)
#   daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[1]),"distance"]
#   paste0(sister.group(daughters,daughters.dist),root)
# }

############################## v3



# sister.group<-function(sis.names,sis.dist){
#   n<-length(sis.dist)
#   distances<-paste0(rep(":",n),sis.dist)
#   sis<-paste0(sis.names,distances)
#   res<-paste0(sis,collapse=",")
#   res<-paste0("(",res,")")
#   res
# }
# 
# ToPhylo<-function(data){
#   k<-dim(data)[1]
#   data.2<-data
#   data.2$repr<-data$spp
#   sisters<-levels(as.factor(data$spp))
#   mothers<-levels(as.factor(data$ancestor))
#   tips<-setdiff(sisters,mothers)
#   root<-setdiff(mothers,sisters)
#   data.2[k+1,"spp"] <- root
#   data.2[k+1,"ancestor"] <- root
#   data.2[k+1,"distance"] <- 0
#   foc.nodes<-unique(data[which(data$spp%in%tips),"ancestor"])
#   n<-length(foc.nodes)
#   data.2$repr[data.2$spp%in%tips]<-data.2$repr[data.2$spp%in%tips]
#   while(n>1){
#     #print(foc.nodes)
#     foc.nodes2<-unique(data.2[which(data.2$spp%in%foc.nodes),"ancestor"])
#     for(i in 1:n){
#       daughters<-data.2[which(data.2$ancestor==foc.nodes[i]),"repr"]
#       #print(daughters)
#       daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[i]),"distance"]
#       data.2$repr[data.2$spp==foc.nodes[i]]<-paste0(sister.group(daughters,daughters.dist),foc.nodes[i])
#     }
#     tips<-foc.nodes
#     foc.nodes<-foc.nodes2
#     n<-length(foc.nodes)
#   }
#   #print(foc.nodes)
#   data.2<-data.2[-(k+1),]
#   daughters<-data.2[which(data.2$ancestor==foc.nodes[1]),"repr"]
#   #print(daughters)
#   daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[1]),"distance"]
#   paste0(sister.group(daughters,daughters.dist),root)
# }





## Order alhpabetically rows and columns of dataframe

order_col.row_names <- function(df) {
  
  df_ordered <- df[order(names(df)) , order(names(df))]
  
  return(df_ordered)
  
}




#####################
#Compute NMI
####################


## NMI is for similarity. This function makes it a distance metric (1 - MNI)

nmi_to_distance <- function(value) {
  
  result <- 1 - value
  
  return(result)
  
}


compute_nmi_igraph_pred <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- compare(
      matrix[,id[i,1]],
      matrix[,id[i,2]],
      method = "nmi")
    
    matrix_nmi[id[i,2], id[i,1]] <- compare(
      matrix[,id[i,1]],
      matrix[,id[i,2]],
      method = "nmi")
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}



compute_nmi_igraph_prey <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- compare(
      matrix[id[i,1],],
      matrix[id[i,2],],
      method = "nmi")
    
    matrix_nmi[id[i,2], id[i,1]] <- compare(
      matrix[id[i,1],],
      matrix[id[i,2],],
      method = "nmi")
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}


compute_jaccard_pred <- function(matrix) {
  
  
  
  
  
}



compute_nmi_aricode_pred <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    # upper triangle
    
    matrix_nmi[id[i,1], id[i,2]] <- NMI(
      matrix[,id[i,1]],
      matrix[,id[i,2]])
    
    # lower triangle
    
    matrix_nmi[id[i,2], id[i,1]] <- matrix_nmi[id[i,1], id[i,2]]
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}




compute_nmi_aricode_prey <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    # upper triangle
    
    matrix_nmi[id[i,1], id[i,2]] <- NMI(
      matrix[id[i,1],],
      matrix[id[i,2],])
    
    # lower triangle
    
    matrix_nmi[id[i,2], id[i,1]] <- matrix_nmi[id[i,1], id[i,2]]
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}





# Compute mean matrix from a pair of matrix 

compute_mean_two_mat_from_list<- function(list){
  y<-apply(array(unlist(list), c(dim(list[[1]]), dim(list[[2]]), length(list))), 
           c(1,2), mean)
  colnames(y)<-colnames(list[[1]])
  rownames(y)<-rownames(list[[1]])
  return(y)
}



## Check that colnames (number and order) coincide in two lists


check_colnames_lists <- function(list1, list2) {
  
  vec_check <- c()
  
  for (i in 1:length(list1)) {
    
    if(length(which((colnames(list1[[i]]) == colnames(list2[[i]])) == FALSE)) == 0) {
      
      vec_check[i] <- "good"
      
    } else {
      
      vec_check[i] <- "PROBLEM"
      
      
    }
    
  }
  
  return(vec_check)
  
}


# Set diagonal elements of matrix to 1



diag_to0 <- function(matrix) {
  
  diag(matrix)<-0
  
  return(matrix)
  
}



# convert NaN to 1


convert_nan_to_0 <- function(vector) {
  
  vector[is.nan(vector)] <- 0
  
  return(vector)
  
}



# convert nan to 1 for all columns or rows


convet_nan_to_0_matrix <- function(matrix, marg) {
  
  m <- apply(matrix,
        MARGIN = marg,
        FUN = convert_nan_to_0)
  
  return(m)
  
}


check_all.finite <- function(matrix) {
  
  result <- all( is.finite( matrix ) )
  
  return(result)
  
}




# eliminate_basals <- function(mat, nbasals) {
#   #cat("Checking eliminate_basals: class =", class(mat), " | dim =", dim(mat), "\n")
#    matrix <- mat[(nbasals+1):nrow(mat),]
#   
#    return(matrix)
# }

eliminate_basals <- function(mat, nbasals) {
  # Skip if matrix is NULL or too small
  # if (is.null(mat) || nrow(mat) <= nbasals || ncol(mat) == 0) {
  #   return(NULL)
  # }
  return(mat[(nbasals + 1):nrow(mat), , drop = FALSE])
}




































