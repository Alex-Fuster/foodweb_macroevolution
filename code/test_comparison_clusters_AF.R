library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(igraph)
library(EnvStats)
library(aricode)
library(vegan)
library(ade4)
library(ape)
library(tidyr)
library(ggraph)
library(cluster)    
library(fpc)  
library(ggtree)
library(gridExtra)

source(here::here("functions_simulation.R")) 
source(here::here("function_to_phylo.R")) 
source(here::here("functions_simulation_inspection.R"))
source(here::here("function_compute_decomposed_matrices.R"))


## Overview of simulation restults for 1 simulation

res_sim <- readRDS(here::here("results_toplot_networks_timesteps.rds"))

nsteps =  141
inspect_simulation_fw(simulation_data = res_sim, nbasals = 5, Smax = 1000)



#### Compute decomposed matrices

#The function *compute_decomposed_matrices()* performs Singular Value Decomposition (SVD) on interaction matrices and phylogenetic distance matrices. It saves the egeinvectors of the phylogenetic distance matrix, and the transposed right singular vectors of the interaction matrix (predator roles of species).


sdv_matrices <- compute_decomposed_matrices(results_simulation = res_sim,
                                            int = "foodweb",
                                            Smax = 1000,
                                            nbasals = 5)

list_svd_pred <- sdv_matrices$list_svd_pred
list_svd_eigen.phy <- sdv_matrices$list_svd_eigen.phy

list_network <- sdv_matrices$list_net_present_spp.letters
list_corrphylo <- sdv_matrices$list_phylo.corr_cropped



##############################
# Method 3 - compare clusters
##############################

spectral_clustering_tree <- function(tree, nb_cluster, normalized=FALSE, correlation = TRUE) {  
  # Ensure non-zero branch lengths
  tree.corr <- tree
  tree.corr$edge.length <- sapply(tree.corr$edge.length, function(x) ifelse(x == 0, 1e-5, x))
  
  # Compute variance-covariance matrix
  tree.vcv <- vcv(tree.corr, corr = correlation)
  diag(tree.vcv) <- 0  # Remove self-similarities
  
  # Compute Laplacian matrix
  if (normalized) {
    L <- diag(rep(1, dim(tree.vcv)[1])) - diag(1 / rowSums(tree.vcv)) %*% tree.vcv
  } else {
    L <- diag(rowSums(tree.vcv)) - tree.vcv
  }
  
  ## Generates indices of last (smallest) K vectors
  selected <- rev(1:ncol(L))[1:nb_cluster] 
  ## Extract n normalized eigen-vectors
  U <- eigen(L)$vectors[, selected, drop = FALSE]  # spectral decomposition
  U <- sweep(U, 1, sqrt(rowSums(U^2)), '/')    
  ## Perform k-means
  res <- kmeans(U, nb_cluster, nstart = 40)$cl
  
  res
}


laplacian_spectral_gap_tree <- function(tree, normalized = FALSE, correlation = TRUE) {
  # Ensure non-zero branch lengths
  tree.corr <- tree
  tree.corr$edge.length <- sapply(tree.corr$edge.length, function(x) ifelse(x == 0, 1e-5, x))
  
  # Compute variance-covariance matrix
  tree.vcv <- vcv(tree.corr, corr = correlation)
  diag(tree.vcv) <- 0  # Remove self-similarities
  
  # Compute Laplacian matrix
  if (normalized) {
    L <- diag(rep(1, dim(tree.vcv)[1])) - diag(1 / rowSums(tree.vcv)) %*% tree.vcv
  } else {
    L <- diag(apply(tree.vcv,1,sum))-tree.vcv
  }
  
  # Compute eigenvalues
  lambdas <- sort(eigen(L)$values)
  comps <- length(which(lambdas < 1e-12))
  l <- length(lambdas)
  
  # Compute spectral gaps
  s_gaps <- lambdas[-1] - lambdas[-l]
  s_util <- s_gaps[-(1:comps)]
  s_util <- s_util[1:round(l / 2)]
  opt_n <- which.max(s_util) + comps
  
  # Plot eigenvalues & spectral gaps
  # par(mfrow = c(2, 1))
  # plot(lambdas, xlab = "", ylab = "Lambda", type = "l")
  # plot(s_gaps, xlab = "", ylab = "Spectral Gap", type = "l")
  
  return(list("spectral_gaps" = s_gaps, "optim_n" = opt_n))
}



compute_optimal_cluster_metrics_tree <- function(results_simulation, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  n_steps <- length(results_simulation$network_list)
  
  # Identify valid timesteps
  non.valid_timesteps_phylo_distance <- which(rowSums(presence_matrix) < 3)
  final.discarded_timestep <- ifelse(length(non.valid_timesteps_phylo_distance) > 0, 
                                     max(non.valid_timesteps_phylo_distance) + 1, 1)
  
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep + 1):n_steps]
  network_list <- results_simulation$network_list[(final.discarded_timestep + 1):n_steps]
  
  # Convert species names to letters
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  colnames(presence_matrix) <- seq(1:Smax)
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  presence_matrix <- presence_matrix[(final.discarded_timestep + 1):n_steps, ]
  
  # Store results
  skipped_timesteps <- c()  
  original_trees <- list()
  pruned_trees <- list()
  phylo_clusters_list <- list()
  num_clusters_list <- numeric(length(list_anc_dist_letters))
  ari_values <- numeric(length(list_anc_dist_letters))
  nmi_values <- numeric(length(list_anc_dist_letters))
  
  for (i in seq_along(list_anc_dist_letters)) {
    
    cat("Processing timestep:", i, "\n")
    
    # Process phylogenetic tree
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";", newick_tail))
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    
    # Save full tree before pruning
    original_trees[[i]] <- tree
    
    # Get alive species
    alive_species <- colnames(presence_matrix)[which(presence_matrix[i, ] == 1)]
    alive_species_filtered <- intersect(alive_species, tree$tip.label)
    
    # Check for unmatched species
    dropped_species <- setdiff(alive_species, tree$tip.label)
    if (length(dropped_species) > 0) {
      warning("Some species in alive_species are not in the tree at timestep ", i, ": ", paste(dropped_species, collapse = ", "))
    }
    
    # Prune the tree only using the filtered species
    pruned_tree <- keep.tip(tree, alive_species_filtered)
    pruned_trees[[i]] <- pruned_tree
    
    # Compute optimal number of clusters using spectral gap
    optim_n <- laplacian_spectral_gap_tree(pruned_tree)$optim_n
    
    # Handle cases where `optim_n` is invalid
    if (is.null(optim_n) || length(optim_n) == 0 || is.na(optim_n) || optim_n < 2) {
      warning("Skipping timestep ", i, " due to invalid optim_n value: ", optim_n)
      skipped_timesteps <- c(skipped_timesteps, i)
      num_clusters_list[i] <- NA
      ari_values[i] <- NA
      nmi_values[i] <- NA
      next  # Skip this timestep
    }
    
    # Apply spectral clustering on the pruned tree
    phylo_clusters <- spectral_clustering_tree(pruned_tree, optim_n)
    num_clusters_list[i] <- length(unique(phylo_clusters))
    
    # Assign cluster labels to species
    names(phylo_clusters) <- pruned_tree$tip.label
    phylo_clusters_list[[i]] <- phylo_clusters
    
    # **Clustering the Interaction Matrix (SBM)**
    interaction_matrix <- list_networks_sppnames_letters[[i]][alive_species_filtered, alive_species_filtered]
    diag(interaction_matrix) <- 0  
    
    if (sum(interaction_matrix) == 0) {
      cat("Sparse or empty interaction matrix at timestep:", i, "- Assigning all species to one cluster.\n")
      interaction_clusters <- rep(1, length(alive_species_filtered))
    } else {
      sbm_fit <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli")
      interaction_clusters <- sbm_fit$memberships
    }
    
    # Compute ARI and NMI
    ari_values[i] <- compare(interaction_clusters, phylo_clusters, method = "adjusted.rand")
    nmi_values[i] <- compare(interaction_clusters, phylo_clusters, method = "nmi")
  }
  
  return(list(
    results = data.frame(
      timestep = 1:length(num_clusters_list),
      num_clusters = num_clusters_list,
      ARI = ari_values,
      NMI = nmi_values
    ),
    original_trees = original_trees,
    pruned_trees = pruned_trees,
    clusters = phylo_clusters_list,
    skipped_timesteps = skipped_timesteps  
  ))
}






phylo_output <- invisible(compute_optimal_cluster_metrics_tree(results_simulation = res_sim, Smax = 1000, nbasals = 5))





# Is the clustering on the phylogeny working well?

plot_phylo_clusters_at_timestep <- function(phylo_output, timestep) {
  
  # Extract pruned tree and clustering results
  tree <- phylo_output$pruned_trees[[timestep]]
  phylo_clusters <- phylo_output$clusters[[timestep]]
  
  # Ensure species names are properly mapped
  species_names <- names(phylo_clusters)  # Extract species names from clusters
  tree_species <- tree$tip.label          # Extract species names from the pruned tree
  
  # *Check which species in the tree are also in the cluster assignments**
  common_species <- intersect(tree_species, species_names)
  
  # If no common species exist, print a warning and return
  if (length(common_species) == 0) {
    warning("No matching species found between tree and phylo_clusters at timestep", timestep)
    return(NULL)
  }
  
  # **Subset the clustering results to match only species that exist in the tree**
  matched_clusters <- phylo_clusters[common_species]
  
  # *Map cluster IDs to the correct tree nodes**
  cluster_data <- data.frame(label = common_species, cluster = as.factor(matched_clusters))
  
  # *Remove unused cluster levels to avoid extra legend categories**
  cluster_data$cluster <- droplevels(cluster_data$cluster)
  
  #*Ensure colors match the actual number of clusters**
  unique_clusters <- unique(cluster_data$cluster)
  num_clusters <- length(unique_clusters)
  cluster_colors <- scales::hue_pal()(num_clusters)  # Automatically generate distinct colors
  
  # *Plot the tree and color nodes based on their clusters**
  ggtree(tree) %<+% cluster_data +  # Attach cluster data to tree
    geom_tippoint(aes(color = cluster), size = 3) +
    scale_color_manual(values = cluster_colors) +  # Use correct number of colors
    theme_minimal() +
    labs(title = paste("Phylogenetic Clusters at Timestep", timestep),
         color = "Cluster")
}




plot_phylo_clusters_at_timestep(phylo_output, timestep = 20)





#n clusters with time:


ggplot(phylo_output$results, aes(x = timestep, y = num_clusters)) +
  geom_point(color = "black", alpha = 0.7, size = 3) +
  geom_smooth(se = FALSE, method = "gam", span = 0.3, alpha = 0.8, color = "black") +
  labs(
    x = "Time Step", 
    y = "Number of Clusters", 
    title = "Detected Number of Clusters Over Time"
  ) +
  theme_minimal()+
  ylim(0,max(phylo_output$results$num_clusters)+2)



# Correlation interaction - phylogeny clusters


filtered_results <- phylo_output$results %>%
  filter(abs(ARI) > 1e-10, abs(NMI) > 1e-10)

ggplot(filtered_results, aes(x = timestep)) +
  geom_point(aes(y = ARI), color = 'blue', alpha = 0.7, size = 3) +
  geom_smooth(aes(y = ARI), se = FALSE, method = "gam", span = 0.3, alpha = 0.8, color = "blue") +
  geom_point(aes(y = NMI), color = 'red', alpha = 0.7, size = 3) +
  geom_smooth(aes(y = NMI), se = FALSE, method = "gam", span = 0.3, alpha = 0.8, color = "red") +
  labs(
    x = "Time Step", 
    y = "Clustering Metrics", 
    title = "ARI and NMI Over Time"
  ) +
  theme_minimal() +
  # scale_y_continuous(limits = c(0, 1), name = "Metric Value") +
  scale_x_continuous(name = "Time Step") +
  theme(legend.position = "right")
