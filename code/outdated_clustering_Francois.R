## Method 3.1. fixed number of clusters


```{r}

spectral_clustering <- function(graph, nb_cluster, normalized = TRUE) {#from J. Chiquet's git page https://jchiquet.github.io/MAP566/docs/mixture-models/map566-lecture-graph-clustering-part1.html
  
  ## Compute Laplcian matrix
  L <- laplacian_matrix(graph, normalized = normalized) 
  ## Generates indices of last (smallest) K vectors
  selected <- rev(1:ncol(L))[1:nb_cluster] 
  ## Extract an normalized eigen-vectors
  U <- eigen(L)$vectors[, selected, drop = FALSE]  # spectral decomposition
  U <- sweep(U, 1, sqrt(rowSums(U^2)), '/')    
  ## Perform k-means
  res <- kmeans(U, nb_cluster, nstart = 40)$cl
  
  res
}

compute_optimal_cluster_metrics_spectral <- function(results_simulation, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # Number of timesteps
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
  ari_values <- numeric(length(list_networks_sppnames_letters))
  nmi_values <- numeric(length(list_networks_sppnames_letters))
  
  # Store trees and clusters
  original_trees <- list()
  pruned_trees <- list()
  phylo_clusters_list <- list()
  
  for (i in seq_along(list_anc_dist_letters)) {
    
    cat("Processing timestep:", i, "\n")
    
    # Process phylogenetic correlation matrix
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";", newick_tail))
    
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    # Save the full tree before pruning
    original_trees[[i]] <- tree
    
    # Get alive species
    alive_species <- names(which(presence_matrix[i, ] == 1))
    phylo_corr_cropped <- phylo.corr[alive_species, alive_species]
    
    # Prune the tree to only keep alive species
    pruned_tree <- keep.tip(tree, alive_species)
    pruned_trees[[i]] <- pruned_tree
    
    # Clustering interaction matrix using SBM
    interaction_matrix <- list_networks_sppnames_letters[[i]][alive_species, alive_species]
    
    diag(interaction_matrix) <- 0  
    
    if (sum(interaction_matrix) == 0) {
      cat("Sparse or empty interaction matrix at timestep:", i, "- Assigning all species to one cluster.\n")
      interaction_clusters <- rep(1, nrow(interaction_matrix))
    } else {
      sbm_fit <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli")
      interaction_clusters <- sbm_fit$memberships
    }
    
    # **Spectral Clustering for Phylogeny**
    num_species <- nrow(phylo_corr_cropped)
    
    if (num_species > 2) {
      
      # Convert correlation matrix into an adjacency matrix (similarity graph)
      adjacency_matrix <- abs(phylo_corr_cropped) # Ensuring positive similarities
      diag(adjacency_matrix) <- 0  # Remove self-connections
      
      # Convert adjacency matrix into a graph
      phylo_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
      
      # Select the number of clusters
      max_k <- min(4, num_species - 1)
      
      # Apply spectral clustering
      phylo_clusters <- spectral_clustering(phylo_graph, nb_cluster = max_k, normalized = TRUE)
      
    } else {
      phylo_clusters <- rep(1, num_species)
    }
    
    phylo_clusters_list[[i]] <- setNames(phylo_clusters, alive_species)  # Save cluster assignments
    
    # Compare clusters using ARI and NMI
    ari_values[i] <- compare(interaction_clusters, phylo_clusters, method = "adjusted.rand")
    nmi_values[i] <- compare(interaction_clusters, phylo_clusters, method = "nmi")
  }
  
  return(list(
    results = data.frame(
      timestep = 1:length(ari_values),
      ARI = ari_values,
      NMI = nmi_values
    ),
    original_trees = original_trees,
    pruned_trees = pruned_trees,
    clusters = phylo_clusters_list
  ))
}


```



```{r}
phylo_output <- compute_optimal_cluster_metrics_spectral(results_simulation = res_sim, Smax = 1000, nbasals = 5)

```



Is the clustering on the phylogeny working well?
  
  
  ```{r}
plot_phylo_clusters_at_timestep <- function(phylo_output, timestep) {
  
  # Extract pruned tree and clustering results
  tree <- phylo_output$pruned_trees[[timestep]]
  phylo_clusters <- phylo_output$clusters[[timestep]]
  
  # Ensure species names are properly mapped
  species_names <- names(phylo_clusters)  # Extract species names from clusters
  tree_species <- tree$tip.label          # Extract species names from the pruned tree
  
  # Sort species names alphabetically to ensure consistent order
  species_names <- sort(species_names)
  tree_species <- sort(tree_species)
  
  # Match clusters to tree species explicitly
  matched_clusters <- phylo_clusters[match(tree_species, species_names)]
  
  # Ensure there are no mismatches
  if (any(is.na(matched_clusters))) {
    stop("Error: Some species in the pruned tree do not have matching cluster assignments.")
  }
  
  # Debugging: Print if species names still mismatch
  if (any(is.na(matched_clusters))) {
    warning("Some species in the pruned tree do not have matching cluster assignments!")
    cat("Mismatched species:", tree_species[is.na(matched_clusters)], "\n")
    
    # Remove species with missing assignments
    valid_idx <- !is.na(matched_clusters)
    tree_species <- tree_species[valid_idx]
    matched_clusters <- matched_clusters[valid_idx]
  }
  
  # Prepare data for plotting
  cluster_data <- data.frame(label = tree_species, cluster = as.factor(matched_clusters))
  
  # Plot the tree with colored clusters
  ggtree(tree) %<+% cluster_data +  # Attach cluster data to tree
    geom_tippoint(aes(color = cluster), size = 3) +
    theme_minimal() +
    labs(title = paste("Phylogenetic Clusters at Timestep", timestep),
         color = "Cluster")
}

# Run the function for a given timestep
plot_phylo_clusters_at_timestep(phylo_output, timestep = 20)
```



Plot results


```{r}
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


```



## Mehtod 3.2. Unsupervised n clusters



```{r}
laplacian_spectral_gap<- function(graph){
  L <- laplacian_matrix(graph, normalized = TRUE)
  comps<-count_components(graph)
  lambdas <- sort(eigen(L)$values)
  l<-length(lambdas)
  s_gaps<-lambdas[-1]-lambdas[-l]
  s_util<-s_gaps[-(1:comps)]
  s_util<-s_util[1:round(l/2)]
  opt_n<-which.max(s_util)+comps
  
  # par(mfrow=c(2,1))
  # plot(lambdas,xlab="",ylab="lambda",type="l")
  # plot(s_gaps,xlab="",ylab="spectral gap",type="l")
  
  list("spectral_gaps"=s_gaps,"optim_n"=opt_n)
}

```


```{r}
spectral_clustering <- function(graph, normalized = TRUE) {
  ## Compute Laplacian matrix
  L <- laplacian_matrix(graph, normalized = normalized)
  
  ## Check if L contains NA/NaN/Inf
  if (any(is.na(L) | is.infinite(L))) {
    warning("Laplacian matrix contains NA/Inf values. Returning single cluster.")
    return(list(clusters = rep(1, nrow(L)), k = 1))
  }
  
  ## Compute Optimal Number of Clusters (`k`) Using Spectral Gap
  spectral_results <- laplacian_spectral_gap(graph)
  optimal_k <- spectral_results$optim_n
  
  ## Ensure `k` is valid
  if (is.na(optimal_k) || optimal_k < 2) {
    optimal_k <- 2  # Default to at least 2 clusters
  }
  
  ## Perform Eigen Decomposition
  eig <- eigen(L, symmetric = TRUE)
  eig_vectors <- eig$vectors
  
  ## Extract eigenvectors associated with `optimal_k`
  U <- eig_vectors[, 1:optimal_k, drop = FALSE]  
  U <- sweep(U, 1, apply(U, 1, norm, "2"), '/')  # Normalize rows
  
  ## Ensure no NaN or Inf before k-means
  if (any(is.na(U) | is.infinite(U))) {
    warning("Eigenvectors contain NA/Inf. Returning single cluster.")
    return(list(clusters = rep(1, nrow(U)), k = 1))
  }
  
  ## Perform k-means clustering
  set.seed(123)
  kmeans_result <- tryCatch({
    kmeans(U, centers = optimal_k, nstart = 40)
  }, error = function(e) {
    warning("K-means failed. Returning single cluster.")
    return(list(cluster = rep(1, nrow(U))))
  })
  
  return(list(clusters = kmeans_result$cluster, k = optimal_k))
}

```



```{r}
compute_optimal_cluster_metrics_spectral <- function(results_simulation, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # Number of timesteps
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
  ari_values <- numeric(length(list_networks_sppnames_letters))
  nmi_values <- numeric(length(list_networks_sppnames_letters))
  
  # Store trees and clusters
  original_trees <- list()
  pruned_trees <- list()
  phylo_clusters_list <- list()
  num_clusters_list <- numeric(length(list_anc_dist_letters))  # Store detected number of clusters
  
  for (i in seq_along(list_anc_dist_letters)) {
    
    cat("Processing timestep:", i, "\n")
    
    # Process phylogenetic correlation matrix
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";", newick_tail))
    
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    # Save the full tree before pruning
    original_trees[[i]] <- tree
    
    # Get alive species
    alive_species <- colnames(presence_matrix)[which(presence_matrix[i, ] == 1)]
    phylo_corr_cropped <- phylo.corr[alive_species, alive_species]
    
    # Prune the tree to only keep alive species
    pruned_tree <- keep.tip(tree, alive_species)
    pruned_trees[[i]] <- pruned_tree
    
    # Clustering interaction matrix using SBM
    interaction_matrix <- list_networks_sppnames_letters[[i]][alive_species, alive_species]
    
    diag(interaction_matrix) <- 0  
    
    if (sum(interaction_matrix) == 0) {
      cat("Sparse or empty interaction matrix at timestep:", i, "- Assigning all species to one cluster.\n")
      interaction_clusters <- rep(1, nrow(interaction_matrix))
    } else {
      sbm_fit <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli")
      interaction_clusters <- sbm_fit$memberships
    }
    
    # **Spectral Clustering for Phylogeny**
    num_species <- nrow(phylo_corr_cropped)
    
    if (num_species > 2) {
      
      # Convert correlation matrix into an adjacency matrix (similarity graph)
      adjacency_matrix <- abs(phylo_corr_cropped) # Ensuring positive similarities
      diag(adjacency_matrix) <- 0  # Remove self-connections
      
      # Convert adjacency matrix into a graph
      phylo_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
      
      cat("Before clustering, alive species:", alive_species, "\n")
      
      # Apply spectral clustering with automatic `k` selection
      spectral_result <- spectral_clustering(phylo_graph, normalized = TRUE)
      phylo_clusters <- spectral_result$clusters
      num_clusters_list[i] <- length(phylo_clusters)  # Store detected number of clusters
      
    } else {
      phylo_clusters <- rep(1, num_species)
      num_clusters_list[i] <- 1
    }
    
    # Ensure all alive species have an assigned cluster
    matched_clusters <- rep(NA, length(alive_species))
    names(matched_clusters) <- alive_species
    
    # Assign clusters to species that exist in clustering results
    assigned_species <- names(phylo_clusters)  # Species present in clustering results
    matched_clusters[assigned_species] <- phylo_clusters[assigned_species]
    
    # Handle missing species
    missing_species <- is.na(matched_clusters)
    if (any(missing_species)) {
      warning("Some species were missing in clustering results at timestep", i)
      
      # Assign each missing species to its own singleton cluster
      next_cluster_id <- max(phylo_clusters, na.rm = TRUE) + 1
      matched_clusters[missing_species] <- seq(next_cluster_id, next_cluster_id + sum(missing_species) - 1)
    }
    
    # Save corrected cluster assignments
    phylo_clusters_list[[i]] <- matched_clusters
    
    # Compare clusters using ARI and NMI
    ari_values[i] <- compare(interaction_clusters, phylo_clusters, method = "adjusted.rand")
    nmi_values[i] <- compare(interaction_clusters, phylo_clusters, method = "nmi")
  }
  
  cat("After clustering, phylo_clusters names:", names(phylo_clusters), "\n")
  
  return(list(
    results = data.frame(
      timestep = 1:length(ari_values),
      ARI = ari_values,
      NMI = nmi_values,
      num_clusters = num_clusters_list
    ),
    original_trees = original_trees,
    pruned_trees = pruned_trees,
    clusters = phylo_clusters_list
  ))
}

```


```{r}
phylo_output <- compute_optimal_cluster_metrics_spectral(results_simulation = res_sim, Smax = 1000, nbasals = 5)

```



```{r}
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


```