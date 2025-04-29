## Function to compute distance matrices

# Taking as input the results from the simulation, this function does the following:
#   
#   1) Prepares data
# - Identifies and removes timesteps where phylogenetic distances can't be computed.
# - Converts species identifiers from numbers to letters for consistency.
# 
# 2) Processes phylogenetic data:
# 
# - Constructs phylogenetic distance matrices for living species at each timestep.
# - Ensures species consistency between phylogenetic and presence matrices.
# - Egeinvectors of the cropped phylogenetic correlation matrix are computed using:
# list_svd_eigen.phy[[i]] <- eigen(list_phylo.corr_cropped[[i]], symmetric = T)$vec
# 
# 3) Processes interaction Data:
# 
# - Filters interaction matrices to include only present species.
# - Performs Singular Value Decomposition (SVD) on interaction matrices and keeps the right singular vectors (predator roles of species). The transposed right singular vectors are kept.
# 
# 
# It saves the computed eigenvectors from the phylogenetic correlation matricesand the ones from the interaction matrices. We will use this output to compute correlations between them.





compute_decomposed_matrices <- function(results_simulation, int, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # Number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  # Identify invalid timesteps (fewer than 3 species alive)
  non.valid_timesteps_phylo_distance <- which(rowSums(presence_matrix) < 3)
  final.discarded_timestep <- (non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]) + 1
  
  # Keep only valid timesteps
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep + 1):length(results_simulation$list_anc_dist)]
  network_list_full <- results_simulation$network_list[(final.discarded_timestep + 1):length(results_simulation$list_anc_dist)]
  presence_matrix <- presence_matrix[(final.discarded_timestep + 1):length(results_simulation$list_anc_dist), ]
  
  # --- Checks ---
  if (length(which(is.null(network_list_full))) > 0) {
    cat("PROBLEM - null network somewhere\n")
  }
  if (length(list_anc_dist) != length(network_list_full)) {
    cat("PROBLEM - list_anc_dist and network_list_full have different lengths\n")
  }
  if (nrow(presence_matrix) != length(network_list_full)) {
    cat("PROBLEM - presence_matrix rows and network_list_full length mismatch\n")
  }
  # ---------------
  
  # === Main decomposition ===
  
  list_svd_eigen.phy <- list()
  list_phylo.corr_cropped <- list()
  list_net_present_spp.letters <- list()
  
  for (i in seq_along(list_anc_dist)) {
    cat("step", i, "\n")
    
    alive_species <- names(which(presence_matrix[i, ] == 1))
    
    ## Build full phylogeny (including dead species)
    newick <- ToPhylo2(list_anc_dist[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root", ";", newick_tail))
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    # Identify missing species
    missing_species <- setdiff(alive_species, colnames(phylo.corr))
    if (length(missing_species) > 0) {
      cat("Warning: Some alive species not found in phylo tree at step", i, "\n")
      cat("Missing species:", paste(missing_species, collapse = ", "), "\n")
      alive_species <- setdiff(alive_species, missing_species)
    }
    
    # Subset network to alive species
    network_alive <- network_list_full[[i]][alive_species, , drop = FALSE]
    network_alive <- network_alive[, alive_species, drop = FALSE]
    
    # Save network with **basals kept** for SVD
    list_net_present_spp.letters[[i]] <- network_alive
    
    # Subset phylogeny to alive species (excluding basals later)
    list_phylo.corr_cropped[[i]] <- phylo.corr[alive_species, alive_species]
    
    # Save eigen decomposition of phylogeny
    list_svd_eigen.phy[[i]] <- eigen(list_phylo.corr_cropped[[i]], symmetric = TRUE)$vectors
  }
  
  # === Identify basal species ===
  all_species <- unique(unlist(lapply(list_net_present_spp.letters, rownames)))
  is_basal <- grepl("^Basal", all_species)
  basal_species <- all_species[is_basal]
  
  # === SVD decomposition of interaction matrices ===
  list_svd_pred <- list()
  for (i in seq_along(list_net_present_spp.letters)) {
    svd_result <- svd(list_net_present_spp.letters[[i]])
    kept_axes <- ncol(list_net_present_spp.letters[[i]])
    V_kept <- svd_result$v[, 1:kept_axes]
    list_svd_pred[[i]] <- t(V_kept)
  }
  
  # === Consistency check ===
  if (length(list_svd_pred) != length(list_svd_eigen.phy)) {
    cat("PROBLEM - length mismatch between interaction and phylo matrices\n")
  }
  
  # === Species lifespan summary ===
  all_anc_dist <- results_simulation$list_anc_dist
  all_species <- unique(unlist(lapply(all_anc_dist, function(x) x$spp)))
  
  species_life_summary <- data.frame(
    spp_name = all_species,
    timestep_birth = NA,
    timestep_extinct = NA,
    n_descendants = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(species_life_summary))) {
    spp_name <- species_life_summary$spp_name[i]
    if (spp_name %in% colnames(presence_matrix)) {
      presence_vector <- presence_matrix[, spp_name]
      if (any(presence_vector == 1)) {
        species_life_summary$timestep_birth[i] <- which(presence_vector == 1)[1]
        species_life_summary$timestep_extinct[i] <- tail(which(presence_vector == 1), 1)
      }
    }
  }
  
  # === Descendant counts ===
  all_relationships <- do.call(rbind, lapply(all_anc_dist, function(x) x[, c("spp", "ancestor")]))
  all_relationships <- all_relationships[all_relationships$ancestor != "0", ]
  
  unique_descendants <- split(all_relationships$spp, all_relationships$ancestor)
  desc_count <- sapply(unique_descendants, function(x) length(unique(x)))
  
  species_life_summary$n_descendants <- ifelse(
    species_life_summary$spp_name %in% names(desc_count),
    desc_count[species_life_summary$spp_name],
    0
  )
  
  # === Return all elements ===
  list(
    list_svd_pred = list_svd_pred,
    list_svd_eigen.phy = list_svd_eigen.phy,
    list_net_present_spp.letters = list_net_present_spp.letters,
    list_phylo.corr_cropped = list_phylo.corr_cropped,
    species_life_summary = species_life_summary,
    presence_matrix = presence_matrix,
    basal_species = basal_species # Save basal species explicitly
  )
}




compute_decomposed_matrices_roles <- function(results_simulation, int, Smax, nbasals) {
  presence_matrix <- results_simulation$presence_matrix
  n_steps <- length(results_simulation$network_list)
  
  non.valid_timesteps <- which(rowSums(presence_matrix) < 3)
  final.discarded_timestep <- (non.valid_timesteps[length(non.valid_timesteps)]) + 1
  
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep + 1):n_steps]
  network_list_full <- results_simulation$network_list[(final.discarded_timestep + 1):n_steps]
  presence_matrix <- presence_matrix[(final.discarded_timestep + 1):n_steps, ]
  
  list_svd_eigen.phy <- list()
  list_phylo.corr_cropped <- list()
  list_net_present_spp.letters <- list()
  
  list_svd_pred <- list()
  list_svd_prey <- list()
  list_svd_both <- list()
  
  for (i in seq_along(list_anc_dist)) {
    cat("step", i, "\n")
    
    alive_species <- names(which(presence_matrix[i, ] == 1))
    
    newick <- ToPhylo2(list_anc_dist[[i]])
    tree <- read.tree(text = sub("A root", ";", paste(newick, "root")))
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    missing_species <- setdiff(alive_species, colnames(phylo.corr))
    alive_species <- setdiff(alive_species, missing_species)
    
    network_alive <- network_list_full[[i]][alive_species, alive_species, drop = FALSE]
    # Identify basal species
    basal_species <- alive_species[grepl("^Basal", alive_species)]
    
    list_net_present_spp.letters[[i]] <- network_alive
    
    list_phylo.corr_cropped[[i]] <- phylo.corr[alive_species, alive_species]
    list_svd_eigen.phy[[i]] <- eigen(list_phylo.corr_cropped[[i]], symmetric = TRUE)$vectors
    
    # Decompose
    svd_result <- svd(network_alive)
    U <- svd_result$u
    V <- svd_result$v
    
    # Correct computation of Ssqrt
    Ssqrt <- diag(sqrt(svd_result$d))
    
    # Traits
    traits_pred <- U %*% Ssqrt
    traits_prey <- V %*% Ssqrt
    traits_both <- (traits_pred + traits_prey) / 2
    
    # Row names for matching
    rownames(traits_pred) <- rownames(network_alive)
    rownames(traits_prey) <- rownames(network_alive)
    rownames(traits_both) <- rownames(network_alive)
    
    # Remove basal species after
    species_nonbasal <- setdiff(rownames(network_alive), basal_species)
    
    list_svd_pred[[i]] <- traits_pred[species_nonbasal, , drop = FALSE]
    list_svd_prey[[i]] <- traits_prey[species_nonbasal, , drop = FALSE]
    list_svd_both[[i]] <- traits_both[species_nonbasal, , drop = FALSE]
  }
  
  all_anc_dist <- results_simulation$list_anc_dist
  all_species <- unique(unlist(lapply(all_anc_dist, function(x) x$spp)))
  
  species_life_summary <- data.frame(
    spp_name = all_species,
    timestep_birth = NA,
    timestep_extinct = NA,
    n_descendants = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(species_life_summary))) {
    spp <- species_life_summary$spp_name[i]
    if (spp %in% colnames(presence_matrix)) {
      presence_vector <- presence_matrix[, spp]
      if (any(presence_vector == 1)) {
        species_life_summary$timestep_birth[i] <- which(presence_vector == 1)[1]
        species_life_summary$timestep_extinct[i] <- tail(which(presence_vector == 1), 1)
      }
    }
  }
  
  all_relationships <- do.call(rbind, lapply(all_anc_dist, function(x) x[, c("spp", "ancestor")]))
  all_relationships <- all_relationships[all_relationships$ancestor != "0", ]
  
  unique_descendants <- split(all_relationships$spp, all_relationships$ancestor)
  desc_count <- sapply(unique_descendants, function(x) length(unique(x)))
  
  species_life_summary$n_descendants <- ifelse(
    species_life_summary$spp_name %in% names(desc_count),
    desc_count[species_life_summary$spp_name],
    0
  )
  
  list(
    list_svd_pred = list_svd_pred,
    list_svd_prey = list_svd_prey,
    list_svd_both = list_svd_both,
    list_svd_eigen.phy = list_svd_eigen.phy,
    list_net_present_spp.letters = list_net_present_spp.letters,
    list_phylo.corr_cropped = list_phylo.corr_cropped,
    species_life_summary = species_life_summary,
    presence_matrix = presence_matrix
  )
}
