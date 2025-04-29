
## PROCRUSTES SIGNAL


compute_procrustes_signal_roles <- function(sdv_matrices, threshold = 0.8) {
  list_svd_pred <- sdv_matrices$list_svd_pred
  list_svd_prey <- sdv_matrices$list_svd_prey
  list_svd_both <- sdv_matrices$list_svd_both
  list_svd_eigen_phy <- sdv_matrices$list_svd_eigen.phy
  list_phylo_corr <- sdv_matrices$list_phylo.corr_cropped
  
  roles <- c("pred", "prey", "both")
  results_all <- list()
  
  for (role in roles) {
    cat("Processing role:", role, "\n")
    
    if (role == "pred") trait_list <- list_svd_pred
    if (role == "prey") trait_list <- list_svd_prey
    if (role == "both") trait_list <- list_svd_both
    
    timestep <- vector()
    S <- vector()
    cor <- vector()
    d_phylo <- vector()
    d_network <- vector()
    
    for (i in seq_along(trait_list)) {
      cat("  step", i, "\n")
      
      traits_i <- trait_list[[i]]
      
      if (is.null(traits_i) || nrow(traits_i) < 3) {
        cat("Warning: less than 3 species at step", i, "for role", role, "- skipping\n")
        next
      }
      
      phylo_i <- list_phylo_corr[[i]]
      phylo_i <- phylo_i[rownames(traits_i), rownames(traits_i)]
      
      eig_phylo <- eigen(phylo_i, symmetric = TRUE)
      eig_vectors_phylo <- eig_phylo$vectors
      
      # Compute cumulative variance explained
      singular_values_traits <- svd(traits_i)$d
      cum_var_traits <- cumsum(singular_values_traits^2) / sum(singular_values_traits^2)
      
      eigenvalues_phylo <- eig_phylo$values
      cum_var_phylo <- cumsum(eigenvalues_phylo) / sum(eigenvalues_phylo)
      
      num_axes_traits <- suppressWarnings(min(which(cum_var_traits >= threshold)))
      num_axes_phylo <- suppressWarnings(min(which(cum_var_phylo >= threshold)))
      
      if (is.infinite(num_axes_traits) || is.infinite(num_axes_phylo)) {
        cat("Warning: no valid axes at step", i, "- skipping\n")
        next
      }
      
      num_axes_to_keep <- max(num_axes_traits, num_axes_phylo)
      
      traits_kept <- traits_i[, 1:num_axes_to_keep, drop = FALSE]
      phylo_kept <- eig_vectors_phylo[, 1:num_axes_to_keep, drop = FALSE]
      
      proc <- protest(traits_kept, phylo_kept)
      
      timestep <- c(timestep, i)
      S <- c(S, nrow(traits_i))
      cor <- c(cor, proc$t0)
      d_phylo <- c(d_phylo, num_axes_phylo)
      d_network <- c(d_network, num_axes_traits)
    }
    
    if (length(timestep) > 0) {
      results_all[[role]] <- data.frame(
        role = role,
        timestep = timestep,
        S = S,
        cor = cor,
        d_phylo = d_phylo,
        d_network = d_network
      )
    } else {
      cat("Warning: no valid steps for role", role, "- skipping this role\n")
    }
  }
  
  if (length(results_all) == 0) {
    warning("Warning: no valid roles processed - returning empty data.frame")
    return(data.frame())
  }
  
  bind_rows(results_all)
}


# MANTEL SIGNAL


compute_mantel_signal_roles <- function(sdv_matrices, threshold = 0.8) {
  list_svd_pred <- sdv_matrices$list_svd_pred
  list_svd_prey <- sdv_matrices$list_svd_prey
  list_svd_both <- sdv_matrices$list_svd_both
  list_svd_eigen_phy <- sdv_matrices$list_svd_eigen.phy
  list_phylo_corr <- sdv_matrices$list_phylo.corr_cropped
  
  roles <- c("pred", "prey", "both")
  results_all <- list()
  
  for (role in roles) {
    cat("Processing role:", role, "\n")
    
    if (role == "pred") trait_list <- list_svd_pred
    if (role == "prey") trait_list <- list_svd_prey
    if (role == "both") trait_list <- list_svd_both
    
    timestep <- vector()
    S <- vector()
    cor_mantel <- vector()
    d_phylo <- vector()
    d_network <- vector()
    
    for (i in seq_along(trait_list)) {
      cat("  step", i, "\n")
      
      traits_i <- trait_list[[i]]
      
      if (is.null(traits_i) || nrow(traits_i) < 3) {
        cat("Warning: less than 3 species at step", i, "for role", role, "- skipping\n")
        next
      }
      
      phylo_i <- list_phylo_corr[[i]]
      phylo_i <- phylo_i[rownames(traits_i), rownames(traits_i)]
      
      eig_phylo <- eigen(phylo_i, symmetric = TRUE)
      eig_vectors_phylo <- eig_phylo$vectors
      
      # Compute cumulative variance explained
      singular_values_traits <- svd(traits_i)$d
      cum_var_traits <- cumsum(singular_values_traits^2) / sum(singular_values_traits^2)
      
      eigenvalues_phylo <- eig_phylo$values
      cum_var_phylo <- cumsum(eigenvalues_phylo) / sum(eigenvalues_phylo)
      
      num_axes_traits <- suppressWarnings(min(which(cum_var_traits >= threshold)))
      num_axes_phylo <- suppressWarnings(min(which(cum_var_phylo >= threshold)))
      
      if (is.infinite(num_axes_traits) || is.infinite(num_axes_phylo)) {
        cat("Warning: no valid axes at step", i, "- skipping\n")
        next
      }
      
      num_axes_to_keep <- max(num_axes_traits, num_axes_phylo)
      
      traits_kept <- traits_i[, 1:num_axes_to_keep, drop = FALSE]
      phylo_kept <- eig_vectors_phylo[, 1:num_axes_to_keep, drop = FALSE]
      
      # Mantel test
      dist_pred <- dist(traits_kept)
      dist_phylo <- dist(phylo_kept)
      
      mantel_result <- vegan::mantel(dist_phylo, dist_pred, permutations = 999)
      
      timestep <- c(timestep, i)
      S <- c(S, nrow(traits_i))
      cor_mantel <- c(cor_mantel, abs(mantel_result$statistic))
      d_phylo <- c(d_phylo, num_axes_phylo)
      d_network <- c(d_network, num_axes_traits)
    }
    
    if (length(timestep) > 0) {
      results_all[[role]] <- data.frame(
        role = role,
        timestep = timestep,
        S = S,
        cor_mantel = cor_mantel,
        d_phylo = d_phylo,
        d_network = d_network
      )
    } else {
      cat("Warning: no valid steps for role", role, "- skipping this role\n")
    }
  }
  
  if (length(results_all) == 0) {
    warning("Warning: no valid roles processed - returning empty data.frame")
    return(data.frame())
  }
  
  bind_rows(results_all)
}
