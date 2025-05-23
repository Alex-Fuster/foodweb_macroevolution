
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
    #  cat("  step", i, "\n")
      
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
    #  cat("  step", i, "\n")
      
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



# JACCARD SIGNAL

compute_jaccard_phylo_correlation <- function(sdv_matrices) {
  list_network <- sdv_matrices$list_net_present_spp.letters
  list_corrphylo <- sdv_matrices$list_phylo.corr_cropped
  
  results <- list()
  
  for (i in seq_along(list_network)) {
    #cat("Processing step", i, "\n")
    
    interaction_matrix <- list_network[[i]]
    phylo_corr <- list_corrphylo[[i]]
    
    common_species <- intersect(rownames(interaction_matrix), rownames(phylo_corr))
    if (length(common_species) < 3) {
      cat("Warning: less than 3 species at step", i, "- skipping\n")
      next
    }
    
    interaction_matrix <- interaction_matrix[common_species, common_species, drop = FALSE]
    phylo_corr <- phylo_corr[common_species, common_species, drop = FALSE]
    
    jaccard_pred <- 1 - as.matrix(vegan::vegdist(t(interaction_matrix), method = "jaccard", binary = TRUE))
    jaccard_prey <- 1 - as.matrix(vegan::vegdist(interaction_matrix, method = "jaccard", binary = TRUE))
    jaccard_both <- (jaccard_pred + jaccard_prey) / 2
    
    upper_tri <- function(mat) mat[upper.tri(mat)]
    
    # Extract vectors
    pred_vec <- upper_tri(jaccard_pred)
    prey_vec <- upper_tri(jaccard_prey)
    both_vec <- upper_tri(jaccard_both)
    phylo_vec <- upper_tri(phylo_corr)
    
    # Check if enough non-NA pairs exist
    if (sum(complete.cases(pred_vec, phylo_vec)) < 3 ||
        sum(complete.cases(prey_vec, phylo_vec)) < 3 ||
        sum(complete.cases(both_vec, phylo_vec)) < 3) {
      cat("Warning: too few complete pairs at step", i, "- skipping\n")
      next
    }
    
    # Compute Spearman correlations
    df_step <- data.frame(
      timestep = i,
      S = length(common_species),
      predator = suppressWarnings(cor(pred_vec, phylo_vec, use = "complete.obs", method = "spearman")),
      prey = suppressWarnings(cor(prey_vec, phylo_vec, use = "complete.obs", method = "spearman")),
      both = suppressWarnings(cor(both_vec, phylo_vec, use = "complete.obs", method = "spearman"))
    )
    
    results[[length(results) + 1]] <- df_step
  }
  
  if (length(results) == 0) {
    warning("No valid steps, returning empty data frame")
    return(data.frame())
  }
  
  dplyr::bind_rows(results)
}


# CLUSTERS SIGNAL

compute_clustering_metrics_roles <- function(sdv_matrices, results_simulation, Smax = 1000, nbasals = 5) {
  
  list_networks_sppnames_letters <- sdv_matrices$list_net_present_spp.letters
  list_basal_species <- sdv_matrices$list_basal_species
  presence_matrix <- results_simulation$presence_matrix
  list_anc_dist <- results_simulation$list_anc_dist
  
  VI <- numeric()
  NMI <- numeric()
  timestep <- numeric()
  
  valid_timesteps <- sdv_matrices$list_valid_timesteps
  valid_timesteps <- valid_timesteps[valid_timesteps <= length(list_anc_dist)]
  
  for (i in seq_along(valid_timesteps))  {  # <-- only loop over existing matrices
    
    step <- valid_timesteps[i]
    cat("Processing timestep:", step, "(index", i, ")\n")
    
    
    basal_species <- list_basal_species[[i]]
    
    newick <- ToPhylo2(list_anc_dist[[step]])
    tree <- read.tree(text = sub("A root", ";", paste(newick, "root")))
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    
    alive_species_network <- rownames(list_networks_sppnames_letters[[i]])
    alive_species_tree <- intersect(tree$tip.label, alive_species_network)
    
    if (length(alive_species_tree) < 3) {
      cat("Warning: less than 3 species after matching at step", i, "- skipping\n")
     # ARI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    pruned_tree <- keep.tip(tree, alive_species_tree)
    
    optim_result <- tryCatch(
      laplacian_spectral_gap_tree(pruned_tree),
      error = function(e) list(spectral_gaps = NA, optim_n = NA)
    )
    
    # Defensive extraction
    optim_n <- if (!is.null(optim_result) && "optim_n" %in% names(optim_result) && !is.null(optim_result$optim_n)) optim_result$optim_n else NA
    
    if (is.na(optim_n) || length(optim_n) != 1 || !is.numeric(optim_n) || optim_n < 2) {
      cat("Warning: invalid or NA optim_n at step", i, "- skipping\n")
    #  ARI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    phylo_clusters <- spectral_clustering_tree(pruned_tree, optim_n)
    names(phylo_clusters) <- pruned_tree$tip.label
    
    interaction_matrix <- list_networks_sppnames_letters[[i]]
    diag(interaction_matrix) <- 0
    
    if (all(is.na(interaction_matrix)) || nrow(interaction_matrix) < 3 || sum(interaction_matrix, na.rm = TRUE) == 0) {
      cat("Warning: too few interactions at step", i, "- skipping\n")
      VI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    # sbm_fit <- tryCatch(
    #   sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli"),
    #   error = function(e) NULL
    # )
    
    sbm_fit <- tryCatch({
      suppressMessages(
        suppressWarnings({
          capture.output(
            model <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli"),
            file = NULL
          )
          model  # return the model object, not the output
        })
      )
    }, error = function(e) NULL)
    
    
    if (is.null(sbm_fit) || !("memberships" %in% names(sbm_fit))) {
      # Optional: if (!quiet) cat("Warning: SBM fitting failed or invalid object at step", i, "- skipping\n")
      VI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    interaction_clusters <- sbm_fit$memberships
    
    names(interaction_clusters) <- rownames(interaction_matrix)
    
    # Match species, remove basals
    species_to_keep <- intersect(names(phylo_clusters), names(interaction_clusters))
    species_to_keep <- setdiff(species_to_keep, basal_species)
    
    if (length(species_to_keep) < 3) {
      cat("Warning: less than 3 nonbasal overlapping species at step", i, "- skipping\n")
      VI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    # Filter cluster assignments to shared species
    # Extract cluster assignments for shared species
    phylo_clusters_nonbasal <- phylo_clusters[species_to_keep]
    interaction_clusters_nonbasal <- interaction_clusters[species_to_keep]
    
    # Defensive: keep only species with non-NA clusters in both
    valid_idx <- complete.cases(phylo_clusters_nonbasal, interaction_clusters_nonbasal)
    
    if (sum(valid_idx) < 3) {
      VI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    # Convert to consecutive integer factors starting from 1
    phylo_vec <- as.integer(factor(phylo_clusters_nonbasal[valid_idx]))
    interact_vec <- as.integer(factor(interaction_clusters_nonbasal[valid_idx]))
    
    # Extra safety: ensure all cluster indices are <= n (this is where your crash happens)
    n_vals <- length(phylo_vec)
    if (any(phylo_vec > n_vals) || any(interact_vec > n_vals)) {
      cat("Warning: Invalid cluster indices at step", i, "- skipping\n")
      VI[i] <- NA
      NMI[i] <- NA
      timestep[i] <- i
      next
    }
    
    # Now safe to compare
    # Adjusted Rand Index from mclust
   # ARI[i] <- adjustedRandIndex(interact_vec, phylo_vec)
    
    # Normalized Mutual Information from clue
    # NMI[i] <- cl_agreement(as.cl_partition(interact_vec),
    #                        as.cl_partition(phylo_vec),
    #                        method = "NMI")
    # 
    # VI[i] <- cl_agreement(as.cl_partition(interact_vec),
    #                       as.cl_partition(phylo_vec),
    #                       method = "VI")
    
    
    
    # Convert to igraph communities
    # comm_phylo <- make_clusters(graph.empty(n = length(phylo_vec)), membership = phylo_vec)
    # comm_interact <- make_clusters(graph.empty(n = length(interact_vec)), membership = interact_vec)
    
    # Defensive checks
    if (length(phylo_vec) != length(interact_vec)) {
      cat("Warning: Mismatched vector lengths at step", step, "- skipping\n")
      NMI[i] <- NA
      VI[i] <- NA
      timestep[i] <- step
      next
    }
    if (max(phylo_vec) > length(phylo_vec) || max(interact_vec) > length(interact_vec)) {
      cat("Warning: Cluster ID exceeds vector length at step", step, "- skipping\n")
      NMI[i] <- NA
      VI[i] <- NA
      timestep[i] <- step
      next
    }
    
    n_species <- length(phylo_vec)
    
    # Use correct number of vertices
    comm_phylo <- tryCatch({
      make_clusters(graph.empty(n = n_species), membership = phylo_vec)
    }, error = function(e) {
      cat("make_clusters error (phylo) at step", step, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    comm_interact <- tryCatch({
      make_clusters(graph.empty(n = n_species), membership = interact_vec)
    }, error = function(e) {
      cat("make_clusters error (interact) at step", step, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(comm_phylo) || is.null(comm_interact)) {
      NMI[i] <- NA
      VI[i] <- NA
      timestep[i] <- step
      next
    }
    
    # Compute metrics
    NMI[i] <- compare(comm_phylo, comm_interact, method = "nmi")
    VI[i] <- compare(comm_phylo, comm_interact, method = "vi")
    
    
    
    
    
    timestep[i] <- step
    
  }
  
  data.frame(timestep = timestep, 
             VI = VI, 
             NMI = NMI)
}





laplacian_spectral_gap_tree <- function(tree, normalized = FALSE, correlation = TRUE) {
  tree.corr <- tree
  tree.corr$edge.length <- sapply(tree.corr$edge.length, function(x) ifelse(x == 0, 1e-5, x))
  
  tree.vcv <- tryCatch(vcv(tree.corr, corr = correlation), error = function(e) NULL)
  
  if (is.null(tree.vcv) || !is.matrix(tree.vcv) || is.na(nrow(tree.vcv)) || nrow(tree.vcv) < 3) {
    cat("Warning: tree.vcv is null or invalid at this step\n")
    return(list("spectral_gaps" = NA, "optim_n" = NA))
  }
  
  diag(tree.vcv) <- 0
  
  if (normalized) {
    L <- diag(rep(1, nrow(tree.vcv))) - diag(1 / rowSums(tree.vcv)) %*% tree.vcv
  } else {
    L <- diag(rowSums(tree.vcv)) - tree.vcv
  }
  
  lambdas <- sort(eigen(L)$values)
  comps <- length(which(lambdas < 1e-12))
  l <- length(lambdas)
  s_gaps <- lambdas[-1] - lambdas[-l]
  s_util <- s_gaps[-(1:comps)]
  s_util <- s_util[1:round(l / 4)]
  opt_n <- which.max(s_util) + comps
  
  list("spectral_gaps" = s_gaps, "optim_n" = opt_n)
}




spectral_clustering_tree <- function(tree, nb_cluster, normalized = FALSE, correlation = TRUE) {  
  tree.corr <- tree
  tree.corr$edge.length <- sapply(tree.corr$edge.length, function(x) ifelse(x == 0, 1e-5, x))
  
  tree.vcv <- vcv(tree.corr, corr = correlation)
  diag(tree.vcv) <- 0
  
  # Skip if too few species
  if (nrow(tree.vcv) < 3 || nb_cluster >= nrow(tree.vcv)) {
    return(rep(1, nrow(tree.vcv)))  # assign all species to the same cluster
  }
  
  if (normalized) {
    L <- diag(rep(1, nrow(tree.vcv))) - diag(1 / rowSums(tree.vcv)) %*% tree.vcv
  } else {
    L <- diag(rowSums(tree.vcv)) - tree.vcv
  }
  
  selected <- rev(1:ncol(L))[1:nb_cluster]
  U <- eigen(L)$vectors[, selected, drop = FALSE]
  U <- sweep(U, 1, sqrt(rowSums(U^2)), '/')
  
  res <- kmeans(U, nb_cluster, nstart = 40)$cluster
  
  return(res)
}


