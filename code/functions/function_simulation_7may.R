
########################################
# Function to generate a new set of traits for ancestors
# Each species is characterized by a set of 3 traits: n, o and r
rand_traits_anc = function(pars) {
  with(as.list(pars),{
    # n = runif(1, 0, 1)
    n = runif(1, 0.1, 0.3)
    r = av_r
    o = runif(1,0,n)
    traits = c(n = n, r = r, o = o)
    traits
  })
}

########################################
rand_traits_mut <- function(traits_anc, pars, direction = "random") {
  with(as.list(c(traits_anc, pars)), {
    
    # Compute the alpha parameter for the beta distribution
    a <- beta_n * n / (1 - n)
    
    # Determine the direction of divergence
    if (direction == "random") {
      direction <- ifelse(runif(1) < 0.5, "greater", "lesser")    }
    
    if (direction == "greater") {
      # Generate mutant trait greater than ancestor trait
      n_m <- rbeta(1, shape1 = a, shape2 = beta_n)
      n_m <- n + abs(n_m - n) * 0.7  # Ensure n_m is on the upper side of n
      n_m <- min(n_m, 1)  # Make sure n_m does not exceed 1
    } else if (direction == "lesser") {
      # Generate mutant trait less than ancestor trait
      n_m <- rbeta(1, shape1 = a, shape2 = beta_n)
      # Ensure n_m is on the lower side of n and within a reasonable range
      n_m <- n - abs(n_m - n) * 0.7  # Adjust this factor (0.5) to control how much lesser
      n_m <- max(n_m, 0)  # Make sure n_m is not below 0
    }
    
    # Calculate the optimum (o_m) for the mutant
    o_m <- n_m / 2
    
    # Combine the new traits
    traits_mut <- c(n = n_m, r = r, o = o_m)
    
    traits_mut
  })
}

########################################
# # Function to compute the interaction network from a set of traits
#  get_L_mat = function(basal, pars, traits_mat) {
#    with(as.list(pars),{
#      L = matrix(0, nr = Smax+Sbasal, nc = Smax)
#      
#      # Lower boundary
#      low = traits_mat$o - traits_mat$r
#      low_mat = matrix(low, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)
#      
#      # Upper boundary
#      high = traits_mat$o + traits_mat$r
#      high_mat = matrix(high, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)	
#      S = nrow(traits_mat)
#      
#      # Matrix of niche positions
#      n_mat = matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
#      
#      # Add the basal species
#      n_basal = matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
#      n_mat = rbind(n_basal, n_mat)
#      
#      # Test interactions
#      L[n_mat > low_mat & n_mat < high_mat] = 1
#      if(Smax > 1) diag(L[(Sbasal+1):(Sbasal+Smax),]) = 0
#      L
#    })
#  }


get_L_mat <- function(basal, pars, traits_mat) {
  with(as.list(pars), {
    # Initialize the interaction matrix
    L <- matrix(0, nr = Smax + Sbasal, nc = Smax)
    
    # Lower and upper boundaries for niches
    low <- traits_mat$o - traits_mat$r
    low_mat <- matrix(low, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
    high <- traits_mat$o + traits_mat$r
    high_mat <- matrix(high, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
    S <- nrow(traits_mat)
    
    # Matrix of niche positions
    n_mat <- matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
    
    # Add the basal species
    n_basal <- matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
    n_mat <- rbind(n_basal, n_mat)
    
    # Combine basal species with other species as potential prey
    n_prey <- c(basal, traits_mat$n)  # Include basal species as prey
    n_predator <- traits_mat$n  # Only non-basal species are potential predators
    
    # Define the probability function (e.g., Gaussian probability)
    prob_interaction <- function(distance, sigma = 0.04) {
      exp(- (distance^2) / (2 * sigma^2))
    }
    
    # Calculate distances and interaction probabilities
    for (i in 1:(Smax + Sbasal)) {
      for (j in 1:Smax) {
        # Check for NA values before comparing
        if (!is.na(n_mat[i, j]) && !is.na(low_mat[i, j]) && !is.na(high_mat[i, j])) {
          if (n_mat[i, j] > low_mat[i, j] && n_mat[i, j] < high_mat[i, j]) {
            # Compute the distance from the optimal niche
            distance <- abs(n_mat[i, j] - traits_mat$o[j])
            # Compute probability of interaction
            interaction_prob <- prob_interaction(distance)
            
            # Apply penalty if prey has a higher niche position than the predator
            if (n_prey[i] > n_predator[j]) {
              interaction_prob <- interaction_prob * 0.1 
            }
            
            # Assign interaction based on probability
            L[i, j] <- rbinom(1, 1, interaction_prob)  # Binomial draw: 1 interaction with probability 'interaction_prob'
          }
        }
      }
    }
    
    # Set diagonal to 0 (no self-interaction)
    if (Smax > 1) diag(L[(Sbasal + 1):(Sbasal + Smax), ]) <- 0
    L
  })
}





########################################
# Function to compute the interactions of a given species
get_L_vec = function(basal, pars, traits_mat, traits_mut) {
  with(as.list(pars),{
    L_vec = numeric(Smax+Sbasal)
    
    # Lower boundary
    low = traits_mut["o"] - traits_mut["r"]
    
    # Upper boundary
    high = traits_mut["o"] + traits_mut["r"]
    
    # Vector of niche positions
    n_vec = c(basal, traits_mat$n)
    
    # Test interactions
    L_vec[n_vec > as.numeric(low) & n_vec < as.numeric(high)] = 1
    L_vec
  })
}


generate_species_names <- function(n) {
  letters1 <- LETTERS
  letters2 <- as.vector(outer(LETTERS, LETTERS, paste0))
  letters3 <- as.vector(outer(letters2, LETTERS, paste0))
  
  all_names <- c(letters1, letters2, letters3)
  
  if (n > length(all_names)) {
    stop("Requested more names than available.")
  }
  
  return(all_names[1:n])
}




sim_model <- function(seed, pars, nsteps) {
  with(pars, {
    set.seed(seed)
    
    # === Generate consistent species names ===
    species_names <- generate_species_names(Smax)
    
    # Draw producer traits
    basal <- runif(Sbasal, 0, 0.005)
    
    # Draw initial species traits
    traits_mat <- matrix(nrow = Smax, ncol = 3)
    traits_mat[1, ] <- rand_traits_anc(pars)
    traits_mat <- as.data.frame(traits_mat)
    names(traits_mat) <- c("n", "r", "o")
    rownames(traits_mat) <- species_names
    
    # Initialize presence/absence matrix
    pres <- matrix(0, nrow = nsteps, ncol = Smax)
    colnames(pres) <- species_names
    rownames(pres) <- 1:nsteps
    pres[1, 1] <- 1
    
    # Initialize ancestry tables
    anc <- data.frame(step = rep(NA, Smax), ancestor = rep(NA, Smax), descendant = rep(NA, Smax))
    rownames(anc) <- species_names
    
    dist_anc <- data.frame(
      spp = species_names,
      ancestor = rep(NA, Smax),
      A_E = rep(NA, Smax),
      distance = rep(0, Smax)
    )
    dist_anc$ancestor[1] <- NA
    dist_anc$A_E[1] <- "A"
    
    list_dist_anc <- list()
    L_list <- list()
    L_cropped_list <- list()
    df_guild <- data.frame(step = integer(), npredators = integer(), nherbivores = integer(), nomnivores = integer())
    
    S <- 1
    
    for (step in 2:nsteps) {
      ActualS <- sum(pres[step - 1, ])
      extinct_trait_disp <- rep(FALSE, Smax)
      
      for (i in 1:Smax) {
        if (S >= Smax) break
        if (pres[step - 1, i] == 1 && dist_anc$A_E[i] != "E") {
          dist_anc$A_E[i] <- "A"
          dist_anc$distance[i] <- dist_anc$distance[i] + 1
          pres[step, i] <- 1
          
          test_number <- runif(1)
          speciation_prob <- u_max / (1 + exp(d * (ActualS - I_max)))
          
          if (test_number < speciation_prob) {
            traits_anc <- traits_mat[i, ]
            first_mutant_direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
            traits_mut1 <- rand_traits_mut(traits_anc, pars, direction = first_mutant_direction)
            
            I <- get_L_vec(basal, pars, traits_mat, traits_mut1)
            sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
            
            estab_prob <- ifelse(sum_I > 0,
                                 SN * estab_prob_neutral + (1 - SN) * (u_0pos + u_1pos * exp(-a_upos * sum_I)),
                                 0
            )
            
            if (runif(1) < estab_prob) {
              S <- S + 1
              traits_mat[S, ] <- traits_mut1
              pres[step, S] <- 1
              anc[S, ] <- c(step, species_names[i], species_names[S])
              dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
              
              if (abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
                pres[step, i] <- 0
                dist_anc$A_E[i] <- "E"
                extinct_trait_disp[i] <- TRUE
              }
              if (S >= Smax) break
            }
            
            second_mutant_direction <- ifelse(first_mutant_direction == "greater", "lesser", "greater")
            traits_mut2 <- rand_traits_mut(traits_anc, pars, direction = second_mutant_direction)
            
            I <- get_L_vec(basal, pars, traits_mat, traits_mut2)
            sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
            
            estab_prob <- ifelse(sum_I > 0,
                                 SN * estab_prob_neutral + (1 - SN) * (u_0pos + u_1pos * exp(-a_upos * sum_I)),
                                 0
            )
            
            
            
            if (runif(1) < estab_prob) {
              S <- S + 1
              traits_mat[S, ] <- traits_mut2
              pres[step, S] <- 1
              anc[S, ] <- c(step, species_names[i], species_names[S])
              dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
              
              if (abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
                pres[step, i] <- 0
                dist_anc$A_E[i] <- "E"
                extinct_trait_disp[i] <- TRUE
              }
              if (S >= Smax) break
            }
          }
        }
      }
      
      if (S >= Smax) break
      
      # --- EXTINCTION DYNAMICS ---
      pres_vec <- pres[step, ]
      cooc <- outer(pres_vec, pres_vec, "*")
      L <- get_L_mat(basal, pars, traits_mat)
      L[(Sbasal + 1):(Sbasal + Smax), ] <- L[(Sbasal + 1):(Sbasal + Smax), ] * cooc
      rownames(L) <- c(paste0("Basal", 1:Sbasal), species_names)
      colnames(L) <- species_names
      L_list[[step]] <- L
      
      present_non_basal <- which(pres[step, ] == 1)
      
      if (length(present_non_basal) > 1) {
        # === GUILD CLASSIFICATION (needs basal rows) ===
        L_for_guilds <- L[c(1:Sbasal, Sbasal + present_non_basal), present_non_basal, drop=FALSE]
        rownames(L_for_guilds) <- c(paste0("Basal", 1:Sbasal), species_names[present_non_basal])
        colnames(L_for_guilds) <- species_names[present_non_basal]
        
        herbivores <- (colSums(L_for_guilds[1:Sbasal, ] > 0) > 0) & (colSums(L_for_guilds[(Sbasal + 1):nrow(L_for_guilds), ] > 0) == 0)
        predators <- (colSums(L_for_guilds[1:Sbasal, ] == 0) & (colSums(L_for_guilds[(Sbasal + 1):nrow(L_for_guilds), ] > 0) > 0))
        omnivores <- (colSums(L_for_guilds[1:Sbasal, ] > 0) > 0) & (colSums(L_for_guilds[(Sbasal + 1):nrow(L_for_guilds), ] > 0) > 0)
        
        df_guild <- rbind(df_guild, data.frame(step = step,
                                               npredators = sum(predators),
                                               nherbivores = sum(herbivores),
                                               nomnivores = sum(omnivores)))
        
        # === INTERACTION MATRIX FOR EXTINCTION: no basal ===
        idx_nonbasal_rows <- Sbasal + present_non_basal
        idx_nonbasal_cols <- present_non_basal
        L_cropped <- L[idx_nonbasal_rows, idx_nonbasal_cols, drop=FALSE]
        rownames(L_cropped) <- species_names[present_non_basal]
        colnames(L_cropped) <- species_names[present_non_basal]
        L_cropped_list[[step]] <- L_cropped
        
        svd_result <- svd(L_cropped)
        similarity_matrix <- svd_result$v %*% diag(svd_result$d) %*% t(svd_result$v)
        diag(similarity_matrix) <- 0
        avg_similarity <- rowMeans(similarity_matrix) * (ncol(similarity_matrix) / (ncol(similarity_matrix) - 1))
        avg_similarity[avg_similarity < exp(-5)] <- 0
        
        out_I_cropped <- rowSums(L_cropped)  # number of predators per prey
        ext_prob_topdown <- e_0neg + e_1neg * (1 - exp(-a_eneg * out_I_cropped))
        
        transformed_similarity <- 1 - exp(-competition_coefficient * avg_similarity)
        
        # --- identify predators/omnivores with no prey ---
        spp_noresources <- colSums(L_cropped) == 0  # no prey → predators & omnivores extinct
        
        # --- initialize extinction probability vector ---
        ext_prob_sel_full <- numeric(Smax)
        
        # === HERBIVORES extinction ===
        ext_prob_sel_full[present_non_basal[herbivores]] <- ext_prob_topdown[herbivores]
        
        # === PREDATORS extinction ===
        ext_prob_pred <- transformed_similarity[predators]
        ext_prob_pred[spp_noresources[predators]] <- 1
        ext_prob_sel_full[present_non_basal[predators]] <- ext_prob_pred
        
        # === OMNIVORES extinction ===
        ext_prob_omn <- beta_ext * ext_prob_topdown[omnivores] + (1 - beta_ext) * transformed_similarity[omnivores]
        ext_prob_omn[spp_noresources[omnivores]] <- 1
        ext_prob_sel_full[present_non_basal[omnivores]] <- ext_prob_omn
        
        # --- apply multiplier to selection-based probabilities ---
        ext_prob_sel_full[present_non_basal] <- ext_prob_sel_full[present_non_basal] * multiplier
        
        # --- blend neutral and selection-based extinction ---
        # Define neutral extinction probability
        ext_prob_neutral <- 0.005  
        
        # Apply the SN interpolation
        ext_prob_sel_full[present_non_basal] <- SN * ext_prob_neutral + 
          (1 - SN) * ext_prob_sel_full[present_non_basal]
        
        
        # --- apply extinction ---
        random_number <- runif(length(present_non_basal))
        pres[step, present_non_basal] <- ifelse(random_number < ext_prob_sel_full[present_non_basal], 0, 1)
      }
      
      
      
      
      for (i in 1:Smax) {
        if (!extinct_trait_disp[i] && !is.na(pres[step, i]) && pres[step, i] != pres[step - 1, i] && pres[step, i] == 0) {
          dist_anc$A_E[i] <- "E"
        }
      }
      
      list_dist_anc[[step]] <- na.omit(dist_anc)
      
      if (step == 10 && sum(pres[step, ]) == 0) return(NULL)
      if (step > 10 && sum(pres[step, ]) == 0) return(NULL)
      if (step == 50 && sum(pres[step, ]) < 10) return(NULL)
      if (step == 30 && sum(pres[step, ]) <= 3) return(NULL)
      if (step == nsteps && sum(pres[step, ]) < 40) return(NULL)
    }
    
    list(
      pres = pres,
      traits = traits_mat,
      anc = anc,
      L_list = L_list,
      L_cropped_list = L_cropped_list,
      basal = basal,
      dist_anc = dist_anc,
      list_dist_anc = list_dist_anc,
      df_guild = df_guild
    )
  })
}
