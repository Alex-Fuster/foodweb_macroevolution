
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


species_names <- generate_species_names(1000)
species_names


########################################

# sim_model <- function(seed, pars, nsteps) {
#   with(pars, {
#     set.seed(seed)
#     
#     # === Generate consistent species names ===
#     species_names <- generate_species_names(Smax)
#     
#     # Debug: print species names and check for NA
#    # cat("Generated species names:\n")
#     #print(species_names)
#     if (any(is.na(species_names))) stop("❌ species_names contains NA!")
#     
#     # Draw producer traits
#     basal <- runif(Sbasal, 0, 0.005)
#     
#     # Initialize traits matrix
#     traits_mat <- matrix(nrow = Smax, ncol = 3)
#     traits_mat[1, ] <- rand_traits_anc(pars)
#     traits_mat <- as.data.frame(traits_mat)
#     names(traits_mat) <- c("n", "r", "o")
#     rownames(traits_mat) <- species_names
#     
#     # Initialize presence matrix
#     pres <- matrix(0, nrow = nsteps, ncol = Smax)
#     colnames(pres) <- species_names
#     rownames(pres) <- 1:nsteps
#     pres[1, 1] <- 1
#     
#     # Initialize ancestry tracking
#     dist_anc <- data.frame(
#       spp = species_names,
#       ancestor = rep(NA, Smax),
#       A_E = rep(NA, Smax),
#       distance = rep(0, Smax)
#     )
#     dist_anc$ancestor[1] <- NA
#     dist_anc$A_E[1] <- "A"
#     
#     list_dist_anc <- list()
#     network_list_full <- list()
#     
#     S <- 1  # initial species count
#     
#     for (step in 2:nsteps) {
#       ActualS <- sum(pres[step - 1, ])
#       extinct_trait_disp <- rep(FALSE, Smax)
#       
#       for (i in 1:Smax) {
#         if (S >= Smax) {
#           cat("⚠️ Reached maximum species capacity (Smax =", Smax, ") at step", step, "\n")
#           break
#         }
#         
#         if (pres[step - 1, i] == 1 && dist_anc$A_E[i] != "E") {
#           dist_anc$A_E[i] <- "A"
#           dist_anc$distance[i] <- dist_anc$distance[i] + 1
#           pres[step, i] <- 1
#           
#           # speciation attempt
#           test_number <- runif(1)
#           speciation_prob <- u_max / (1 + exp(d * (ActualS - I_max)))
#           
#           if (test_number < speciation_prob) {
#             traits_anc <- traits_mat[i, ]
#             first_mutant_direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
#             traits_mut1 <- rand_traits_mut(traits_anc, pars, direction = first_mutant_direction)
#             
#             I <- get_L_vec(basal, pars, traits_mat, traits_mut1)
#             sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
#             
#             estab_prob <- ifelse(sum_I > 0,
#                                  SN * estab_prob_neutral + (1 - SN) * (u_0pos + u_1pos * exp(-a_upos * sum_I)),
#                                  0)
#             
#             if (runif(1) < estab_prob) {
#               S <- S + 1
#               if (S > Smax) stop(paste("❌ S exceeded Smax at step", step, "S =", S, "Smax =", Smax))
#               traits_mat[S, ] <- traits_mut1
#               pres[step, S] <- 1
#               dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
#               
#               if (abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
#                 pres[step, i] <- 0
#                 dist_anc$A_E[i] <- "E"
#                 extinct_trait_disp[i] <- TRUE
#               }
#             }
#             
#             second_mutant_direction <- ifelse(first_mutant_direction == "greater", "lesser", "greater")
#             traits_mut2 <- rand_traits_mut(traits_anc, pars, direction = second_mutant_direction)
#             
#             I <- get_L_vec(basal, pars, traits_mat, traits_mut2)
#             sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
#             
#             estab_prob <- ifelse(sum_I > 0,
#                                  SN * estab_prob_neutral + (1 - SN) * (u_0pos + u_1pos * exp(-a_upos * sum_I)),
#                                  0)
#             
#             if (runif(1) < estab_prob) {
#               S <- S + 1
#               if (S > Smax) stop(paste("❌ S exceeded Smax at step", step, "S =", S, "Smax =", Smax))
#               traits_mat[S, ] <- traits_mut2
#               pres[step, S] <- 1
#               dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
#               
#               if (abs(traits_anc["n"] - traits_mut2["n"]) <= 0.05) {
#                 pres[step, i] <- 0
#                 dist_anc$A_E[i] <- "E"
#                 extinct_trait_disp[i] <- TRUE
#               }
#             }
#           }
#         }
#       }
#       
#       # Explicit extinction updates
#       for (i in 1:Smax) {
#         if (pres[step - 1, i] == 1 && pres[step, i] == 0) {
#           dist_anc$A_E[i] <- "E"
#         }
#       }
#       
#       # Build interaction matrix
#       pres_vec <- pres[step, ]
#       cooc <- outer(pres_vec, pres_vec, "*")
#       L <- get_L_mat(basal, pars, traits_mat)
#       L[(Sbasal + 1):(Sbasal + Smax), ] <- L[(Sbasal + 1):(Sbasal + Smax), ] * cooc
#       rownames(L) <- c(paste0("Basal", 1:Sbasal), species_names)
#       colnames(L) <- species_names
#       
#       # Determine alive non-basal species
#       alive_species <- species_names[pres[step, ] == 1]
#       rows_keep <- c(paste0("Basal", 1:Sbasal), alive_species)
#       cols_keep <- alive_species
#       
#       # Filter interaction matrix: keep basal rows + alive species rows + alive species cols
#       L_filtered <- L[rows_keep, cols_keep, drop = FALSE]
#       network_list_full[[step]] <- L_filtered
#       
#       # Debug check for NA rows marked A
#       na_rows <- which(is.na(dist_anc$spp) & dist_anc$A_E == "A")
#       if (length(na_rows) > 0) {
#         cat("⚠️ At step", step, "found", length(na_rows), "rows with spp==NA and A_E=='A'\n")
#         print(dist_anc[na_rows, ])
#       }
#       
#       list_dist_anc[[step]] <- dist_anc
#       
#      # cat("==== Step", step, " list_dist_anc snapshot ====\n")
#      # print(list_dist_anc[[step]])
#       
#      # cat("Row indices where A_E == 'A':\n")
#      # print(which(list_dist_anc[[step]]$A_E == "A"))
#       
#      # cat("spp values at those rows:\n")
#     #  print(list_dist_anc[[step]]$spp[which(list_dist_anc[[step]]$A_E == "A")])
#       
#       # Consistency check
#       alive_in_pres <- species_names[pres[step, ] == 1]
#       idx_alive <- which(list_dist_anc[[step]]$A_E == "A")
#       alive_in_dist_anc <- list_dist_anc[[step]]$spp[idx_alive]
#       
#     #  cat("alive_in_pres:", alive_in_pres, "\n")
#     #  cat("alive_in_dist_anc:", alive_in_dist_anc, "\n")
#       
#       if (!setequal(alive_in_pres, alive_in_dist_anc)) {
#         cat("⚠️ Mismatch at step", step, "\n")
#         cat("In pres but not dist_anc:", setdiff(alive_in_pres, alive_in_dist_anc), "\n")
#         cat("In dist_anc but not pres:", setdiff(alive_in_dist_anc, alive_in_pres), "\n")
#       }else {
#         cat("✅ Step", step, "passed consistency check.\n")
#       }
#       
#       if (step == 10 && sum(pres[step, ]) == 0) return(NULL)
#       if (step > 10 && sum(pres[step, ]) == 0) return(NULL)
#       if (step == 50 && sum(pres[step, ]) < 10) return(NULL)
#       if (step == 30 && sum(pres[step, ]) <= 3) return(NULL)
#       if (step == nsteps && sum(pres[step, ]) < 40) return(NULL)
#     }
#     
#     return(list(
#       presence_matrix = pres,
#       traits = traits_mat,
#       dist_anc = dist_anc,
#       list_dist_anc = list_dist_anc,
#       network_list_full = network_list_full
#     ))
#   })
# }


sim_model <- function(seed, pars, nsteps) {
  with(pars, {
    set.seed(seed)
    
    # === Generate species names ===
    species_names <- generate_species_names(Smax)
    if (any(is.na(species_names))) stop("❌ species_names contains NA!")
    
    # Producer traits
    basal <- runif(Sbasal, 0, 0.005)
    
    # Traits matrix
    traits_mat <- matrix(nrow = Smax, ncol = 3)
    traits_mat[1, ] <- rand_traits_anc(pars)
    traits_mat <- as.data.frame(traits_mat)
    names(traits_mat) <- c("n", "r", "o")
    rownames(traits_mat) <- species_names
    
    # Presence matrix
    pres <- matrix(0, nrow = nsteps, ncol = Smax)
    colnames(pres) <- species_names
    pres[1, 1] <- 1
    
    # Ancestry tracking
    dist_anc <- data.frame(
      spp = species_names,
      ancestor = rep(NA, Smax),
      A_E = rep(NA, Smax),
      distance = rep(0, Smax)
    )
    dist_anc$spp[1] <- species_names[1]
    dist_anc$ancestor[1] <- NA
    dist_anc$A_E[1] <- "A"
    
    # Store ancestry over time
    list_dist_anc <- list()
    
    # Networks
    L_list <- list()
    L_cropped_list <- list()
    
    # Guild dataframe
    df_guild <- data.frame(step = integer(), npredators = integer(), nherbivores = integer(), nomnivores = integer())
    
    S <- 1
    
    for (step in 2:nsteps) {
      ActualS <- sum(pres[step - 1, ])
      
      # --- SPECIES UPDATE ---
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
            
            # --- First mutant ---
            dir1 <- ifelse(runif(1) < 0.5, "greater", "lesser")
            traits_mut1 <- rand_traits_mut(traits_anc, pars, direction = dir1)
            I <- get_L_vec(basal, pars, traits_mat, traits_mut1)
            sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
            estab_prob_sel <- ifelse(sum_I == 0, 0, u_0pos + u_1pos * exp(-a_upos * sum_I))
            estab_prob <- SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
            
            if (runif(1) < estab_prob) {
              S <- S + 1
              if (S > Smax) stop("S exceeded Smax")
              traits_mat[S, ] <- traits_mut1
              pres[step, S] <- 1
              dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
              if (abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
                pres[step, i] <- 0
                dist_anc$A_E[i] <- "E"
              }
            }
            
            # --- Second mutant ---
            dir2 <- ifelse(dir1 == "greater", "lesser", "greater")
            traits_mut2 <- rand_traits_mut(traits_anc, pars, direction = dir2)
            I <- get_L_vec(basal, pars, traits_mat, traits_mut2)
            sum_I <- sum(I * c(rep(1, Sbasal), pres[step, ]))
            estab_prob_sel <- ifelse(sum_I == 0, 0, u_0pos + u_1pos * exp(-a_upos * sum_I))
            estab_prob <- SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
            
            if (runif(1) < estab_prob) {
              S <- S + 1
              if (S > Smax) stop("S exceeded Smax")
              traits_mat[S, ] <- traits_mut2
              pres[step, S] <- 1
              dist_anc[S, ] <- list(species_names[S], species_names[i], "A", 1)
              if (abs(traits_anc["n"] - traits_mut2["n"]) <= 0.05) {
                pres[step, i] <- 0
                dist_anc$A_E[i] <- "E"
              }
            }
          }
        }
      }
      
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
        L_cropped <- L[c(1:Sbasal, Sbasal + present_non_basal), present_non_basal]
        rownames(L_cropped) <- c(paste0("Basal", 1:Sbasal), species_names[present_non_basal])
        colnames(L_cropped) <- species_names[present_non_basal]
        L_cropped_list[[step]] <- L_cropped
        
        
        svd_result <- svd(L_cropped)
        sim_mat <- svd_result$v %*% diag(svd_result$d) %*% t(svd_result$v)
        diag(sim_mat) <- 0
        avg_sim <- rowMeans(sim_mat) * (ncol(sim_mat) / (ncol(sim_mat) - 1))
        avg_sim[avg_sim < exp(-5)] <- 0
        
        out_I_cropped <- rowSums(L_cropped)[(Sbasal + 1):nrow(L_cropped)]
        ext_prob_topdown <- e_0neg + e_1neg * (1 - exp(-a_eneg * out_I_cropped))
        transformed_sim <- 1 - exp(-competition_coefficient * avg_sim)
        ext_prob_sel <- beta_ext * ext_prob_topdown + (1 - beta_ext) * transformed_sim
        ext_prob_sel_full <- numeric(Smax)
        ext_prob_sel <- ext_prob_sel * multiplier
        ext_prob_sel_full[present_non_basal] <- ext_prob_sel
        
        herbivores <- (colSums(L_cropped[1:Sbasal, ] > 0) > 0) & (colSums(L_cropped[(Sbasal + 1):nrow(L_cropped), ] > 0) == 0)
        predators <- (colSums(L_cropped[1:Sbasal, ] > 0) == 0) & (colSums(L_cropped[(Sbasal + 1):nrow(L_cropped), ] > 0) > 0)
        omnivores <- (colSums(L_cropped[1:Sbasal, ] > 0) > 0) & (colSums(L_cropped[(Sbasal + 1):nrow(L_cropped), ] > 0) > 0)
        df_guild <- rbind(df_guild, data.frame(step, npredators = sum(predators), nherbivores = sum(herbivores), nomnivores = sum(omnivores)))
        
        ext_prob <- SN * ext_prob_neutral + (1 - SN) * ext_prob_sel_full
        test_extprob <- runif(Smax)
        pres[step, present_non_basal] <- ifelse(test_extprob[present_non_basal] < ext_prob[present_non_basal], 0, 1)
      }
      
      # Update dist_anc extinction status
      dist_anc$A_E[pres[step, ] == 0] <- "E"
      
      
      # --- Consistency check ---
      idx_alive <- which(dist_anc$A_E == "A")
      alive_in_dist_anc <- dist_anc$spp[idx_alive]
      alive_in_pres <- species_names[pres[step, ] == 1]
      if (!setequal(alive_in_pres, alive_in_dist_anc)) {
        cat("⚠️ Mismatch at step", step, "\n")
      }
      
      dist_anc_clean <- na.omit(dist_anc)
      dist_anc_clean$distance <- as.numeric(dist_anc_clean$distance)
      list_dist_anc[[step]] <- dist_anc_clean
      
      
      # Early exits
      if (step == 10 && sum(pres[step, ]) == 0) return(NULL)
      if (step > 10 && sum(pres[step, ]) == 0) return(NULL)
      if (step == 50 && sum(pres[step, ]) < 10) return(NULL)
      if (step == 30 && sum(pres[step, ]) <= 3) return(NULL)
      if (step == nsteps && sum(pres[step, ]) < 40) return(NULL)
    }
    
    return(list(
      pres = pres,
      traits = traits_mat,
     # anc = anc,
      L_list = L_list,
      L_cropped_list = L_cropped_list,
      basal = basal,
      dist_anc = dist_anc,
      list_dist_anc = list_dist_anc,
      df_guild = df_guild
    ))
    
  })
}



