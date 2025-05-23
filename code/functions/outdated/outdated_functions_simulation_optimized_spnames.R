
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
      direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
    }
    
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

sim_model = function(seed, pars, nsteps) {
  with(pars, {
    set.seed(seed)
    
    # === Generate consistent species names ===
    species_names <- generate_species_names(Smax)
    
    # Draw the traits of the producers
    basal = runif(Sbasal, 0, 0.005)
    
    # Draw the first species traits
    traits_mat = matrix(nrow = Smax, ncol = 3)
    traits_mat[1,] = rand_traits_anc(pars)
    traits_mat = as.data.frame(traits_mat)
    names(traits_mat) = c("n", "r", "o")
    rownames(traits_mat) = species_names
    
    # Set the presence/absence matrix
    pres = matrix(0, nrow = nsteps, ncol = Smax)
    colnames(pres) <- species_names
    rownames(pres) <- 1:nsteps
    pres[1, 1] = 1  # The first species is present initially
    
    # Set the ancestry object
    anc = data.frame(step = rep(NA, Smax), ancestor = rep(NA, Smax), descendant = rep(NA, Smax))
    rownames(anc) <- species_names
    
    # Set the ancestry distance table
    dist_anc = data.frame(
      spp = species_names,
      ancestor = rep(NA, Smax),
      A_E = rep(NA, Smax),
      distance = rep(0, Smax)
    )
    dist_anc$ancestor[1] = NA
    dist_anc$A_E[1] = "A"
    
    # To record dist_anc across steps
    list_dist_anc <- list()
    
    # To record interaction matrices
    L_list = list()
    L_cropped_list <- list()
    
    # Guild dataframe
    df_guild <- data.frame(step = integer(), npredators = integer(), nherbivores = integer(), nomnivores = integer())
    
    # Number of species currently existing
    S = 1
    
    ##########################
    # MAIN LOOP
    for(step in 2:nsteps) {
      ActualS = sum(pres[step - 1,])
      
      for(i in 1:Smax) {
        if(S >= Smax) break
        
        if(pres[step - 1, i] == 1) {
          dist_anc$A_E[i] <- "A"
          dist_anc$distance[i] <- dist_anc$distance[i] + 1
          pres[step, i] = 1
          
          test_number = runif(1)
          speciation_prob = u_max / (1 + exp(d * (ActualS - I_max)))
          
          if(test_number < speciation_prob) {
            traits_anc <- traits_mat[i, ]
            first_mutant_direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
            traits_mut1 <- rand_traits_mut(traits_anc, pars, direction = first_mutant_direction)
            
            I = get_L_vec(basal, pars, traits_mat, traits_mut1)
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))
            
            if(sum_I > 0) {
              estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
              estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
            } else {
              estab_prob = 0
            }
            
            if(runif(1) < estab_prob) {
              S = S + 1
              traits_mat[S, ] = traits_mut1
              pres[step, S] = 1
              anc[S,] = c(step, species_names[i], species_names[S])
              dist_anc$spp[S] = species_names[S]
              dist_anc$ancestor[S] = species_names[i]
              dist_anc$A_E[S] = "A"
              dist_anc$distance[S] = 1
              
              if(abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
                pres[step, i] = 0
                dist_anc$A_E[i] = "E"
              }
              
              if(S >= Smax) break
            }
            
            second_mutant_direction <- ifelse(first_mutant_direction == "greater", "lesser", "greater")
            traits_mut2 <- rand_traits_mut(traits_anc, pars, direction = second_mutant_direction)
            
            I = get_L_vec(basal, pars, traits_mat, traits_mut2)
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))
            
            if(sum_I > 0) {
              estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
              estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
            } else {
              estab_prob = 0
            }
            
            if(runif(1) < estab_prob) {
              S = S + 1
              traits_mat[S, ] = traits_mut2
              pres[step, S] = 1
              anc[S,] = c(step, species_names[i], species_names[S])
              dist_anc$spp[S] = species_names[S]
              dist_anc$ancestor[S] = species_names[i]
              dist_anc$A_E[S] = "A"
              dist_anc$distance[S] = 1
              
              if(abs(traits_anc["n"] - traits_mut1["n"]) <= 0.05) {
                pres[step, i] = 0
                dist_anc$A_E[i] = "E"
              }
              
              if(S >= Smax) break
            }
          }
        }
      }
      
      if(S >= Smax) break
      
      # EXTINCTION DYNAMICS
      pres_vec = pres[step, ]
      cooc = outer(pres_vec, pres_vec, FUN = "*")
      L = get_L_mat(basal, pars, traits_mat)
      L[(Sbasal + 1):(Sbasal + Smax), ] = L[(Sbasal + 1):(Sbasal + Smax), ] * cooc
      rownames(L) <- c(paste0("Basal", 1:Sbasal), species_names)
      colnames(L) <- species_names
      L_list[[step]] = L
      
      present_non_basal = which(pres[step, ] == 1)
      if(length(present_non_basal) > 1) {
        L_cropped = L[c(1:Sbasal, Sbasal + present_non_basal), present_non_basal]
        rownames(L_cropped) <- c(paste0("Basal", 1:Sbasal), species_names[present_non_basal])
        colnames(L_cropped) <- species_names[present_non_basal]
        L_cropped_list[[step]] <- L_cropped
      }
      
      # EXTINCTION PROBABILITY
      # [... keep extinction calculation logic ...]
      
      # Update extinct species
      for(i in 1:Smax) {
        if(!is.na(pres[step, i]) && pres[step, i] != pres[step-1, i] && pres[step, i] == 0) {
          dist_anc$A_E[i] = "E"
        }
      }
      
      # Save dist_anc at this step
      list_dist_anc[[step]] = na.omit(dist_anc)
      
      # EARLY STOP CRITERIA
      if(step == 10 && sum(pres[step, ]) == 0) return(NULL)
      if(step > 10 && sum(pres[step, ]) == 0) return(NULL)
      if(step == 50 && sum(pres[step, ]) < 10) return(NULL)
      if(step == 30 && sum(pres[step, ]) <= 3) return(NULL)
      if(step == nsteps && sum(pres[step, ]) < 40) return(NULL)
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


