######## compute phylogenetic signal from simulation results


check_matrix.with.values <- function(matrix) {
  result <- any(matrix != 0)
  return(result)
}


eliminate_basals <- function(matrix, nbasals) {
  
  mat <- matrix[(nbasals+1):nrow(matrix),]
  
  return(mat)
}


compute_df_phylosignal_fromsimulation <- function(results_simulation, int, Smax, sim, input_matrix_procrustes) {
  
  
  
  presence_matrix <- results_simulation$presence_matrix
  
  # number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  #### Crop presence matrix to n_steps
  presence_matrix <- presence_matrix[1:n_steps,]
  
  #### Identify timesteps where phylogenetic distances cant be calculated
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  # until what timestep need to discard:
  final.discarded_timestep <- non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]
  
  #### homogenize elements to start from valid timesteps
  # ancestry table
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep+1):length(results_simulation$list_anc_dist)]
  # network list
  network_list <- results_simulation$network_list[(final.discarded_timestep+1):length(results_simulation$list_anc_dist)]
  
  
  ## ------------------- Check
  if (length(which(is.null(results_simulation$network_list))) > 0) {
    print("PROBLEM - null network somewhere")
  } 
  if(length(list_anc_dist) != length(network_list)){
    print("PROBLEM - length list_anc_dist != length(network_list)")
  }
  ## -------------------
  
  
  
  if(int == "foodweb"){
    #### Eliminate basal species
    Sbasals <- length(results_simulation$basal)
    network_list <- lapply(network_list, eliminate_basals, nbasals = Sbasals)
  }
  
  #### Convert spp names from numbers to letters
  ## ancestry-distances table
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  ## Network list 
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  
  
  ## ------------------- Check
  if(length(list_anc_dist_letters) != length(list_networks_sppnames_letters)){
    print("PROBLEM - list_anc_dist_letters != list_networks_sppnames_letters")
  } 
  ## -------------------
  
  
  ## Loop for obtaining phylogenetic distances:
  list_newick <- list()
  list_trees <- list()
  list_newick_tails <- list()
  list_dist.phylo <- list()
  
  for (i in 1:length(list_anc_dist_letters)) {
    list_newick[[i]] <- ToPhylo(list_anc_dist_letters[[i]])
    list_newick_tails[[i]] <- paste(list_newick[[i]], "root")
    list_trees[[i]] <- read.tree(text = sub("A root",";",list_newick_tails[[i]]))
    list_dist.phylo[[i]] <- cophenetic.phylo(list_trees[[i]])
  }
  
  
  ## ------------------- Check
  if (length(which(unlist(lapply(list_dist.phylo, is.null)) == TRUE)) != 0) {
    print("PROBLEM - nulls in list_dist.phylo")
  }
  if (length(list_dist.phylo) != length(list_networks_sppnames_letters)){
    print("PROBLEM - length(list_dist.phylo) != length(list_networks_sppnames_letters)")
  }
  ## -------------------
  
  
  #### convert spp names to letters in presence matrix
  colnames(presence_matrix) <- seq(1:Smax)
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  
  
  ## ------------------- Check
  
  if(length(which(colnames(presence_matrix) != colnames(list_networks_sppnames_letters[[1]]))) != 0){
    print("PROBLEM - spp names in presence_matrix dont correspond to spp names in list_networks_sppnames_letterst")
  }
  ## -------------------
  
  
  
  # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
  presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(results_simulation$list_anc_dist),]
  
  
  
  ## ------------------- Check
  
  if(nrow(presence_matrix) != length(list_networks_sppnames_letters) | 
     nrow(presence_matrix) != length(list_dist.phylo)){
    print("PROBLEM - presence matrix, list phylo dist and list interaction networks dont have the same n_steps")
  }
  ## -------------------
  
  
  #### Retain only present species in phylogenetic distance matrices
  list_dist.phylo_pres <- list()
  
  for (i in 1:length(list_dist.phylo)) {
    list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
  }
  
  
  ## ------------------- Check
  #Check that phylogenetic distance matrices retain present species:
  vec_error <- c()
  for (i in 1:nrow(presence_matrix)) {
    vec_tf <-  names(presence_matrix[i, which(presence_matrix[i,] == 1)]) == colnames(list_dist.phylo_pres[[i]])
    if(length(which(vec_tf == FALSE)) > 0){
      vec_error[i] <- "error"
    } else if (length(which(vec_tf == FALSE)) == 0 ){
      vec_error[i] <- "g"
    }
  }
  if(length(which(vec_error == "error")) > 0){
    print("PROBLEM - list_dist.phylo_pres dont retain present species")
  }
  ## -------------------
  
  
  #### Retain only present species in network matrices
  list_net_present_spp.letters <- list()
  
  for (i in 1:length(list_networks_sppnames_letters)) {
    list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
  }
  
  
  ## Compute interaction distances (NMI)
  list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_pred)
  list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_prey)
  # set diagonal to 0
  list_interact_distances_pred <- lapply(list_interact_distances_pred, FUN = diag_to0)
  list_interact_distances_prey <- lapply(list_interact_distances_prey, FUN = diag_to0)
  # Assign 0 to NaN
  list_interact_distances_pred_corrected <- list()
  list_interact_distances_prey_corrected <- list()
  
  for (i in 1:length(list_interact_distances_pred)) {
    list_interact_distances_pred_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_pred[[i]], marg = 2) 
  }
  
  for (i in 1:length(list_interact_distances_prey)) {
    list_interact_distances_prey_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_prey[[i]], marg = 1)
  }
  
  
  ## ------------------- Check
  if(length(list_interact_distances_prey_corrected) != length(list_interact_distances_pred_corrected) |
     length(list_interact_distances_prey_corrected) != length(list_dist.phylo_pres)){
    
    print("PROBLEM - lists of phylo dist and interactions dont the same length")
    
  } 
  ## -------------------
  
  
  
  # Compute mean distance (as preys and predators)
  list_interact_distances_mean_corrected <- list()
  
  for (i in 1:length(list_interact_distances_pred_corrected)) {
    list_interact_distances_mean_corrected[[i]] <- (list_interact_distances_pred_corrected[[i]] + list_interact_distances_prey_corrected[[i]]) / 2
  }
  
  
  
  ## ------------------- Check
  if(length(list_interact_distances_mean_corrected) != length(list_dist.phylo_pres)){
    print("PROBLEM - length(list_interact_distances_mean_corrected) != length(list_dist.phylo_pres")
  }
  if(length(list_dist.phylo) != length(list_interact_distances_pred) |
     length(list_dist.phylo) != length(list_interact_distances_prey)) {
    print("PROBLEM - list phylo dist and lists interact dist dont have the same length")
  } 
  
  vec_problems_ncol <- c()
  
  for (i in 1:length(list_interact_distances_mean_corrected)) {
    if(ncol(list_interact_distances_mean_corrected[[i]]) != ncol(list_dist.phylo_pres[[i]])){
      vec_problems_ncol[i] <- "P"
    } else if(ncol(list_interact_distances_mean_corrected[[i]]) == ncol(list_dist.phylo_pres[[i]])){
      vec_problems_ncol[i] <- "g"
    }
  }
  if(length(which(vec_problems_ncol == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same ncols")
  }
  
  vec_problems_sppcomp <- c()
  
  for (i in 1:length(list_dist.phylo_pres)) {
    vec_problems_sppcomp[i] <- identical(sort(colnames(list_interact_distances_mean_corrected[[i]])), sort(colnames(list_dist.phylo_pres[[i]])))
  }
  if(length(which(vec_problems_ncol == "FALSE") > 0)){
    print("PROBLEM -  Interact and phylo dist. matrices dont have the order of colnames")
  }
  
  vec_problems <- c()
  
  for (i in 1:length(list_interact_distances_mean_corrected)) {
    vec_truefalse <- colnames(list_interact_distances_mean_corrected[[i]]) == colnames(list_dist.phylo_pres[[i]])
    if(FALSE %in% vec_truefalse){
      vec_problems[i] <- "P"
    }else{
      vec_problems[i] <- "_"
    }
  }
  if(length(which(vec_problems == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same order of colnames")
  }
  
  ## -------------------
  
  
  
  ## make sure we have matrices
  list_dist_dist.phylo_pres <- list()
  list_dist_interact_distances_pred_corrected <- list()
  list_dist_interact_distances_prey_corrected <- list()
  list_dist_interact_distances_mean_corrected <- list()
  
  for (i in 1:length(list_dist.phylo_pres)) {
    list_dist_dist.phylo_pres[[i]] <- as.matrix(list_dist.phylo_pres[[i]])
    list_dist_interact_distances_pred_corrected[[i]] <- as.matrix(list_interact_distances_pred_corrected[[i]])
    list_dist_interact_distances_prey_corrected[[i]] <- as.matrix(list_interact_distances_prey_corrected[[i]])
    list_dist_interact_distances_mean_corrected[[i]] <- as.matrix(list_interact_distances_mean_corrected[[i]])
  }
  
  
  
  
  # identify what matrices have all interaction distances = 0 and discard them
  
  vec_timesteps_all0 <- which(lapply(list_interact_distances_pred_corrected,check_matrix.with.values) != TRUE)
  
  if(length(vec_timesteps_all0) > 0){
    list_dist_dist.phylo_pres <- list_dist_dist.phylo_pres[-c(vec_timesteps_all0)]
    list_dist_interact_distances_pred_corrected <- list_dist_interact_distances_pred_corrected[-c(vec_timesteps_all0)]
    list_dist_interact_distances_prey_corrected <- list_dist_interact_distances_prey_corrected[-c(vec_timesteps_all0)]
    list_dist_interact_distances_mean_corrected <- list_dist_interact_distances_mean_corrected[-c(vec_timesteps_all0)]
    
  }
  
  
  
  
  
  
  if(input_matrix_procrustes == "pco"){
    
    
    ## compute Principal Coordinate Analyses
    list_pco.phy <- list()
    list_pco.int_pred <- list()
    list_pco.int_prey <- list()
    list_pco.int_mean <- list()
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      pco_phy <- cmdscale(list_dist_dist.phylo_pres[[i]], k = 2)
      list_pco.phy[[i]] <- pco_phy
      pco_pred <- cmdscale(list_dist_interact_distances_pred_corrected[[i]], k = 2)
      list_pco.int_pred[[i]] <- pco_pred
      pco_prey <- cmdscale(list_dist_interact_distances_prey_corrected[[i]], k = 2)
      list_pco.int_prey[[i]] <- pco_prey
      pco_mean <- cmdscale(list_dist_interact_distances_mean_corrected[[i]], k = 2)
      list_pco.int_mean[[i]] <- pco_mean
    }
    
    
    ## Compute correlation - Procrustes test
    
    protest_pred <- list()
    procrustes_pred <- list()
    
    protest_prey <- list()
    procrustes_prey <- list()
    
    protest_mean <- list()
    procrustes_mean <- list()
    
    
    for (i in 1:length(list_pco.phy)) {
      protest_pred[[i]] <- protest(list_pco.int_pred[[i]], list_pco.phy[[i]])
      procrustes_pred[[i]] <- procrustes(list_pco.int_pred[[i]],list_pco.phy[[i]])
    }
    
    for (i in 1:length(list_pco.phy)) {
      protest_prey[[i]] <- protest(list_pco.int_prey[[i]], list_pco.phy[[i]])
      procrustes_prey[[i]] <- procrustes(list_pco.int_prey[[i]],list_pco.phy[[i]])
    }
    
    for (i in 1:length(list_pco.phy)) {
      protest_mean[[i]] <- protest(list_pco.int_mean[[i]], list_pco.phy[[i]])
      procrustes_mean[[i]] <- procrustes(list_pco.int_mean[[i]],list_pco.phy[[i]])
    }
    
  } else if(input_matrix_procrustes == "matrix"){
    
    ## Compute correlation - Procrustes test
    
    protest_pred <- list()
    procrustes_pred <- list()
    
    protest_prey <- list()
    procrustes_prey <- list()
    
    protest_mean <- list()
    procrustes_mean <- list()
    
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      protest_pred[[i]] <- protest(list_dist_interact_distances_pred_corrected[[i]], list_dist_dist.phylo_pres[[i]])
      procrustes_pred[[i]] <- procrustes(list_dist_interact_distances_pred_corrected[[i]],list_dist_dist.phylo_pres[[i]])
    }
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      protest_prey[[i]] <- protest(list_dist_interact_distances_prey_corrected[[i]], list_dist_dist.phylo_pres[[i]])
      procrustes_prey[[i]] <- procrustes(list_dist_interact_distances_prey_corrected[[i]],list_dist_dist.phylo_pres[[i]])
    }
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      protest_mean[[i]] <- protest(list_dist_interact_distances_mean_corrected[[i]], list_dist_dist.phylo_pres[[i]])
      procrustes_mean[[i]] <- procrustes(list_dist_interact_distances_mean_corrected[[i]],list_dist_dist.phylo_pres[[i]])
    }
    
    
    
  }
  
  

  
  
  #### Create dataframe results
  protest_pval_pred <- c()
  protest_corr_pred <- c()
  
  protest_pval_prey <- c()
  protest_corr_prey <- c()
  
  protest_pval_mean <- c()
  protest_corr_mean <- c()
  
  for (i in 1:length(protest_mean)) {
    
    protest_pval_pred[i] <- protest_pred[[i]]$signif
    protest_corr_pred[i] <-protest_pred[[i]]$t0
    
    protest_pval_prey[i] <- protest_prey[[i]]$signif
    protest_corr_prey[i] <-protest_prey[[i]]$t0
    
    protest_pval_mean[i] <- protest_mean[[i]]$signif
    protest_corr_mean[i] <-protest_mean[[i]]$t0
  }
  
  
  ## add timesteps
  if(length(vec_timesteps_all0) > 0){
    timesteps <- (final.discarded_timestep+1):n_steps
    timesteps <- timesteps[-c(vec_timesteps_all0)]
  }else{
    timesteps <- (final.discarded_timestep+1):n_steps
  }
  
  df_phylosign_time_sim <- data.frame("timesteps" = timesteps,
                                   "phylosign_cor_mean" = protest_corr_mean,
                                   "phylosign_p_mean" = protest_pval_mean,
                                   "phylosign_cor_pred" = protest_corr_pred,
                                   "phylosign_p_pred" = protest_pval_pred,
                                   "phylosign_cor_prey" = protest_corr_prey,
                                   "phylosign_p_prey" = protest_pval_prey,
                                   "sim" = rep(sim, times = length(timesteps)),
                                   "interaction" = rep(int, times = length(timesteps)))
  
  ## add S
  n_spp_time <- c()
  
  for (i in 1:length(list_dist_dist.phylo_pres)) {
    n_spp_time[i] <- nrow(list_dist_dist.phylo_pres[[i]])
  }
  
  ## add list of distance matrices
  df_phylosign_time_sim$nspp <- n_spp_time
  
  df_phylosign_time_sim$list_phylo_dist <- list_dist_dist.phylo_pres
  df_phylosign_time_sim$list_dist_interact_distances_mean_corrected <- list_dist_interact_distances_mean_corrected
  df_phylosign_time_sim$list_dist_interact_distances_pred_corrected <- list_dist_interact_distances_pred_corrected
  df_phylosign_time_sim$list_dist_interact_distances_prey_corrected <- list_dist_interact_distances_prey_corrected
  
  
  return(df_phylosign_time_sim)
  
}
