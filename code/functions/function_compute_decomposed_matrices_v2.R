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
  
  # number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  #### Identify timesteps where phylogenetic distances cant be calculated
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  # until what timestep need to discard:
  final.discarded_timestep <- (non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)])+1
  
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
    Sbasals <- nbasals
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
  
  if(nrow(presence_matrix) != length(list_networks_sppnames_letters)){
    print("PROBLEM - presence matrix, list phylo dist and list interaction networks dont have the same n_steps")
  }
  
  ## -------------------
  
  ## Loop for obtaining phylogenetic distances:
  list_svd_eigen.phy <- list()
  list_phylo.corr_cropped <- list()
  
  
  for (i in 1:length(list_anc_dist_letters)) {
    
    print(paste("step", i))
    
    # obtain phylogenetic distances between all tree nodes
    
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";",newick_tail))
    #list_dist.phylo[[i]] <- cophenetic.phylo(tree)
    
    tree$edge.length<-sapply(tree$edge.length,function(x) ifelse(x==0,1e-5,x))
    phylo.vcv<-vcv(tree)
    phylo.corr<-cov2cor(phylo.vcv)
    
    # Crop matrix of phylogenetic distances with those species that are alive
    
    # alive_species <- list_anc_dist_letters[[i]][which(list_anc_dist_letters[[i]]$`A/E` == "A"), "spp"]
    alive_species <- names(which(presence_matrix[i, ] == 1))
    list_phylo.corr_cropped[[i]] <- phylo.corr[alive_species, alive_species]
    
    # Check that species in the phylogenetic distance matrix are the same present in the presence matrix
    
    if (all(names(which(presence_matrix[i, ] == 1)) != colnames(list_phylo.corr_cropped[[i]]))) {
      cat("Species mismatch detected for index", i, "\n")
    }
    
    
    list_svd_eigen.phy[[i]] <-eigen(list_phylo.corr_cropped[[i]], symmetric = T)$vec
  }
  
  
  ## ------------------- Check
  if (length(which(unlist(lapply(list_svd_eigen.phy, is.null)) == TRUE)) != 0) {
    print("PROBLEM - nulls in list_dist.phylo")
  }
  if (length(list_svd_eigen.phy) != length(list_networks_sppnames_letters)){
    print("PROBLEM - length(list_dist.phylo) != length(list_networks_sppnames_letters)")
  }
  ## -------------------
  
  
  
  ## ------------------- Check
  #Check that phylogenetic distance matrices retain present species:
  vec_error <- c()
  for (i in 1:length(list_svd_eigen.phy)) {
    vec_tf <-  names(presence_matrix[i, which(presence_matrix[i,] == 1)]) == colnames(list_svd_eigen.phy[[i]])
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
  
  
  
  
  
  #############################################################################
  
  
  list_svd_pred <- list()
  
  for (i in 1:length(list_net_present_spp.letters)) {
    
    # Perform SVD
    svd_result <- svd(list_net_present_spp.letters[[i]])
    
    # Extract left (U), middle (D), and right (V) matrices
    U <- svd_result$u       # Left matrix (n_prey x kept axes)
    D <- diag(svd_result$d) # Middle matrix (kept axes x kept axes)
    V <- svd_result$v       # Right matrix (n_predators x kept axes)
    
    # Decide on the number of axes to keep
    kept_axes <- ncol(list_net_present_spp.letters[[i]])  
    
    # Select the kept axes from U, D, and V
    U_kept <- U[, 1:kept_axes]   # n_prey x kept axes
    D_kept <- D[1:kept_axes, 1:kept_axes] # kept axes x kept axes
    V_kept <- V[, 1:kept_axes]   # n_predators x kept axes
    
    # Transpose the right matrix V to match the desired output
    list_svd_pred[[i]] <- t(V_kept)  
    
  }
  
  
  
  ## ------------------- Check
  if(length(list_svd_pred) != length(list_svd_eigen.phy)){
    print("PROBLEM - length(list_interact_distances_mean_corrected) != length(list_dist.phylo_pres")
  }
  if(length(list_svd_eigen.phy) != length(list_svd_pred)) {
    print("PROBLEM - list phylo dist and lists interact dist dont have the same length")
  } 
  
  vec_problems_ncol <- c()
  
  for (i in 1:length(list_svd_pred)) {
    if(ncol(list_svd_pred[[i]]) != ncol(list_svd_eigen.phy[[i]])){
      vec_problems_ncol[i] <- "P"
    } else if(ncol(list_svd_pred[[i]]) == ncol(list_svd_eigen.phy[[i]])){
      vec_problems_ncol[i] <- "g"
    }
  }
  if(length(which(vec_problems_ncol == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same ncols")
  }
  
  vec_problems_sppcomp <- c()
  
  for (i in 1:length(list_svd_eigen.phy)) {
    vec_problems_sppcomp[i] <- identical(sort(colnames(list_svd_pred[[i]])), sort(colnames(list_svd_eigen.phy[[i]])))
  }
  if(length(which(vec_problems_ncol == "FALSE") > 0)){
    print("PROBLEM -  Interact and phylo dist. matrices dont have the order of colnames")
  }
  
  vec_problems <- c()
  
  for (i in 1:length(list_svd_pred)) {
    vec_truefalse <- colnames(list_svd_pred[[i]]) == colnames(list_svd_eigen.phy[[i]])
    if(FALSE %in% vec_truefalse){
      vec_problems[i] <- "P"
    }else{
      vec_problems[i] <- "_"
    }
  }
  if(length(which(vec_problems == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same order of colnames")
  }
  
  
  
  result <- list("list_svd_pred" = list_svd_pred,
                 "list_svd_eigen.phy" = list_svd_eigen.phy,
                 "list_net_present_spp.letters" = list_net_present_spp.letters,
                 "list_phylo.corr_cropped" = list_phylo.corr_cropped)
  
  return(result)
}



####################### UPDATE

compute_result_matrices <- function(results_simulation, int, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  #### Identify timesteps where phylogenetic distances cant be calculated
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  # until what timestep need to discard:
  final.discarded_timestep <- (non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)])+1
  
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
    Sbasals <- nbasals
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
  
  if(nrow(presence_matrix) != length(list_networks_sppnames_letters)){
    print("PROBLEM - presence matrix, list phylo dist and list interaction networks dont have the same n_steps")
  }
  
  ## -------------------
  
  ## Loop for obtaining phylogenetic distances:
  list_svd_eigen.phy <- list()
  list_phylo.corr_cropped <- list()
  
  
  list_svd_eigen.phy <- list()
  list_phylo.corr_cropped <- list()
  list_net_present_spp.letters <- list()
  list_svd_pred <- list()
  valid_timesteps <- c()
  
  for (i in 1:length(list_anc_dist_letters)) {
    cat("step", i, "\n")
    
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root", ";", newick_tail))
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    alive_species <- names(which(presence_matrix[i, ] == 1))
    valid_species <- intersect(alive_species, rownames(phylo.corr))
    
    if (length(valid_species) < 3) {
      warning(paste("⚠️ Skipping step", i, "- too few valid species"))
      next
    }
    
    valid_timesteps <- c(valid_timesteps, i)
    
    # Phylogenetic SVD
    phy_corr_crop <- phylo.corr[valid_species, valid_species]
    list_phylo.corr_cropped[[length(valid_timesteps)]] <- phy_corr_crop
    list_svd_eigen.phy[[length(valid_timesteps)]] <- eigen(phy_corr_crop, symmetric = TRUE)$vec
    
    # Network cropping and SVD
    net <- list_networks_sppnames_letters[[i]]
    net_crop <- net[valid_species, valid_species]
    list_net_present_spp.letters[[length(valid_timesteps)]] <- net_crop
    
    svd_net <- svd(net_crop)
    kept_axes <- ncol(net_crop)
    list_svd_pred[[length(valid_timesteps)]] <- t(svd_net$v[, 1:kept_axes])
  }
  
  # Filter original lists to match valid steps only
  list_networks_sppnames_letters <- list_networks_sppnames_letters[valid_timesteps]
  list_anc_dist_letters <- list_anc_dist_letters[valid_timesteps]
  
  # Update the presence matrix to keep only valid timesteps
  presence_matrix <- presence_matrix[valid_timesteps, , drop = FALSE]
  
  
  
  ## ------------------- Check
  if (length(which(unlist(lapply(list_svd_eigen.phy, is.null)) == TRUE)) != 0) {
    print("PROBLEM - nulls in list_dist.phylo")
  }
  if (length(list_svd_eigen.phy) != length(list_networks_sppnames_letters)){
    print("PROBLEM - length(list_dist.phylo) != length(list_networks_sppnames_letters)")
  }
  ## -------------------
  
  
  
  ## ------------------- Check
  #Check that phylogenetic distance matrices retain present species:
  vec_error <- c()
  for (i in 1:length(list_svd_eigen.phy)) {
    vec_tf <-  names(presence_matrix[i, which(presence_matrix[i,] == 1)]) == colnames(list_svd_eigen.phy[[i]])
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
  
  
  
  
  
  #############################################################################
  
  
  list_svd_pred <- list()
  
  for (i in 1:length(list_net_present_spp.letters)) {
    
    # Perform SVD
    svd_result <- svd(list_net_present_spp.letters[[i]])
    
    # Extract left (U), middle (D), and right (V) matrices
    U <- svd_result$u       # Left matrix (n_prey x kept axes)
    D <- diag(svd_result$d) # Middle matrix (kept axes x kept axes)
    V <- svd_result$v       # Right matrix (n_predators x kept axes)
    
    # Decide on the number of axes to keep
    kept_axes <- ncol(list_net_present_spp.letters[[i]])  
    
    # Select the kept axes from U, D, and V
    U_kept <- U[, 1:kept_axes]   # n_prey x kept axes
    D_kept <- D[1:kept_axes, 1:kept_axes] # kept axes x kept axes
    V_kept <- V[, 1:kept_axes]   # n_predators x kept axes
    
    # Transpose the right matrix V to match the desired output
    list_svd_pred[[i]] <- t(V_kept)  
    
  }
  
  
  
  ## ------------------- Check
  if(length(list_svd_pred) != length(list_svd_eigen.phy)){
    print("PROBLEM - length(list_interact_distances_mean_corrected) != length(list_dist.phylo_pres")
  }
  if(length(list_svd_eigen.phy) != length(list_svd_pred)) {
    print("PROBLEM - list phylo dist and lists interact dist dont have the same length")
  } 
  
  vec_problems_ncol <- c()
  
  for (i in 1:length(list_svd_pred)) {
    if(ncol(list_svd_pred[[i]]) != ncol(list_svd_eigen.phy[[i]])){
      vec_problems_ncol[i] <- "P"
    } else if(ncol(list_svd_pred[[i]]) == ncol(list_svd_eigen.phy[[i]])){
      vec_problems_ncol[i] <- "g"
    }
  }
  if(length(which(vec_problems_ncol == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same ncols")
  }
  
  vec_problems_sppcomp <- c()
  
  for (i in 1:length(list_svd_eigen.phy)) {
    vec_problems_sppcomp[i] <- identical(sort(colnames(list_svd_pred[[i]])), sort(colnames(list_svd_eigen.phy[[i]])))
  }
  if(length(which(vec_problems_ncol == "FALSE") > 0)){
    print("PROBLEM -  Interact and phylo dist. matrices dont have the order of colnames")
  }
  
  vec_problems <- c()
  
  for (i in 1:length(list_svd_pred)) {
    vec_truefalse <- colnames(list_svd_pred[[i]]) == colnames(list_svd_eigen.phy[[i]])
    if(FALSE %in% vec_truefalse){
      vec_problems[i] <- "P"
    }else{
      vec_problems[i] <- "_"
    }
  }
  if(length(which(vec_problems == "P") > 0)){
    print("PROBLEM - Interact and phylo dist. matrices dont have the same order of colnames")
  }
  
  
  
  
  
  ##########################################################
  # Create a table of species lifespan and descendants
all_anc_dist <- results_simulation$list_anc_dist

# Vector of all species (original names as characters)
all_species <- unique(unlist(lapply(all_anc_dist, function(x) x$spp)))

# Initialize data.frame
species_life_summary <- data.frame(
  spp_number = all_species,
  spp_letter = NA,
  timestep_birth = NA,
  timestep_extinct = NA,
  n_descendants = 0,
  stringsAsFactors = FALSE
)

# Convert column names of presence_matrix to letters (in sync with earlier conversion)
spp_letters <- colnames(presence_matrix)

# Assign letter names
species_life_summary$spp_letter <- spp_letters[match(species_life_summary$spp_number, chartr("ABCDEFGHIJ", "0123456789", spp_letters))]

# Fill birth and extinction times
for (i in 1:nrow(species_life_summary)) {
  spp_letter <- species_life_summary$spp_letter[i]
  
  if (spp_letter %in% spp_letters) {
    presence_vector <- presence_matrix[, spp_letters == spp_letter]
    
    if (any(presence_vector == 1)) {
      species_life_summary$timestep_birth[i] <- which(presence_vector == 1)[1]
      species_life_summary$timestep_extinct[i] <- tail(which(presence_vector == 1), 1)
    }
  }
}

# Count descendants using the full untrimmed ancestry table
# Count unique descendant species per ancestor
# Gather all ancestor–descendant relationships across timesteps
all_relationships <- do.call(rbind, lapply(all_anc_dist, function(x) {
  x[, c("spp", "ancestor")]
}))

# Remove rows where ancestor is 0 (i.e., origin)
all_relationships <- all_relationships[all_relationships$ancestor != "0", ]

# For each species, count how many **unique spp** list it as their ancestor
unique_descendants <- split(all_relationships$spp, all_relationships$ancestor)
desc_count <- sapply(unique_descendants, function(x) length(unique(x)))

# Match to the summary table
species_life_summary$n_descendants <- ifelse(
  species_life_summary$spp_number %in% names(desc_count),
  desc_count[species_life_summary$spp_number],
  0
)


  
result <- list("list_svd_pred" = list_svd_pred,
               "list_svd_eigen.phy" = list_svd_eigen.phy,
               "list_net_present_spp.letters" = list_net_present_spp.letters,
               "list_phylo.corr_cropped" = list_phylo.corr_cropped,
               "species_life_summary" = species_life_summary,
               "presence_matrix" = presence_matrix)  # ← add this line
  
  return(result)
}