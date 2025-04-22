
#### FUNCTIONS FOR MEASURING NETWORK PROPERTIES



# Get the adjacency matrix from simulation objects (.rds) stored in a folder (file_path)


get_adjancency_from_simulation.objects <- function(file_path) {
  
  
  adj_matrix_list <- list()
  
  for (i in 1:length(file_path)) {
    
    list_files <- readRDS(file_path[[i]])
    
    ## spp names from numbers to letters
    
    net_numbers <- set_sppNames_numbers(list_files$network_list[[length(list_files$network_list)]])
    
    net_letters <- convert_sppnames_toletters(net_numbers)
    
    
    ## Crop network for present spp
    
    
    colnames(list_files$presence_matrix) <- seq(1:1000)
    
    colnames(list_files$presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(list_files$presence_matrix))
    
    list_files$presence_matrix <- list_files$presence_matrix[length(list_files$network_list),]
    
    adj_matrix_list[[i]] <- net_letters[names(which(list_files$presence_matrix == 1)), names(which(list_files$presence_matrix == 1))]
    
    
    
  }
  
  return(adj_matrix_list)
  
} 





get_adjancency_from_simulation.objects_fw <- function(file_path) {
  
  
  adj_matrix_list <- list()
  
  
  
  for (i in 1:length(file_path)) {
    
    list_files <- readRDS(file_path[[i]])
    
    mat_f <- list_files$network_list[[length(list_files$network_list)]]
    
    nbasals <- length(list_files$basal)
    
    # create the columns for basals (all 0)
    
    mat_basals_colums <- matrix(0, nrow = nrow(mat_f), ncol = nbasals)
    
    # Add basal spp as columns to the interaction matrix
    
    mat <- cbind(mat_basals_colums, mat_f)
    
    ## spp names from numbers to letters
    
    net_numbers <- set_sppNames_numbers(mat)
    
    net_letters <- convert_sppnames_toletters(net_numbers)
    
    
    
    # add basals to presence matrix (all 1)
    
    pres_mat <- list_files$presence_matrix
    
    mat_basals_colums_pres <- matrix(1, nrow = nrow(pres_mat), ncol =  nbasals)
    
    mat_press <- cbind(mat_basals_colums_pres, pres_mat)
    
    colnames(mat_press) <- seq(1:ncol(mat_press))
    
    colnames(mat_press) <- chartr("0123456789", "ABCDEFGHIJ", colnames(mat_press))
    
    
    
    # crop interaction network with present species at last timestep
    
    mat_press <- mat_press[length(list_files$network_list),]
    
    
    adj_matrix_list[[i]] <- net_letters[names(which(mat_press == 1)),
                                names(which(mat_press == 1))]
    
    
    
  }
  
  return(adj_matrix_list)
  
} 



get_adjancency_from_simulation.fw_single <- function(list_networks, nbasals, pres_mat) {
  

    
    mat_f <- list_networks[[length(list_networks)]]

    # create the columns for basals (all 0)
    
    mat_basals_colums <- matrix(0, nrow = nrow(mat_f), ncol = nbasals)
    
    # Add basal spp as columns to the interaction matrix
    
    mat <- cbind(mat_basals_colums, mat_f)
    
    ## spp names from numbers to letters
    
    net_numbers <- set_sppNames_numbers(mat)
    
    net_letters <- convert_sppnames_toletters(net_numbers)
    
    
    
    # add basals to presence matrix (all 1)
    
    mat_basals_colums_pres <- matrix(1, nrow = nrow(pres_mat), ncol =  nbasals)
    
    mat_press <- cbind(mat_basals_colums_pres, pres_mat)
    
    colnames(mat_press) <- seq(1:ncol(mat_press))
    
    colnames(mat_press) <- chartr("0123456789", "ABCDEFGHIJ", colnames(mat_press))
    
    
    
    # crop interaction network with present species at last timestep
    
    mat_press <- mat_press[length(list_networks),]
    
    
    adjacency_f <- net_letters[names(which(mat_press == 1)),
                                        names(which(mat_press == 1))]
    
    

  
  return(adjacency_f)
  
} 






get_network_measures_from_adjacency <- function(list_adj_matrix) {
  
  columns <- c("S",
               "Link_density",
               "C",
               "perc_tops",
               "perc_int",
               "perc_basals",
               "perc_cannibals",
               "omnivory",
               "mean_gen",
               "sd_gen",
               "mean_vul",
               "sd_vul")
  
  df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df) = columns
  
  for (i in 1:length(list_adj_matrix)) {
    
    graph <- graph_from_adjacency_matrix(adjmatrix = list_adj_matrix[[i]],
                                         mode = "directed")
    
    # 1. Number of species within the foodweb (S)
    S <- gorder(graph)  # Number of vertices in the graph
    df[i, "S"] <- S
    
    # Skip computation for metrics if S <= 1, as they are undefined or not meaningful
    if (S <= 1) {
      df[i, c("Link_density", "C", "perc_tops", "perc_int", 
              "perc_basals", "perc_cannibals", "omnivory", 
              "mean_gen", "sd_gen", "mean_vul", "sd_vul")] <- NA
      next
    }
    
    # 2. Links density (L/S)
    L <- ecount(graph)  # Number of edges
    df[i, "Link_density"] <- L / S
    
    # 3. Connectance (C = L / S^2)
    df[i, "C"] <- L / S^2
    
    # 4. %T (Top species that have resource species but lack any consumer species)
    # Top predators: Species with in-degree = 0 (they consume others) and out-degree > 0 (no one consumes them)
    n_tops <- length(V(graph)[igraph::degree(graph, mode = 'in') == 0 & igraph::degree(graph, mode = 'out') > 0])
    df[i, "perc_tops"] <- (n_tops / S) * 100
    
    # 5. % I (Intermediate species that have both resource and consumer species)
    # Intermediate species: in-degree > 0 (consumed by others) and out-degree > 0 (they also consume)
    n_int <- length(V(graph)[igraph::degree(graph, mode = 'out') > 0 & igraph::degree(graph, mode = 'in') > 0])
    df[i, "perc_int"] <- (n_int / S) * 100
    
    # 6. % B (Basal species that have consumers but lack resources)
    # Basal species: Species with in-degree > 0 (they are consumed) and out-degree = 0 (they do not consume anything)
    n_basals <- length(V(graph)[igraph::degree(graph, mode = 'in') > 0 & igraph::degree(graph, mode = 'out') == 0])
    df[i, "perc_basals"] <- (n_basals / S) * 100
    
    # 7. %C (Cannibal species that eat themselves)
    edge_list <- as_edgelist(graph, names = TRUE)
    n_cannibals <- length(which(edge_list[, 1] == edge_list[, 2]))
    df[i, "perc_cannibals"] <- (n_cannibals / S) * 100
    
    # 8. Omnivory (based on trophic levels)
    # Compute the trophic level for each node  
    tlnodes <- ego_size(graph, mode = "out")
    
    # Link the trophic level to the interactions
    netmatrix <- as_adjacency_matrix(graph, sparse = FALSE)
    webtl <- netmatrix * as.vector(tlnodes)
    
    # Remove the trophic level when no interactions
    webtl[webtl == 0] <- NA
    
    # Compute the standard deviation of the trophic levels of prey
    omninodes <- apply(webtl, 2, sd, na.rm = TRUE)
    
    # Average the standard deviation over all species (excluding those with < 2 preys)
    df[i, "omnivory"] <- mean(omninodes, na.rm = TRUE)
    
    # 9. Generality (mean and SD of incoming links)
    pred <- igraph::degree(graph, mode = "in") > 0  # Identify predator nodes (with at least one prey)
    df[i, "mean_gen"] <- mean(igraph::degree(graph, mode = "in")[pred], na.rm = TRUE)
    df[i, "sd_gen"] <- sd(igraph::degree(graph, mode = "in")[pred], na.rm = TRUE)
    
    # 10. Vulnerability (mean and SD of outgoing links)
    prey <- igraph::degree(graph, mode = "out") > 0  # Identify prey nodes (with at least one predator)
    df[i, "mean_vul"] <- mean(igraph::degree(graph, mode = "out")[prey], na.rm = TRUE)
    df[i, "sd_vul"] <- sd(igraph::degree(graph, mode = "out")[prey], na.rm = TRUE)
  }
  
  # Add simulation column for record-keeping
  df$simulation <- seq_len(nrow(df))
  
  return(df)
}





# For simulated networks for NEUTRAL scenarios, that come as adjacency matrices:

get_network_measures_from_adjacency_neutral <- function(list_adj_matrix, vec_p_est, vec_p_ext) {
  
  columns <- c("p(est)",
               "p(ext)",
               "S",
               "Link_density",
               "C",
               "perc_tops",
               "perc_int",
               "perc_basals",
               "perc_cannibals",
               "omnivory",
               "mean_gen",
               "sd_gen",
               "mean_vul",
               "sd_vul")
  
  df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df) = columns
  
  
  for (i in 1:length(list_adj_matrix)) {
    
    graph <- graph_from_adjacency_matrix(adjmatrix = list_adj_matrix[[i]],
                                         mode = "directed")
    
    
    # 1. Number of species within the foodweb (S)
    
    S <- gorder(graph)
    
    df[i, "S"] <- S
    
    
    
    
    # 2. links density (L/S)
    
    L <- gsize(graph)
    
    df[i, "Link_density"] <- L/S
    
    
    
    
    # 3. connectance (C = L/S2)
    
    df[i, "C"] <- L/(S*(S-1)/2)
    
    
    
    
    # 4. %T (Top species that have resource species but lack any consumer species)
    
    n_tops <- length(V(graph)[degree(graph, mode = 'in')>0 & degree(graph, mode = 'out') == 0])
    
    df[i, "perc_tops"] <- n_tops/S*100
    
    
    
    
    
    #5. % I (Intermediate species that have both resource and consumer species)
    
    n_int <- length(V(graph)[degree(graph, mode = 'out')>0 & degree(graph, mode = 'in') > 0])
    
    df[i, "perc_int"] <- n_int/S*100
    
    
    
    
    #6. % B (Basal species that have consumer species but lack resources species)
    
    n_basals <- length(V(graph)[degree(graph, mode = 'out')>0 & degree(graph, mode = 'in') == 0])
    
    df[i, "perc_basals"] <- n_basals/S*100
    
    
    
    
    # 7. %C (Cannibal species that eat themselves)
    
    edge_list <- as_edgelist(graph, names = TRUE)
    
    n_cannibals <- length(which(edge_list[,1] == edge_list[,2]))
    
    df[i, "perc_cannibals"] <- n_cannibals/S*100
    
    
    
    
    
    #8. %Omn (species that eat species at different trophic levels) -> changed for omnivory (general omnivory of the foodweb)
    
    # Compute the trophic level for each node  
    tlnodes <- neighborhood.size(graph,mode="out")
    
    # Calculate the average trophic level of the food web
    TL <- mean(tlnodes)
    
    #Omnivory is based on the calculation of trophic levels, and corresponds to the standard deviation of the trophic levels of a species' prey.
    
    netmatrix <- get.adjacency(graph, sparse=F)
    # Link the trophic level to the interactions
    webtl <- netmatrix*as.vector(tlnodes)
    # Remove the trophic level when no interactions
    webtl[webtl==0] <- NA
    
    #Compute the standard of the trophic levels of prey 
    omninodes <- apply(webtl,2,sd, na.rm=TRUE)
    
    # Average the standard deviation over all taxa (with more than 2 preys)  
    df[i, "omnivory"] <- mean(omninodes, na.rm=TRUE)
    
    
    
    
    
    
    #9. Standard deviation of mean generality (GenSD)
    # Generality
    
    pred <- degree(graph, mode="in")>0 # Identify predator nodes, i.e. taxa with at least one prey
    
    mean_gen <- mean(degree(graph, mode="in")[pred])
    
    df[i, "mean_gen"] <- mean(degree(graph, mode="in")[pred])
    df[i, "sd_gen"] <- sd(degree(graph, mode="in")[pred])
    
    
    
    
    
    # 10. Standard deviation of mean vulnerability (VulSD)
    # Generality
    
    prey <- degree(graph, mode="out")>0 
    
    df[i, "mean_vul"] <- mean(degree(graph, mode="out")[prey])
    df[i, "sd_vul"] <- sd(degree(graph, mode="out")[prey])
    
    
  }
  
  
  
  df[,"p(est)"] <- vec_p_est
  df[,"p(ext)"] <- vec_p_ext
  
  
  return(df)
  
}



# For empirical networks that come as incidence matrices:


get_network_measures_from_incidence <- function(list_incidence_matrix, ref_id) {
  
  columns <- c("id",
               "S",
               "Link_density",
               "C",
               "perc_tops",
               "perc_int",
               "perc_basals",
               "perc_cannibals",
               "omnivory",
               "mean_gen",
               "sd_gen",
               "mean_vul",
               "sd_vul")
  
  df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df) = columns
  
  
  for (i in 1:length(list_incidence_matrix)) {
    
    graph <- graph_from_adjacency_matrix(adjmatrix = list_adj_matrix[[i]],
                                         mode = "directed",
                                         weighted = NULL)  # Ensuring it's binary (0/1)
    
    
    # 1. Number of species within the foodweb (S)
    
    S <- vcount(graph)
    
    df[i, "S"] <- S
    
    
    
    
    # 2. links density (L/S)
    
    L <- gsize(graph)
    
    df[i, "Link_density"] <- L/S
    
    
    
    
    # 3. connectance (C = L/S2)
    
    df[i, "C"] <- L/(S*(S-1)/2)
    
    
    
    
    # 4. %T (Top species that have resource species but lack any consumer species)
    
    n_tops <- length(V(graph)[degree(g, mode = 'in')>0 & degree(graph, mode = 'out') == 0])
    
    df[i, "perc_tops"] <- n_tops/S*100
    
    
    
    
    
    #5. % I (Intermediate species that have both resource and consumer species)
    
    n_int <- length(V(graph)[degree(graph, mode = 'out')>0 & degree(graph, mode = 'in') > 0])
    
    df[i, "perc_int"] <- n_int/S*100
    
    
    
    
    #6. % B (Basal species that have consumer species but lack resources species)
    
    n_basals <- length(V(graph)[degree(graph, mode = 'out')>0 & degree(graph, mode = 'in') == 0])
    
    df[i, "perc_basals"] <- n_basals/S*100
    
    
    
    
    # 7. %C (Cannibal species that eat themselves)
    
    edge_list <- as_edgelist(graph, names = TRUE)
    
    n_cannibals <- length(which(edge_list[,1] == edge_list[,2]))
    
    df[i, "perc_cannibals"] <- n_cannibals/S*100
    
    
    
    
    
    #8. %Omn (species that eat species at different trophic levels) -> changed for omnivory (general omnivory of the foodweb)
    
    # Compute the trophic level for each node  
    # tlnodes <- neighborhood.size(graph,mode="out")
    # 
    # # Calculate the average trophic level of the food web
    # TL <- mean(tlnodes)
    # 
    # #Omnivory is based on the calculation of trophic levels, and corresponds to the standard deviation of the trophic levels of a species' prey.
    # 
    # netmatrix <- get.adjacency(graph, sparse=F)
    # # Link the trophic level to the interactions
    # webtl <- netmatrix*as.vector(tlnodes)
    # # Remove the trophic level when no interactions
    # webtl[webtl==0] <- NA
    # 
    # #Compute the standard of the trophic levels of prey 
    # omninodes <- apply(webtl,2,sd, na.rm=TRUE)
    # 
    # # Average the standard deviation over all taxa (with more than 2 preys)  
    # df[i, "omnivory"] <- mean(omninodes, na.rm=TRUE)
    
    trophic_levels <- igraph::estimate_trophic_levels(graph, mode = "in")
    # Omnivory is the standard deviation of the trophic levels of a species' prey
    prey_tl <- matrix(NA, nrow = S, ncol = S)
    for (pred in 1:S) {
      prey_tl[pred, ] <- trophic_levels[which(adjmatrix[, pred] > 0)]
    }
    omnivory <- apply(prey_tl, 1, sd, na.rm = TRUE)
    df[i, "omnivory"] <- mean(omnivory, na.rm = TRUE)
    
    
    
    
    
    
    # #9. Standard deviation of mean generality (GenSD)
    # # Generality
    # 
    # pred <- degree(graph, mode="in")>0 # Identify predator nodes, i.e. taxa with at least one prey
    # 
    # mean_gen <- mean(degree(graph, mode="in")[pred])
    # 
    # df[i, "mean_gen"] <- mean(degree(graph, mode="in")[pred])
    # df[i, "sd_gen"] <- sd(degree(graph, mode="in")[pred])
    # 
    # 
    # 
    # 
    # 
    # # 10. Standard deviation of mean vulnerability (VulSD)
    # # Generality
    # 
    # prey <- degree(graph, mode="out")>0 
    # 
    # df[i, "mean_vul"] <- mean(degree(graph, mode="out")[prey])
    # df[i, "sd_vul"] <- sd(degree(graph, mode="out")[prey])
    
    generality <- degree(graph, mode = "in")  # Number of prey
    vulnerability <- degree(graph, mode = "out")  # Number of predators
    
    df[i, "mean_gen"] <- mean(generality[generality > 0])
    df[i, "sd_gen"] <- sd(generality[generality > 0])
    
    df[i, "mean_vul"] <- mean(vulnerability[vulnerability > 0])
    df[i, "sd_vul"] <- sd(vulnerability[vulnerability > 0])
    
    
  }
  
  
  
  df[,"id"] <- ref_id
  
  
  return(df)
  
  
  
}





# Compute motifs df and plot from adjacency networks (of neutral scenarios)

compute_motifs_from_adjacency <- function(list_adj_matrix) {
  
  
  columns <- as.character(c(1:16))
  df_motifs = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df_motifs) = columns
  
  
  
  for (i in 1:length(list_adj_matrix)) {
    
    graph <- graph_from_adjacency_matrix(adjmatrix = list_adj_matrix[[i]],
                                         mode = "directed")
    
    
    v_motifs <- motifs(graph, 3)
    
    df_motifs[i,] <- v_motifs
    
  }
  
  
  df_motifs <- Filter(function(x)!all(is.na(x)), df_motifs)
  
  df_motifs$scenario <- chartr("123456789", "ABCDEFGHI", c(1:nrow(df_motifs)))
  
  
  df_melted = melt(df_motifs, id.vars = 'scenario')
  
  p <- ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = scenario, group = scenario))+
    theme_bw()+
    my.theme+
    theme(legend.position = "right")+
    xlab("triad isomorphism classes")+
    ylab("frequency")
  
  list_results <- list(df_motifs, p)
  
  return(list_results)
  
}




## Compute motifs from incidence matrix (empirical networks)


compute_motifs_from_incidende <- function(list_incidence_matrix, ref_id) {
  
  columns <- as.character(c(1:16))
  df_motifs = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df_motifs) = columns
  
  
  
  for (i in 1:length(list_files)) {
    
    graph <- graph_from_incidence_matrix(incidence = list_incidence_matrix[[i]],
                                         directed = TRUE)
    
    
    v_motifs <- motifs(graph, 3)
    
    df_motifs[i,] <- v_motifs
    
  }
  
  
  df_motifs <- Filter(function(x)!all(is.na(x)), df_motifs)
  
  df_motifs$scenario <- ref_id
  
  
  df_melted = melt(df_motifs, id.vars = 'scenario')
  
  p <- ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = scenario, group = scenario))+
    theme_bw()+
    my.theme+
    theme(legend.position = "right")+
    xlab("triad isomorphism classes")+
    ylab("frequency")
  
  list_result <- list(df_motifs, p)
  
  return(list_result)
  
}





# Reshape table with metrics


reshape_df_net.metrics  <- function(df_metrics, interaction) {
  
  
  df_mean_net_metrics <- data.frame(names(colMeans(df_metrics)),colMeans(df_metrics))
  rownames(df_mean_net_metrics) <- NULL
  
  sd <- apply(df_metrics, 2, sd)
  
  se <-  sd/sqrt(nrow(df_metrics))
  
  df_mean_net_metrics$se <- se
  
  colnames(df_mean_net_metrics) <- c("metric", interaction, paste("se_", interaction))
  
  return(df_mean_net_metrics)
  
  
}


# Reshape table with motif frequencies





reshape_df_motifs <- function(df_motifs, interaction) {
  
  df_motifs$scenario <- NULL
  
  df_mean_motifs <- data.frame(names(colMeans(df_motifs)),colMeans(df_motifs))
  rownames(df_mean_motifs) <- NULL
  
  colnames(df_mean_motifs) <- c("motif", interaction)
  
  return(df_mean_motifs)
  
}



# Plot metric from summary mean metrics table



plot_metric <- function(metric, y_axis) {
  
  df <- df_metrics[df_metrics$metric == metric,]
  
  df_plot <- df[,c("metric","facilitation", "competition", "foodweb", "neutral", "empirical")]
  
  df_se <- df[1,c("metric","se_ facilitation", "se_ competition", "se_ foodweb", "se_ neutral", "se_ empirical")]
  
  df_plot <- melt(df_plot, id = "metric")
  
  df_se_plot <- melt(df_se, id = "metric")
  
  colnames(df_se_plot)[colnames(df_se_plot) == "value"] <- "se"
  
  df_plot$se <- df_se_plot$se
  
  p <- ggplot(df_plot) +
    geom_bar( aes(x=variable, y=value, fill = variable), stat="identity", position = position_dodge(width = 0.2)) +
    geom_errorbar( aes(x=variable, ymin=value-se, ymax=value+se), width=0.1, alpha=0.9, size=1) +
    ggtitle(paste("metric:",  metric))+
    scale_fill_manual(values = c("gray50", "gray50", "gray50", "gray50","gray50", "gray82"))+
    xlab("")+
    ylab(y_axis)+
    theme_bw()+
    my.theme+
    theme(legend.position = "none")
  
  return(p)
  
  
  
  
  
}




plot_metric_foodweb <- function(df, metric, y_axis) {
  
  df <- df[df$metric == metric,]
  
  df_plot <- df[,c("metric","foodweb", "empirical")]
  
  df_se <- df[1,c("metric","se_ foodweb", "se_ empirical")]
  
  df_plot <- melt(df_plot, id = "metric")
  
  df_se_plot <- melt(df_se, id = "metric")
  
  colnames(df_se_plot)[colnames(df_se_plot) == "value"] <- "se"
  
  df_plot$se <- df_se_plot$se
  
  p <- ggplot(df_plot) +
    geom_bar( aes(x=variable, y=value, fill = variable), stat="identity", position = position_dodge(width = 0.2)) +
    geom_errorbar( aes(x=variable, ymin=value-se, ymax=value+se), width=0.1, alpha=0.9, size=1) +
    ggtitle(paste("metric:",  metric))+
    scale_fill_manual(values = c("gray50", "gray82"))+
    xlab("")+
    ylab(y_axis)+
    theme_bw()+
    my.theme+
    theme(legend.position = "none")
  
  return(p)
  
  
  
  
  
}



# set colnames and rownames as numbers


set_sppNames_numbers <- function(mat) {
  
  rownames(mat) <- seq(1:nrow(mat))
  colnames(mat) <- seq(1:ncol(mat))
  
  return(mat)  
}





#set_sppNames_icolrows <-  function(matrix) {

# mat <- matrix
# rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "")
# colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "")

#  return(mat)

#}



# convert colnames and rownames from numbers to letters

convert_sppnames_toletters <-  function(mat) {
  
  
  rownames(mat) <- chartr("0123456789", "ABCDEFGHIJ", rownames(mat))
  colnames(mat) <- chartr("0123456789", "ABCDEFGHIJ", colnames(mat))
  
  return(mat)
  
}



# convert species names of ancestry table from numbers to letters

change_sppnames_letters_ancdist.table <- function(mat) {
  
  mat[,"ancestor"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"ancestor"])
  
  mat[,"spp"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"spp"])
  
  return(mat)
  
}