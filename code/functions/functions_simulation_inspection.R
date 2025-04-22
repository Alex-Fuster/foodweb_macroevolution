
####### Function to plot species richness over time

plot_spp_richness <- function(steps, pres_mat, title) {
  
  timesteps <- 1:steps
  
  n_spp_time <- c()
  
  for (i in 1:nrow(pres_mat)) {
    
    n_spp_time[i] <- length(which(pres_mat[i,] == 1))
    
  }
  
  
  df_spprich_time <- data.frame(timesteps,n_spp_time)
  
  p.spp <- ggplot(df_spprich_time, aes(x=timesteps, y = n_spp_time)) +
    geom_line(color="black", linetype="twodash") +
    theme_classic()+
    my.theme+
    xlab("timesteps")+
    ylab("N species")+
    ggtitle(title)
  
  return(p.spp)
}


####### function to compute plots for inspecting results from simulation

inspect_simulation_fac_comp <- function(simulation_data, interaction_type) {
  
  n_steps <- length(simulation_data$network_list)
  pres = simulation_data$presence_matrix
  
  
  
  
  # Niche values through time
  
  res_spp = data.frame(matrix(ncol=2,
                              nrow=0, 
                              dimnames=list(NULL, c("step", "n values")))) #makes an empty dataframe
  
  res_mean = data.frame(matrix(ncol=2,
                               nrow=0, 
                               dimnames=list(NULL, c("step", "n values")))) 
  
  
  for (i in 1:nsteps) {
    
    pres_stepi <- which(pres[i,] == 1)
    traits_stepi <- simulation_data$traits_df[pres_stepi, "n"]
    vec_step <- rep(i, times = length(traits_stepi))
    res.add <- data.frame("step" = vec_step, 
                          "n values" = traits_stepi)
    res_spp <- rbind(res_spp, res.add)
    res.add1 <- data.frame("step" = i, 
                           "n values" = mean(traits_stepi))
    res_mean <- rbind(res_mean, res.add1)
    
  }
  
  plot_nichetime <- ggarrange(
    
    ggplot(res_spp, aes(x = step, y = n.values)) +
      geom_point(alpha = 0.2)+
      ggtitle("spp n values")+
      theme_classic()+
      my.theme+
      ylim(0,1),
    
    ggplot(res_mean, aes(x = step, y = n.values)) +
      geom_point(alpha = 0.2)+
      ggtitle("mean n values")+
      theme_classic()+
      my.theme+
      ylim(0,1),
    
    nrow = 2, 
    ncol = 1
  )
  
  
  
  
  
  # final degree distribution
  
  
  net_f <- simulation_data$network_list[[length(simulation_data$network_list)]]
  colnames(net_f) <- seq(1:1000)
  rownames(net_f) <- seq(1:1000)
  
  pres_f <- simulation_data$presence_matrix[length(simulation_data$network_list),]
  names(pres_f) <- seq(1:1000)
  net_f <- net_f[names(which(pres_f == 1)), names(which(pres_f == 1))]
  graph <- graph_from_adjacency_matrix(adjmatrix = net_f,
                                       mode = "directed")
  degree <- igraph::degree(graph, mode="all")
  deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  spp <- 1:length(deg.dist)
  df_degreeplot <- data.frame(deg.dist, spp)
  
  plot_degreedistr <- ggplot(df_degreeplot, aes(x = spp, y = deg.dist))+
    geom_point() +
    theme_classic()+
    my.theme+
    ylab("degree")+
    xlab("spp")
  
  
  
  # species richness
  
  
  plot_richnesstime <- plot_spp_richness(steps = n_steps, pres_mat = simulation_data$presence_matrix[1:n_steps,], title = "foodweb")
  
  
  ##  Speciation and extinction curves
  
  t0 = pres[1:(n_steps-1),]
  t1 = pres[2:n_steps,] 
  spec_mat = pres[1:(n_steps-1),]*0
  ext_mat = pres[1:(n_steps-1),]*0
  spec_mat[t1-t0==1] = 1 
  ext_mat[t1-t0==-1] = 1
  spec = apply(spec_mat ,1,sum)
  ext = apply(ext_mat, 1, sum)
  S = apply(pres,1,sum)[2:n_steps]
  time =  1:(n_steps-1)
  
  df_divrates <- data.frame(spec, ext, S, time)
  df_plot_divrates <- melt(df_divrates, id = c("time", "S"))
  
  
  # plot diversification rates
  
  plot_divrates <- ggplot(data = df_plot_divrates[df_plot_divrates$variable %in% c("spec", "ext"),], aes(x = time, y = value / S, group = variable)) +
    geom_smooth(aes(color = variable), method = "loess") +
    scale_color_manual(values = c("blue4", "red2")) +
    theme_classic() +
    ylab("Rate") +
    labs(color = "")
  
  
  ## diversification-richness dependence
  
  plot_div_richness <- ggplot(data=df_plot_divrates, aes(x=S, y=value/S, group=variable)) +
    geom_point(aes(color = variable), size = 2, alpha = 0.4)+
    scale_color_manual(values = c("blue4", "red2"))+
    theme_classic()+
    my.theme+
    ylab("rate")+
    xlab("Species richness")+
    labs(color = "")
  
  
  ##########################################
  
  arrange_plots <- ggarrange(
    plot_nichetime + ggtitle("evolution niche values"),
    plot_degreedistr + ggtitle("degree distr. final"),
    plot_richnesstime + ggtitle("evolution richness"),
    plot_divrates + ggtitle("evolution div. rates"),
    plot_div_richness +ggtitle("div.rates & richness")
  )
  
  
  arrange_plots <- annotate_figure(arrange_plots, top = text_grob(interaction_type, size = 14, face = "bold"))
  
  return(arrange_plots)
  
  
  
}





####################


inspect_simulation_fw <- function(simulation_data, nbasals, Smax) {
  
  
  # count number of timesteps where there are spp
  
  n_steps <- length(simulation_data$network_list)
  
  pres = simulation_data$presence_matrix
  res_spp = data.frame(matrix(ncol=2,
                              nrow=0, 
                              dimnames=list(NULL, c("step", "n values")))) #makes an empty dataframe
  
  res_mean = data.frame(matrix(ncol=2,
                               nrow=0, 
                               dimnames=list(NULL, c("step", "n values")))) 
  
  
  for (i in 1:nsteps) {
    
    pres_stepi <- which(pres[i,] == 1)
    traits_stepi <- simulation_data$traits_df[pres_stepi, "n"]
    
    vec_step <- rep(i, times = length(traits_stepi))
    res.add <- data.frame("step" = vec_step, 
                          "n values" = traits_stepi)
    res_spp <- rbind(res_spp, res.add)
    res.add1 <- data.frame("step" = i, 
                           "n values" = mean(traits_stepi))
    res_mean <- rbind(res_mean, res.add1)
    
    
  }
  
  plot_nichetime <- ggarrange(
    
    ggplot(res_spp, aes(x = step, y = n.values)) +
      geom_point(alpha = 0.2)+
      ggtitle("spp n values")+
      theme_classic()+
      my.theme+
      ylim(0,1),
    
    ggplot(res_mean, aes(x = step, y = n.values)) +
      geom_point(alpha = 0.2)+
      ggtitle("mean n values")+
      theme_classic()+
      my.theme+
      ylim(0,1),
    
    nrow = 2, 
    ncol = 1
    
    
  )
  
  
  
  net_f <- simulation_data$network_list[[length(simulation_data$network_list)]]
  
  
  net_f <- net_f[(nbasals+1):(Smax+nbasals),]
  
  colnames(net_f) <- seq(1:1000)
  rownames(net_f) <- seq(1:1000)
  pres_f <- simulation_data$presence_matrix[length(simulation_data$network_list),]
  names(pres_f) <- seq(1:1000)
  net_f <- net_f[names(which(pres_f == 1)), names(which(pres_f == 1))]
  graph <- graph_from_adjacency_matrix(adjmatrix = net_f,
                                       mode = "directed")
  degree <- igraph::degree(graph, mode="all")
  deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  spp <- 1:length(deg.dist)
  df_degreeplot <- data.frame(deg.dist, spp)
  
  plot_degreedistr <- ggplot(df_degreeplot, aes(x = spp, y = deg.dist))+
    geom_point() +
    theme_classic()+
    my.theme+
    ylab("degree")+
    xlab("spp")
  
  
  
  ## Speces richness with time
  

  
  plot_richnesstime <- plot_spp_richness(steps = n_steps, pres_mat = simulation_data$presence_matrix[1:n_steps,], title = "foodweb")
  
  
  
  
  ## Speciation and extinction events
  
  
  # Number of speciation and extinction events
  t0 = pres[1:(n_steps-1),]
  t1 = pres[2:n_steps,] 
  spec_mat = pres[1:(n_steps-1),]*0
  ext_mat = pres[1:(n_steps-1),]*0
  spec_mat[t1-t0==1] = 1 
  ext_mat[t1-t0==-1] = 1
  spec = apply(spec_mat ,1,sum)
  ext = apply(ext_mat, 1, sum)
  
  S = apply(pres,1,sum)[2:n_steps]
  
  time =  1:(n_steps-1)
  
  # Compute the diversification rate
  div_rate <- (spec - ext) / S
  
  # Create data frames for plotting
  df_rates<- data.frame(time, S, spec, ext)
  df_plot_rates <- melt(df_rates, id.vars = c("time", "S"))
  
  
  plot_divrates_smoothed <- ggplot(data = df_plot_rates[df_plot_rates$variable %in% c("spec", "ext"),], aes(x = time, y = value / S, group = variable)) +
    geom_smooth(aes(color = variable), method = "loess") +
    scale_color_manual(values = c("blue4", "red2")) +
    theme_classic() +
    ylab("Rate") +
    labs(color = "")
  
  
  
  ## diversification-richness dependence
  
  
  plot_div_richness <- ggplot(data=df_plot_rates, aes(x=S, y=value/S, group=variable)) +
    geom_point(aes(color = variable), size = 2, alpha = 0.4)+
    scale_color_manual(values = c("blue4", "red2"))+
    theme_classic()+
    my.theme+
    ylab("rate")+
    xlab("Species richness")+
    labs(color = "")
  
  
  
  arrange_plots <- ggarrange(
    plot_nichetime + ggtitle("evolution niche values"),
    plot_degreedistr + ggtitle("degree distr. final"),
    plot_richnesstime + ggtitle("evolution richness"),
    plot_divrates_smoothed + ggtitle("evolution div. rates"),
    plot_div_richness +ggtitle("div.rates & richness")
  )
  
  arrange_plots <- annotate_figure(arrange_plots, top = text_grob("foodweb", size = 14, face = "bold"))
  
  return(arrange_plots)
  
}
