# get_BT_data reads BT (BayesTraits) output and discards the first 500 rows
# plot_HPD plots locations using colors that 
# reflect densities and adds a line representing the 95% HPD interval
# plot_HPD_multi does the same as plot_HPD but allows for multiple 
# cases in a single map

# filename <- "abkh1242_AncStates.txt"
# node number is the number used as ID by BT

# this auxiliary function loops through all BT output files and
# plots homelands
allfams <- function() {
  filenames <- dir("BT_outfiles")
  for (i in 1:length(filenames)) {
    run(filenames[i])
  }
}

# this auxiliary function runs the functions below
# uncomment the line for saving plot if needed
run <- function(filename, node=0) {
  node <<- node
  BT_data <- get_BT_data(filename, node)
  plot_HPD(BT_data)
  # save_plot(filename)
}

get_BT_data <- function(filename, node=0) {
  # read the BT output file
  rl <- readLines(paste0("BT_outfiles/", filename))
  w_node <- grep("^Node", rl)
  rl2 <- rl[-w_node]
  writeLines(rl2, con="tmp.txt")
  BT_table_full <- read.table(file="tmp.txt", header=TRUE, sep="\t")
  file.remove("tmp.txt")
  # discard the first 500 rows
  BT_table <- BT_table_full[501:1000,]
  # get columns corresponding to the node desired
  node <- as.character(node)
  L <- nchar(node)
  zeros <- 5 - L
  node_name <- paste0(paste(rep(0, zeros),collapse=""),node)
  iteration <- BT_table$Itter
  likelihood <- BT_table$Lh
  lat_col_name <- paste0("Node.", node_name, "...Lat")
  w_lat_col_name <- which(names(BT_table)==lat_col_name)
  latitude <- BT_table[,w_lat_col_name]
  lon_col_name <- paste0("Node.", node_name, "...Long")
  w_lon_col_name <- which(names(BT_table)==lon_col_name)
  longitude <- BT_table[,w_lon_col_name]
  BT_data <- data.frame(iteration, likelihood, latitude, longitude)
  return(BT_data)
}


library(MASS)  # kde2d()
library(maps)  # map_data()
library(ggplot2)  # ggplot()

# plot a single homeland
plot_HPD <- function(BT_data) {
  longitude <- BT_data$longitude
  latitude <- BT_data$latitude
  posterior_samples <- data.frame(longitude, latitude)
  
  # if all coordinates are identical kernel density estimation will not work
  if(length(unique(paste(longitude, latitude)))==1) {
    stop("There is only one coordinate, cannot compute an HPD interval")
  }

  # Kernel density estimation
  density_estimate <- kde2d(
    x = posterior_samples$longitude,
    y = posterior_samples$latitude,
    n = 200 # Resolution of the grid
  )

  # Compute density values for the sample points
  # density_values <- MASS::kde2d(
  #  posterior_samples$longitude, 
  #  posterior_samples$latitude, 
  #  n = 200
  # )$z

  # Map density values back to the points
  posterior_samples$density <- apply(posterior_samples, 1, function(row) {
    x_idx <- which.min(abs(density_estimate$x - row["longitude"]))
    y_idx <- which.min(abs(density_estimate$y - row["latitude"]))
    density_estimate$z[x_idx, y_idx]
  })

  # Convert density to a vector
  densities <- as.vector(density_estimate$z)

  # Sort densities in descending order
  sorted_densities <- sort(densities, decreasing = TRUE)

  # Compute cumulative probabilities
  cumulative_prob <- cumsum(sorted_densities) / sum(sorted_densities)

  # Find threshold density for the 95% HPD region
  threshold_density <- sorted_densities[which.max(cumulative_prob >= 0.95)]

  # Identify grid points in the HPD region
  # hpd_region <- which(density_estimate$z >= threshold_density, arr.ind = TRUE)

  # Prepare grid for contour plotting
  grid_data <- expand.grid(
    longitude = density_estimate$x,
    latitude = density_estimate$y
  )
  grid_data$density <- as.vector(density_estimate$z)
  grid_data$in_hpd <- grid_data$density >= threshold_density

  # Compute the map boundaries with padding
  xlim <- range(posterior_samples$longitude, na.rm = TRUE) + c(-5, 5)  # Add padding
  ylim <- range(posterior_samples$latitude, na.rm = TRUE) + c(-5, 5)   # Add padding
  
  # Base map using the 'world' map from the maps package
  world_map <- map_data("world")
  
  # Plot the results
  ggplot() +
    # Plot the base map
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
    # Plot points with color based on density
    geom_point(data = posterior_samples, aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
    # Set map limits to focus on the region
    coord_fixed(xlim = xlim, ylim = ylim) +
    # Add a contour for the 95% HPD region
    geom_contour(data = grid_data, aes(x = longitude, y = latitude, z = density), 
                 breaks = threshold_density, color = "red") +
    # Add color scale for density
    scale_color_viridis_c(name = "Density") +
    # Labels and theme
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()

}

save_plot <- function(filename) {
  family <- strsplit(filename, "_")[[1]][1]
  plotfile <- paste0("BT_HPD_homeland_plots/", family, ".png")
  ggsave(
    plotfile,
    device="png",
    create.dir=TRUE
  )
}

# takes as input a list of BT_data-type objects of length > 1; these are prepared
# separately using get_BT_data and combining those in a list
# lines inside the function bracket may have to be run
# since calling the function sometimes (at least) doesn't work
plot_HPD_multi <- function(BTL) {  # BTL stands for something like BT_data_list
  
  all_lon <- c()
  all_lat <- c()
  for (i in 1:length(BTL)) {
    all_lon <- c(all_lon, BTL[[i]]$longitude)
    all_lat <- c(all_lat, BTL[[i]]$latitude)
  }
  
  # Compute the map boundaries with padding
  xlim <- range(all_lon, na.rm = TRUE) + c(-3, 3)  # Add padding
  ylim <- range(all_lat, na.rm = TRUE) + c(-3, 3)   # Add padding
  
  posterior_samples_list <- list()
  grid_data_list <- list()
  
  for (i in 1:length(BTL)) {
    BT_data <- BTL[[i]]
    longitude <- BT_data$longitude
    latitude <- BT_data$latitude
    posterior_samples <- data.frame(longitude, latitude)
    
    # if all coordinates are identical kernel density estimation will not work
    if(length(unique(paste(longitude, latitude)))==1) {
      stop("There is only one coordinate for case", i, "cannot compute an HPD interval\n")
    }
    
    # Kernel density estimation
    density_estimate <- kde2d(
      x = posterior_samples$longitude,
      y = posterior_samples$latitude,
      n = 200 # Resolution of the grid
    )
    
    # Map density values back to the points
    posterior_samples$density <- apply(posterior_samples, 1, function(row) {
      x_idx <- which.min(abs(density_estimate$x - row["longitude"]))
      y_idx <- which.min(abs(density_estimate$y - row["latitude"]))
      density_estimate$z[x_idx, y_idx]
    })
    
    # Convert density to a vector
    densities <- as.vector(density_estimate$z)
    
    # Sort densities in descending order
    sorted_densities <- sort(densities, decreasing = TRUE)
    
    # Compute cumulative probabilities
    cumulative_prob <- cumsum(sorted_densities) / sum(sorted_densities)
    
    # Find threshold density for the 95% HPD region
    threshold_density <- sorted_densities[which.max(cumulative_prob >= 0.95)]
    
    # Identify grid points in the HPD region
    # hpd_region <- which(density_estimate$z >= threshold_density, arr.ind = TRUE)
    
    # Prepare grid for contour plotting
    grid_data <- expand.grid(
      longitude = density_estimate$x,
      latitude = density_estimate$y
    )
    grid_data$density <- as.vector(density_estimate$z)
    grid_data$in_hpd <- grid_data$density >= threshold_density

    # put the data in the lists
    grid_data_list[[i]] <- grid_data
    posterior_samples_list[[i]] <- posterior_samples
    
  }
  
  # Base map using the 'world' map from the maps package
  world_map <- map_data("world")
  
  # Plot the results
  if (length(BTL)==2) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==3) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==4) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "blue") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "blue") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==5) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[5]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[5]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==6) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[5]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[6]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[5]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[6]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==7) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[5]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[6]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[7]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[5]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[6]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[7]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==8) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[5]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[6]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[7]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[8]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[5]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[6]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[7]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[8]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
  if (length(BTL)==9) {
    ggplot() +
      # Plot the base map
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
      # Set map limits to focus on the region
      coord_fixed(xlim = xlim, ylim = ylim) +
      # Plot points with color based on density
      geom_point(data = posterior_samples_list[[1]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[2]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[3]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[4]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[5]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[6]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[7]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[8]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      geom_point(data = posterior_samples_list[[9]], aes(x = longitude, y = latitude, color = density), size = 1, alpha = 0.8) +
      # Add a contour for the 95% HPD region
      geom_contour(data = grid_data_list[[1]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[2]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[3]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[4]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[5]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[6]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[7]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[8]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      geom_contour(data = grid_data_list[[9]], aes(x = longitude, y = latitude, z = density), 
                   breaks = threshold_density, color = "red") +
      # Add color scale for density
      scale_color_viridis_c(name = "Density") +
      # Labels and theme
      # labs(
      #  x = "Longitude",
      #  y = "Latitude"
      # ) +
      theme_minimal() +
      theme(legend.position = "none")  # Remove the legend
  }
}

# Plotting Indo-European
# source("get_node.R")  # get_node_id() function
get_node_id("Indo-European,ClassicalIndo-European", "indo1319_AncStates.txt")  # 19
dat1 <- get_BT_data("indo1319_AncStates.txt", 19)
plot_HPD(dat1)
get_node_id("Indo-European,ClassicalIndo-European,Indo-Iranian", "indo1319_AncStates.txt")  # 269
dat2 <- get_BT_data("indo1319_AncStates.txt", 269)
plot_HPD(dat2)
get_node_id("Indo-European,ClassicalIndo-European,Balto-Slavic", "indo1319_AncStates.txt")  # 30
dat3 <- get_BT_data("indo1319_AncStates.txt", 30)
plot_HPD(dat3)
get_node_id("Indo-European,ClassicalIndo-European,Albanian", "indo1319_AncStates.txt")  # 20
dat4 <- get_BT_data("indo1319_AncStates.txt", 20)
plot_HPD(dat4)
comb <- list(dat1, dat2)
plot_HPD_multi(comb)

# Plotting Atlantic-Congo
# here "Atlantic-Congo" and "Atlantic-Congo,North-CentralAtlantic" will not overlap
get_node_id("Atlantic-Congo", "atla1278_AncStates.txt")  # 0
dat1 <- get_BT_data("atla1278_AncStates.txt", 0)
plot_HPD(dat1)

# get_node_id("Atlantic-Congo,Mel", "atla1278_AncStates.txt")  # 5
# dat2 <- get_BT_data("atla1278_AncStates.txt", 5)
# plot_HPD(dat2)  # overlaps with Atlantic-Congo

get_node_id("Atlantic-Congo,North-CentralAtlantic", "atla1278_AncStates.txt")  # 20
dat3 <- get_BT_data("atla1278_AncStates.txt", 20)
plot_HPD(dat3)  # appears not to overlap with Atlantic-Congo

# all four of the following are non-overlapping
get_node_id("Atlantic-Congo,Volta-Congo", "atla1278_AncStates.txt")  # 96
dat4 <- get_BT_data("atla1278_AncStates.txt", 96)
plot_HPD(dat4)  # appears not to overlap with Atlantic-Congo

# get_node_id("Atlantic-Congo,Volta-Congo,NorthVolta-Congo", "atla1278_AncStates.txt")  # 1731
# dat5 <- get_BT_data("atla1278_AncStates.txt", 1731)
# plot_HPD(dat5)  # may overlap with Atlantic-Congo,Volta-Congo

# get_node_id("Atlantic-Congo,Volta-Congo,KwaVolta-Congo", "atla1278_AncStates.txt")  # 1731
# dat6 <- get_BT_data("atla1278_AncStates.txt", 1608)
# plot_HPD(dat6)  # no overlap with Atlantic-Congo,Volta-Congo

get_node_id("Atlantic-Congo,Volta-Congo,Benue-Congo", "atla1278_AncStates.txt")  # 97
dat7 <- get_BT_data("atla1278_AncStates.txt", 97)
plot_HPD(dat7)  # apparently some  overlap with Atlantic-Congo,Volta-Congo

# plotting Afro-Asiatic
# the three homelands will not overlap
get_node_id("Afro-Asiatic", "afro1255_AncStates.txt")  # 0
dat8 <- get_BT_data("afro1255_AncStates.txt", 0)
plot_HPD(dat8)  #
get_node_id("Afro-Asiatic,Chadic", "afro1255_AncStates.txt")  # 36
dat9 <- get_BT_data("afro1255_AncStates.txt", 36)
plot_HPD(dat9)  # no overlap with Afro-Asiatic
get_node_id("Afro-Asiatic,Cushitic", "afro1255_AncStates.txt")  # 367
dat10 <- get_BT_data("afro1255_AncStates.txt", 367)
plot_HPD(dat10)  # no overlap with Afro-Asiatic
get_node_id("Afro-Asiatic,Berber", "afro1255_AncStates.txt")  # 1
dat11 <- get_BT_data("afro1255_AncStates.txt", 1)
plot_HPD(dat11)  # no overlap with Afro-Asiatic


# plotting Arawakan
get_node_id("Arawakan", "araw1281_AncStates.txt")  # 0
dat1 <- get_BT_data("araw1281_AncStates.txt", 0)
plot_HPD(dat1)

get_node_id("Arawakan,Central-EasternMaipuran", "araw1281_AncStates.txt")  # 22
dat2 <- get_BT_data("araw1281_AncStates.txt", 22)
plot_HPD(dat2)  # small overlap with Arawakan

get_node_id("Arawakan,CaribbeanArawakan", "araw1281_AncStates.txt")  # 6
dat3 <- get_BT_data("araw1281_AncStates.txt", 6)
plot_HPD(dat3)

get_node_id("Arawakan,CaribbeanArawakan,AntilleanArawakan", "araw1281_AncStates.txt")  # 7
dat4 <- get_BT_data("araw1281_AncStates.txt", 7)
plot_HPD(dat4)  # small overlap with Arawakan,CaribbeanArawakan

get_node_id("Arawakan,CaribbeanArawakan,Guajiro-Paraujano", "araw1281_AncStates.txt")  # 16
dat5 <- get_BT_data("araw1281_AncStates.txt", 16)
plot_HPD(dat5)  # small overlap with Arawakan,CaribbeanArawakan

get_node_id("Arawakan,SouthernMaipuran", "araw1281_AncStates.txt")  # 83
dat6 <- get_BT_data("araw1281_AncStates.txt", 83)
plot_HPD(dat6)

get_node_id("Arawakan,SouthernMaipuran,BolivianArawakan", "araw1281_AncStates.txt")  # 84
dat7 <- get_BT_data("araw1281_AncStates.txt", 84)
plot_HPD(dat7)  # no overlap with the previous

get_node_id("Arawakan,SouthernMaipuran,Kampa-Amuesha", "araw1281_AncStates.txt")  # 93
dat8 <- get_BT_data("araw1281_AncStates.txt", 93)
plot_HPD(dat8)

get_node_id("Arawakan,SouthernMaipuran,Kampa-Amuesha,Pre-AndineMaipuran", "araw1281_AncStates.txt")  # 94
dat9 <- get_BT_data("araw1281_AncStates.txt", 94)
plot_HPD(dat9) # no overlap with the previous

# plotting Chibchan
get_node_id("Chibchan,CoreChibchan", "chib1249_AncStates.txt")  # not identified

get_node_id("Chibchan,CoreChibchan,Magdalenic", "chib1249_AncStates.txt")  # 17
dat1 <- get_BT_data("araw1281_AncStates.txt", 17)
plot_HPD(dat1)  # only one coordinate, no HPD interval

get_node_id("Chibchan,CoreChibchan,IsthmicChibchan,EasternIsthmicChibchan", "chib1249_AncStates.txt")  # not identified

# plotting Tupian
get_node_id("Tupian", "tupi1275_AncStates.txt")  # 0
dat1 <- get_BT_data("tupi1275_AncStates.txt", 0)
plot_HPD(dat1)

get_node_id("Tupian,EasternTupian", "tupi1275_AncStates.txt")  # 13
dat2 <- get_BT_data("tupi1275_AncStates.txt", 13)
plot_HPD(dat2)  # 
