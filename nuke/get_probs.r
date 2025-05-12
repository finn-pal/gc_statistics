# Load necessary libraries
library(hdf5r)
library(ks)

get_it_id <- function(it) {
  it_id <- sprintf("it%03d", it)
  return(it_id)
}

get_snap_id <- function(snapshot) {
  snap_id <- sprintf("snap%03d", snapshot)
  return(snap_id)
}

get_et_lz_norm_cov <- function(h5file, it, snapshot) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  # Access the data
  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et <- iteration_snapshot[["et"]][]
  lz <- iteration_snapshot[["lz"]][]

  data_matrix <- cbind(lz, et)

  num_data <- nrow(data_matrix)
  cov_matrix <- cov(data_matrix)

  # Set off-diagonal elements to 0
  cov_matrix[lower.tri(cov_matrix)] <- 0
  cov_matrix[upper.tri(cov_matrix)] <- 0

  h_bandwidth <- cov_matrix / sqrt(num_data)

  return(h_bandwidth)
}

get_all_h_bandwidths <- function(h5file, it_lst, snapshot) {
  lapply(it_lst, function(it) get_et_lz_norm_cov(h5file, it, snapshot))
}

get_average_matrix <- function(h_bandwidth_list) {
  num_matrices <- length(h_bandwidth_list)

  # Ensure the list is not empty
  if (num_matrices == 0) {
    stop("h_bandwidth_list is empty")
  }

  # Initialize the sum matrix with zeros of the same size as the first matrix
  sum_matrix <- Reduce(`+`, h_bandwidth_list)

  # Compute the average matrix
  avg_matrix <- sum_matrix / num_matrices

  return(avg_matrix)
}

plot_lz_et_iterations <- function(h_bandwidth_list) {
  # Extract the (1,1) elements from each matrix
  lz_values <- sapply(h_bandwidth_list, function(mat) sqrt(mat[1, 1]))
  et_values <- sapply(h_bandwidth_list, function(mat) sqrt(mat[2, 2]))

  # Set up the plotting area to have 1 row and 2 columns
  par(mfrow = c(1, 2))

  # Plot the distribution of delta L_z values without the default x-axis
  hist(lz_values,
    main = expression(paste(Delta, "L"[z], " Values")),
    xlab = expression(paste(Delta, "L"[z], "(10"^3, " kpc km s"^-1, ")")),
    breaks = 20,
    xaxt = "n", # Suppress the default x-axis
    xlim = c(0, max(lz_values) * 1.1) # Extend x-axis limit for scaling
  )

  # Get the current x-axis tick positions
  lz_ticks <- axTicks(1)

  # Scale the tick positions and labels by 10^3
  scaled_lz_ticks <- lz_ticks / 1000
  scaled_lz_labels <- scaled_lz_ticks

  # Set the x-axis with the scaled ticks and labels
  axis(1, at = lz_ticks, labels = scaled_lz_labels)

  # Plot the distribution of delta E_t values without the default x-axis
  hist(et_values,
    main = expression(paste(Delta, "E"[t], " Values")),
    xlab = expression(paste(Delta, "E"[t], "(10"^5, " km"^2, "s"^-2, ")")),
    breaks = 20,
    xaxt = "n", # Suppress the default x-axis
    xlim = c(1e4, max(et_values) * 1.1) # Extend x-axis limit for scaling
  )

  # Get the current x-axis tick positions
  et_ticks <- axTicks(1)

  # Scale the tick positions and labels by 10^5
  scaled_et_ticks <- et_ticks / 100000
  scaled_et_labels <- scaled_et_ticks

  # Set the x-axis with the scaled ticks and labels
  axis(1, at = et_ticks, labels = scaled_et_labels)

  # Reset the plotting layout to the default
  par(mfrow = c(1, 1))
}

plot_it_cont <- function(h5file, it, snapshot, group, h_bandwidth, cont_level) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  # Access the data
  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et <- iteration_snapshot[["et"]][]
  lz <- iteration_snapshot[["lz"]][]

  plot(lz, et,
    main = "Plot of lz vs et",
    xlab = expression(paste("L"[z], " Values")),
    ylab = expression(paste("E"[t], " Values")),
    pch = 19, # Use solid points
    col = "blue"
  ) # Optional color for the points


  group_id <- iteration_snapshot[["group_id"]][]
  et_0 <- et[group_id == group]
  lz_0 <- lz[group_id == group]

  kde <- kde(x = cbind(lz_0, et_0), H = h_bandwidth)

  plot(kde,
    cont = c(cont_level), xlab = "lz", ylab = "et",
    drawlabels = FALSE, col = "red",
    main = "KDE Overlay of lz vs et", add = TRUE
  )
}

get_grid_bounds <- function(h5file, it_lst, snapshot) {
  max_et_abs <- -Inf
  max_lz_abs <- -Inf

  for (it in it_lst) {
    t_id <- get_it_id(it)
    snap_id <- get_snap_id(snapshot)

    # Access the data
    iteration <- h5file[[t_id]]
    iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

    et_abs <- abs(iteration_snapshot[["et"]][])
    lz_abs <- abs(iteration_snapshot[["lz"]][])

    # Update maximum absolute values
    max_et_abs <- max(max_et_abs, max(et_abs, na.rm = TRUE))
    max_lz_abs <- max(max_lz_abs, max(lz_abs, na.rm = TRUE))
  }

  et_bound <- ceiling(max_et_abs / 100) * 100
  lz_bound <- ceiling(max_lz_abs / 100) * 100

  return(list(et_bound = et_bound, lz_bound = lz_bound))
}

# Example usage #####################################################################

# Specify the path to your HDF5 file
file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# Open the HDF5 file
h5file <- H5File$new(file_path, mode = "r")

it_lst <- seq(0, 100, by = 1)
snapshot <- 600

h_bandwidth_list <- get_all_h_bandwidths(h5file, it_lst, snapshot)
average_h_bandwidth <- get_average_matrix(h_bandwidth_list)


# TESTING #######################################################################

# it <- 0
# snapshot <- 600
# group <- 8580896
# # group <- 0

# n_grid <- 1000

# # it_lst <- seq(0, 2, by = 1)

# bounds <- get_grid_bounds(h5file, it_lst, snapshot)

# min_bounds <- c(-bounds$lz_bound, -bounds$et_bound)
# max_bounds <- c(bounds$lz_bound, 0)
# gridsize <- c(n_grid, n_grid)


# it_id <- get_it_id(it)
# snap_id <- get_snap_id(snapshot)

# # Access the data
# iteration <- h5file[[it_id]]
# iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

# et <- iteration_snapshot[["et"]][]
# lz <- iteration_snapshot[["lz"]][]

# group_id <- iteration_snapshot[["group_id"]][]
# et_0 <- et[group_id == group]
# lz_0 <- lz[group_id == group]

# kde <- kde(x = cbind(lz_0, et_0), H = average_h_bandwidth)
# kde_density <- kde$estimate

# # Normalize to ensure it's still a proper probability density function
# dx <- diff(kde$eval.points[[1]])[1]
# dy <- diff(kde$eval.points[[2]])[1]
# cell_area <- dx * dy
# kde_density <- kde_density / sum(kde_density * cell_area) # Normalize

# # Compute 75% cumulative probability threshold
# sorted_densities <- sort(as.vector(kde_density), decreasing = TRUE)
# cumulative_probs <- cumsum(sorted_densities * cell_area)
# threshold_cumulative <- sorted_densities[which(cumulative_probs >= 0.75)[1]]

# # Compute the density level at the 75% cumulative probability
# density_at_threshold <- threshold_cumulative

# # Verify normalization: The sum should be close to 1
# total_density <- sum(kde1_density * cell_area)
# print(total_density) # Should be approximately 1

# plot(lz, et,
#   main = "Plot of lz vs et",
#   xlab = expression(paste("L"[z], " Values")),
#   ylab = expression(paste("E"[t], " Values")),
#   pch = 19, # Use solid points
#   cex = 0.5, # Increase point size
#   col = "blue",
#   xlim = c(-bounds$lz_bound, bounds$lz_bound) * 1.1,
#   ylim = c(-bounds$et_bound, 0) * 1.1
# ) # Optional color for the points


# plot(kde, cont = c(99), col = "red", add = TRUE, drawlabels = FALSE, lwd = 2)

# # Add the contour based on the density at 75% cumulative probability (red)
# contour(kde$eval.points[[1]], kde$eval.points[[2]], kde_density,
#   levels = density_at_threshold, drawlabels = FALSE,
#   lwd = 1, col = "green", add = TRUE
# )



#####################################################################

# Combined function to process KDEs from multiple iterations and snapshots
process_kde_data <- function(
    h5file, it_lst, snapshot,
    group, n_grid = 1000, h_bandwidth) {
  # Get grid bounds
  bounds <- get_grid_bounds(h5file, it_lst, snapshot)

  min_bounds <- c(-bounds$lz_bound, -bounds$et_bound)
  max_bounds <- c(bounds$lz_bound, 0)
  gridsize <- c(n_grid, n_grid)

  # Initialize an empty density array for convolved density
  convolved_density <- matrix(0, nrow = gridsize[1], ncol = gridsize[2])

  # Initialize a list to store individual KDE results
  kde_results <- list()

  snap_id <- get_snap_id(snapshot)
  # Loop through each iteration in it_lst
  for (it in it_lst) {
    it_id <- get_it_id(it)
    print(it_id)

    # Access the data for the current iteration and snapshot
    iteration <- h5file[[it_id]]
    iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

    et <- iteration_snapshot[["et"]][]
    lz <- iteration_snapshot[["lz"]][]
    group_id <- iteration_snapshot[["group_id"]][]

    # Filter data based on the specified group
    et_0 <- et[group_id == group]
    lz_0 <- lz[group_id == group]

    # Compute KDE for the current iteration
    kde_it <- kde(
      x = cbind(lz_0, et_0), H = average_h_bandwidth,
      xmin = min_bounds, xmax = max_bounds, gridsize = gridsize
    )

    # Add the current KDE to the convolved density
    convolved_density <- convolved_density + kde_it$estimate

    # Store the KDE result for this iteration
    kde_results[[paste("it", it, sep = "_")]] <- kde_it
  }

  # Normalize the convolved density
  dx <- diff(kde$eval.points[[1]])[1]
  dy <- diff(kde$eval.points[[2]])[1]
  cell_area <- dx * dy
  convolved_density <- convolved_density / sum(convolved_density * cell_area)

  # Compute the 75% probability contour level
  sorted_densities <- sort(as.vector(convolved_density), decreasing = TRUE)
  cumulative_probs <- cumsum(sorted_densities * cell_area)
  threshold <- sorted_densities[which(cumulative_probs >= 0.25)[1]]

  # Verify normalization: The sum should be close to 1
  total_density <- sum(convolved_density * cell_area)
  print(total_density) # Should be approximately 1

  return(list(
    convolved_density = convolved_density,
    threshold = threshold,
    kde_results = kde_results
  ))
}

# Example structure for it_lst (list of iterations)
it_lst <- c(0, 1)

# Call the function with the necessary parameters
kde_results <- process_kde_data(h5file, it_lst, snapshot = 600, group = 0, n_grid = 1000, average_h_bandwidth)

# Access the results
convolved_density <- kde_results$convolved_density
threshold <- kde_results$threshold
individual_kdes <- kde_results$kde_results

# Define the bounds and grid size
min_bounds <- c(-bounds$lz_bound, -bounds$et_bound)
max_bounds <- c(bounds$lz_bound, 0)
gridsize <- c(n_grid, n_grid)

# Create the grid
x_grid <- seq(min_bounds[1], max_bounds[1], length.out = gridsize[1])
y_grid <- seq(min_bounds[2], max_bounds[2], length.out = gridsize[2])

# Plot heatmap (image) and overlay the 75% contour
image(x_grid, y_grid, convolved_density,
  col = terrain.colors(50), main = "75% Contour of Convolved Density",
  xlab = "X", ylab = "Y"
)

# Overlay the 75% contour
contour(x_grid, y_grid, convolved_density,
  levels = threshold, drawlabels = FALSE, lwd = 2, col = "red", add = TRUE
)
