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

get_initial_kde <- function(
    h5file, it, snapshot, bounds,
    group_id = FALSE, print_plot = FALSE, group_exl = FALSE) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  # Access the data
  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et <- iteration_snapshot[["et"]][]
  lz <- iteration_snapshot[["lz"]][]

  if (identical(group_id, FALSE)) {
    if (identical(group_exl, FALSE)) {
      et_grp <- et
      lz_grp <- lz
    } else {
      it_grp_id <- iteration_snapshot[["group_id"]][]
      it_et <- iteration_snapshot[["et"]][]
      it_lz <- iteration_snapshot[["lz"]][]

      et_grp <- it_et[it_grp_id != group_exl]
      lz_grp <- it_lz[it_grp_id != group_exl]
    }
  } else {
    it_grp_id <- iteration_snapshot[["group_id"]][]
    it_et <- iteration_snapshot[["et"]][]
    it_lz <- iteration_snapshot[["lz"]][]

    et_grp <- it_et[it_grp_id == group_id]
    lz_grp <- it_lz[it_grp_id == group_id]
  }

  xmin <- c(-bounds$lz_bound, -bounds$et_bound)
  xmax <- c(bounds$lz_bound, 0)
  gridsize <- c(100, 100)

  kde_grp <- kde(
    x = cbind(lz_grp, et_grp), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )

  if (identical(print_plot, TRUE)) {
    # Plot KDE and scatter plot
    par(mfrow = c(1, 1))

    # Set up an empty plot with appropriate limits
    plot(kde_grp,
      cont = c(75), xlab = "lz", ylab = "et", drawlabels = FALSE, col = "blue",
      xlim = c(-bounds$lz_bound, bounds$lz_bound),
      ylim = c(-bounds$et_bound, 0)
    )

    points(lz_grp, et_grp,
      pch = 19,
      col = adjustcolor("blue", alpha.f = 0.5), cex = 0.5
    )
  }

  return(kde_grp)
}

get_kl <- function(kde_p, kde_q, bounds) {
  # Define the common grid
  xmin <- c(-bounds$lz_bound, -bounds$et_bound)
  xmax <- c(bounds$lz_bound, 0)
  gridsize <- c(100, 100)

  # Evaluate both KDEs on the common grid
  grid_points <- expand.grid(
    seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2])
  )

  # Evaluate the KDEs at grid points
  P <- predict(kde_p, x = grid_points)
  Q <- predict(kde_q, x = grid_points)

  # Normalize densities
  # P <- P / sum(P)
  # Q <- Q / sum(Q)

  # print(sum(P))
  # print(sum(Q))

  # Compute KL divergence (avoiding division by zero)
  mask <- (P > 0) & (Q > 0)

  # Normalize the masked P and Q
  P_masked <- P[mask]
  Q_masked <- Q[mask]

  # Normalize both masked distributions to sum to 1
  P_masked <- P_masked / sum(P_masked)
  Q_masked <- Q_masked / sum(Q_masked)

  print(sum(P_masked))
  print(sum(Q_masked))

  kl_pq <- sum(P_masked * log(P_masked / Q_masked))

  return(kl_pq)
}

# SET UP ###########################################################################

# Specify the path to your HDF5 file
file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# Open the HDF5 file
h5file <- H5File$new(file_path, mode = "r")

it_lst <- seq(0, 100, by = 1)
snapshot <- 600

h_bandwidth_list <- get_all_h_bandwidths(h5file, it_lst, snapshot)
average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

bounds <- get_grid_bounds(h5file, it_lst, snapshot)

# group_id <- FALSE
# group_id <- 8580896
it <- 0
snapshot <- 600
print_plot <- FALSE

# group_exl <- FALSE
# group_exl <- 8580896

group_interest <- 8580896

# TESTING ###########################################################################

kde_rest <- get_initial_kde(
  h5file, it, snapshot, bounds,
  group_id = FALSE, print_plot = FALSE, group_exl = group_interest
)

kde_interest <- get_initial_kde(
  h5file, it, snapshot, bounds,
  group_id = group_interest, print_plot = FALSE, group_exl = FALSE
)


kl_pq <- get_kl(kde_interest, kde_rest, bounds)

print(kl_pq)


### MUCK ##################################

# it_id <- get_it_id(it)
# snap_id <- get_snap_id(snapshot)

# iteration <- h5file[[it_id]]
# iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

# et <- iteration_snapshot[["et"]][]
# lz <- iteration_snapshot[["lz"]][]

# par(mfrow = c(1, 1))

# plot(kde_rest,
#   cont = c(75), xlab = "lz", ylab = "et", drawlabels = FALSE, col = "red",
#   xlim = c(-bounds$lz_bound, bounds$lz_bound),
#   ylim = c(-bounds$et_bound, 0)
# )

# plot(kde_interest,
#   cont = c(75), xlab = "lz", ylab = "et", drawlabels = FALSE, col = "green",
#   add = TRUE
# )

# points(lz, et,
#   pch = 19,
#   col = adjustcolor("blue", alpha.f = 0.5), cex = 0.5
# )
