# Load necessary libraries
library(hdf5r)
library(ks)
library("RColorBrewer")
library(ggplot2)
library(patchwork)

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

  # print(sum(P_masked))
  # print(sum(Q_masked))

  kl_pq <- sum(P_masked * log(P_masked / Q_masked))

  return(kl_pq)
}

get_kde <- function(h5file, it, snapshot, group_id, bounds, bandwidth) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et_full <- iteration_snapshot[["et"]][]
  lz_full <- iteration_snapshot[["lz"]][]

  xmin <- c(-bounds$lz_bound, -bounds$et_bound)
  xmax <- c(bounds$lz_bound, 0)
  gridsize <- c(1000, 1000)

  it_grp_id <- iteration_snapshot[["group_id"]][]

  et_grp <- et_full[it_grp_id == group_id]
  lz_grp <- lz_full[it_grp_id == group_id]

  # Step 1: Perform KDE for the group (no filtering of points at the start)
  kde_return <- kde(
    x = cbind(lz_grp, et_grp), H = bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )

  return(kde_return)
}

get_centroid_distances <- function(
    h5file, it, snapshot,
    group_a, group_b = 0, value = "et") {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  val_full <- iteration_snapshot[[value]][]

  it_grp_id <- iteration_snapshot[["group_id"]][]

  val_a <- val_full[it_grp_id == group_a]
  val_b <- val_full[it_grp_id == group_b]

  cent_a_val <- mean(val_a)
  cent_b_val <- mean(val_b)

  # distance <- sqrt((cent_b_val - cent_a_val)^2)
  distance <- (cent_a_val - cent_b_val)

  return(distance)
}

get_kde_split <- function(h5file, it, snapshot, group_id, bounds, bandwidth) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et_full <- iteration_snapshot[["et"]][]
  lz_full <- iteration_snapshot[["lz"]][]

  xmin <- c(-bounds$lz_bound, -bounds$et_bound)
  xmax <- c(bounds$lz_bound, 0)
  gridsize <- c(1000, 1000)

  it_grp_id <- iteration_snapshot[["group_id"]][]

  et_grp <- et_full[it_grp_id == group_id]
  lz_grp <- lz_full[it_grp_id == group_id]

  et_oth <- et_full[it_grp_id != group_id]
  lz_oth <- lz_full[it_grp_id != group_id]

  kde_p <- kde(
    x = cbind(lz_grp, et_grp), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )

  kde_q <- kde(
    x = cbind(lz_oth, et_oth), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )

  return(list(kde_p = kde_p, kde_q = kde_q))
}


# SET UP ##################################################################

# Specify the path to your HDF5 file
file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# Open the HDF5 file
h5file <- H5File$new(file_path, mode = "r")

it_lst <- seq(0, 100, by = 1)
snapshot <- 600

h_bandwidth_list <- get_all_h_bandwidths(h5file, it_lst, snapshot)
average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

bounds <- get_grid_bounds(h5file, it_lst, snapshot)


# it <- 56
# it_id <- get_it_id(it)
# snap_id <- get_snap_id(snapshot)

# iteration <- h5file[[it_id]]
# iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

# it_grp_id <- iteration_snapshot[["group_id"]][]

# # Remove duplicates
# it_grp_id_unique <- unique(it_grp_id)
# it_grp_id_unique_no_first <- it_grp_id_unique[-1]

# # Output the result
# print(it_grp_id_unique_no_first)

# # Count occurrences of each value in it_grp_id
# group_counts <- table(it_grp_id)

# # Initialize a vector to store the counts
# count_list <- numeric(length(it_grp_id_unique_no_first))
# lz_distance_list <- numeric(length(it_grp_id_unique_no_first))
# et_distance_list <- numeric(length(it_grp_id_unique_no_first))
# kl_list <- numeric(length(it_grp_id_unique_no_first))

# # Iterate over each value in it_grp_id_unique_no_first
# for (i in seq_along(it_grp_id_unique_no_first)) {
#   group_value <- it_grp_id_unique_no_first[i]
#   # Count occurrences of group_value in it_grp_id and store it in the count_vector
#   count_list[i] <- sum(it_grp_id == group_value)

#   kde_result <- get_kde_split(
#     h5file, it, snapshot,
#     group_value, bounds, average_h_bandwidth
#   )

#   kl_list[i] <- get_kl(kde_result$kde_p, kde_result$kde_q, bounds)

#   lz_distance_list[i] <- get_centroid_distances(h5file, it,
#     snapshot, group_value,
#     value = "lz"
#   )
#   et_distance_list[i] <- get_centroid_distances(h5file, it,
#     snapshot, group_value,
#     value = "et"
#   )
# }

# Output the count_vector
# print(count_list)
# print(distance_list)
# print(kl_list)

# df <- data.frame(
#   count = as.numeric(unlist(count_list)),
#   lz_distance = as.numeric(unlist(lz_distance_list)),
#   et_distance = as.numeric(unlist(et_distance_list)), # Corrected et_distance
#   kl = as.numeric(unlist(kl_list))
# )

# Print first few rows to check the data
# print(head(df))

# Define a common color scale
# color_scale <- scale_color_gradient(low = "blue", high = "red")

# # Plot 1: x = lz_distance
# plot1 <- ggplot(df, aes(x = lz_distance, y = count, color = kl)) +
#   geom_point(size = 3) +
#   color_scale + # Apply the common color scale
#   labs(
#     x = "Lz Distance From In-Situ Centroid",
#     y = "Number of GCs", color = "KL Divergence"
#   ) +
#   theme_minimal()

# # Plot 2: x = et_distance
# plot2 <- ggplot(df, aes(x = et_distance, y = count, color = kl)) +
#   geom_point(size = 3) +
#   color_scale + # Apply the same color scale
#   labs(x = "Et Distance From In-Situ Centroid", y = "Number of GCs", color = "KL Divergence") +
#   theme_minimal()

# # Combine both plots with a shared color bar
# combined_plot <- (plot1 + plot2) + plot_layout(guides = "collect")

# Print the combined plot
# print(combined_plot)



group_id <- 8580896
it_lst <- seq(0, 100, by = 1)
snapshot <- 600

kl_list <- list() # Initialize an empty list

for (it in it_lst) {
  print(it)
  kde_result <- get_kde_split(
    h5file, it, snapshot,
    group_id, bounds, average_h_bandwidth # Fixed group_id
  )

  kl_val <- get_kl(kde_result$kde_p, kde_result$kde_q, bounds)
  print(kl_val)
  # Append KL divergence value to kl_list
  kl_list[[length(kl_list) + 1]] <- kl_val
}

# Print the KL divergence list
# print(kl_list)

kl_vec <- as.numeric(unlist(kl_list))


# # Create a histogram
p <- ggplot(data.frame(kl = kl_vec), aes(x = kl)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(x = "KL Divergence", y = "Frequency", title = "Histogram of KL Divergence") +
  theme_minimal()

print(p)
