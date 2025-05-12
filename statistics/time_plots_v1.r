# Load necessary libraries
library(hdf5r)
library(jsonlite)
library(ks)
library("RColorBrewer")


get_it_id <- function(it) {
  sprintf("it%03d", it)
}

get_snap_id <- function(snapshot) {
  sprintf("snap%03d", snapshot)
}

get_et_lz_norm_cov <- function(h5file, it, snapshot) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  # Access the data
  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et_norm <- iteration_snapshot[["et_norm"]][]
  lz_norm <- iteration_snapshot[["lz_norm"]][]

  data_matrix <- cbind(lz_norm, et_norm)

  num_data <- nrow(data_matrix)
  cov_matrix <- cov(data_matrix)

  # Set off-diagonal elements to 0
  cov_matrix[lower.tri(cov_matrix)] <- 0
  cov_matrix[upper.tri(cov_matrix)] <- 0
  cov_matrix / sqrt(num_data)
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
  sum_matrix / num_matrices
}
##########################################################################

get_kde_norm <- function(h5file, it, snapshot, group_id, bandwidth) {
  it_id <- get_it_id(it)
  snap_id <- get_snap_id(snapshot)

  iteration <- h5file[[it_id]]
  iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

  et_full <- iteration_snapshot[["et_norm"]][]
  lz_full <- iteration_snapshot[["lz_norm"]][]

  xmin <- c(-1, -1)
  xmax <- c(1, 0)
  gridsize <- c(1000, 1000)

  it_grp_id <- iteration_snapshot[["group_id"]][]

  et_grp <- et_full[it_grp_id == group_id]
  lz_grp <- lz_full[it_grp_id == group_id]

  et_rst <- et_full[it_grp_id != group_id]
  lz_rst <- lz_full[it_grp_id != group_id]

  kde_res <- kde(
    x = cbind(lz_grp, et_grp), H = bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )

  list(
    group_id = group_id,
    et_grp = et_grp, lz_grp = lz_grp,
    et_rst = et_rst, lz_rst = lz_rst,
    kde_res = kde_res
  )
}

get_contour_area <- function(
    kde_res, cont_level,
    xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
  # Extract grid spacing from kde output
  dx <- diff(kde_res$eval.points[[1]])[1]
  dy <- diff(kde_res$eval.points[[2]])[1]
  cell_area <- dx * dy

  den_matrix <- kde_res$estimate
  # Normalize density matrix to a proper probability distribution
  den_matrix <- den_matrix / sum(den_matrix * cell_area, na.rm = TRUE)

  # Flatten the density matrix and sort values (excluding NA)
  sort_den <- sort(na.omit(as.vector(den_matrix)), decreasing = TRUE)

  # Compute cumulative probability
  cum_probs <- cumsum(sort_den * cell_area)

  # Find the threshold density level at 75% cumulative probability
  threshold_index <- which(cum_probs >= cont_level)[1]
  density_at_threshold <- sort_den[threshold_index]

  # Create a mask for all grid cells with density >= threshold
  mask_contour <- den_matrix >= density_at_threshold

  # Compute the area of the cells inside the contour
  cont_area <- sum(mask_contour * cell_area, na.rm = TRUE)
  full_area <- (xmax[1] - xmin[1]) * (xmax[2] - xmin[2])
  cont_area_per <- cont_area / full_area

  # Verify normalization
  # total_density <- sum(den_matrix * cell_area, na.rm = TRUE)
  # print(paste("Total density (should be approx 1):", total_density))

  return(list(
    cont_level = cont_level,
    cont_area = cont_area, cont_area_per = cont_area_per,
    den_matrix = den_matrix, density_at_threshold = density_at_threshold
  ))
}

get_kde_from_vals <- function(
    x, y, bandwidth,
    xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
  # return kde
  kde(
    x = cbind(x, y), H = bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
  )
}

get_kl <- function(
    kde_p, kde_q,
    xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
  # Evaluate both KDEs on the common grid
  grid_points <- expand.grid(
    seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2])
  )

  # Evaluate the KDEs at grid points
  p <- predict(kde_p, x = grid_points)
  q <- predict(kde_q, x = grid_points)

  # Compute KL divergence (avoiding division by zero)
  mask <- (p > 0) & (q > 0)

  # Normalize the masked P and Q
  p_masked <- p[mask]
  q_masked <- q[mask]

  # Normalize both masked distributions to sum to 1
  p_masked <- p_masked / sum(p_masked)
  q_masked <- q_masked / sum(q_masked)

  # print(sum(P_masked))
  # print(sum(Q_masked))

  sum(p_masked * log(q_masked / p_masked))
}

append_to_list <- function(lst, key, value) {
  lst[[as.character(key)]] <- append(lst[[as.character(key)]], list(value))
}

##########################################################################

# # Specify the path to your HDF5 file
# file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# # Open the HDF5 file
# proc_data <- H5File$new(file_path, mode = "r")

# h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snapshot)
# average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

# it_lst <- seq(0, 100, by = 1)

# group_id <- 8580896
# # group_id <- 0
# snapshot <- 600
# it <- 0

# kde_res_lst <- get_kde_norm(proc_data, it, snapshot, group_id, average_h_bandwidth)
# kde_res <- kde_res_lst$kde_res

# cont_res <- get_contour_area(kde_res_lst$kde_res, 0.75)

# contour(kde_res$eval.points[[1]], kde_res$eval.points[[2]], cont_res$den_matrix,
#   levels = cont_res$density_at_threshold, drawlabels = FALSE,
#   col = "blue", lwd = 2,
#   xlab = "lz_norm", ylab = "et_norm"
# )

# points(kde_res_lst$lz_grp, kde_res_lst$et_grp, col = "blue", pch = 16, cex = 0.5)
# points(kde_res_lst$lz_rst, kde_res_lst$et_rst, col = "red", pch = 16, cex = 0.5)

##########################################################################

# sim <- "m12i"

# it <- 0
# it_id <- get_it_id(it)

# data_dir <- "/Users/z5114326/Documents/GitHub/gc_statistics/data/"
# kin_data <- file.path(
#   data_dir, "kin_dicts", sim,
#   paste0(it_id, "_kin_dict.json")
# )

# # Load the JSON file into an R list
# kin_dict <- fromJSON(kin_data)

# print(names(kin_dict))

##########################################################################

# Specify the path to your HDF5 file
file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# Open the HDF5 file
proc_data <- H5File$new(file_path, mode = "r")

h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snapshot)
average_h_bandwidth <- get_average_matrix(h_bandwidth_list)



it_lst <- seq(0, 100, by = 1)

it <- 0
it_id <- get_it_id(it)
src_dat <- proc_data[[it_id]][["source"]]
it_grp_id <- src_dat[["group_id"]][]
group_ids <- unique(it_grp_id[it_grp_id >= 0])

it_dict <- list() # Create an empty list to store the data

i <- 1
for (group_id in group_ids[1:2]) { # Iterating through all group_ids
  cat(i, "/", length(group_ids), "\n")
  i <- i + 1

  # Initialize an empty list for each group_id in the it_dict
  group_dict <- it_dict[[as.character(group_id)]]
  group_dict <- list()

  src_mask <- ((src_dat[["group_id"]][] == group_id) & (src_dat[["analyse_flag"]][] == 1)) # nolint

  if (group_id == 0) {
    # if not accreted (formed in-situ)
    snap_base <- min(src_dat[["snap_zform"]][src_mask], na.rm = TRUE)
  } else {
    snap_base <- min(src_dat[["snap_acc"]][src_mask], na.rm = TRUE)
  }

  key_lst <- c()
  for (key in names(proc_data[[it_id]][["snapshots"]])) {
    snap_grp_lst <- proc_data[[it_id]][["snapshots"]][[key]][["group_id"]][]
    if (group_id %in% snap_grp_lst) {
      key_num <- as.integer(substr(key, 5, nchar(key)))
      if (!is.na(key_num) && key_num >= snap_base && key_num >= 46) {
        key_lst <- c(key_lst, key)
      }
    }
  }

  # Create an empty data frame for storing results
  result_df <- data.frame(
    key_num = integer(),
    den_thresh = numeric(),
    cont_area_per = numeric(),
    kl = numeric(),
    med_lz_norm_grp = numeric(),
    med_et_norm_grp = numeric(),
    avg_lz_norm_grp = numeric(),
    avg_et_norm_grp = numeric(),
    avg_vx_grp = numeric(),
    avg_vy_grp = numeric(),
    avg_vz_grp = numeric(),
    sig_vx_grp = numeric(),
    sig_vy_grp = numeric(),
    sig_vz_grp = numeric(),
    num_gc = integer(),
    stringsAsFactors = FALSE
  )

  # Adjust iteration here as needed (1:2 is just for testing)
  for (key in key_lst[1:2]) {
    # get snap value
    key_num <- as.integer(substr(key, 5, nchar(key)))
    print(key_num)

    snap_dat <- proc_data[[it_id]][["snapshots"]][[key]]

    snap_grp_lst <- snap_dat[["group_id"]][]

    lz_norm <- snap_dat[["lz_norm"]][]
    et_norm <- snap_dat[["et_norm"]][]

    vx <- snap_dat[["vx"]][]
    vy <- snap_dat[["vy"]][]
    vz <- snap_dat[["vz"]][]

    lz_norm_grp <- lz_norm[snap_grp_lst == group_id]
    et_norm_grp <- et_norm[snap_grp_lst == group_id]

    lz_norm_rst <- lz_norm[snap_grp_lst != group_id]
    et_norm_rst <- et_norm[snap_grp_lst != group_id]

    vx_grp <- vx[snap_grp_lst == group_id]
    vy_grp <- vy[snap_grp_lst == group_id]
    vz_grp <- vz[snap_grp_lst == group_id]

    kde_p <- get_kde_from_vals(lz_norm_grp, et_norm_grp, average_h_bandwidth)
    kde_q <- get_kde_from_vals(lz_norm_rst, et_norm_rst, average_h_bandwidth)

    cont_res <- get_contour_area(kde_p, 0.75)

    # Values to return
    den_thresh <- cont_res$density_at_threshold
    cont_area_per <- cont_res$cont_area_per

    # Create a temporary row to append
    temp_row <- data.frame(
      key_num = key_num,
      den_thresh = den_thresh,
      cont_area_per = cont_area_per,
      kl = get_kl(kde_p, kde_q),
      med_lz_norm_grp = median(lz_norm_grp, na.rm = TRUE),
      med_et_norm_grp = median(et_norm_grp, na.rm = TRUE),
      avg_lz_norm_grp = mean(lz_norm_grp, na.rm = TRUE),
      avg_et_norm_grp = mean(et_norm_grp, na.rm = TRUE),
      avg_vx_grp = mean(vx_grp, na.rm = TRUE),
      avg_vy_grp = mean(vy_grp, na.rm = TRUE),
      avg_vz_grp = mean(vz_grp, na.rm = TRUE),
      sig_vx_grp = sd(vx_grp, na.rm = TRUE),
      sig_vy_grp = sd(vy_grp, na.rm = TRUE),
      sig_vz_grp = sd(vz_grp, na.rm = TRUE),
      num_gc = length(lz_norm_grp)
    )

    # Append the result to the data frame
    result_df <- rbind(result_df, temp_row)
  }

  # After iterating over the keys, append the result_df to the list for the group_id
  group_dict <- append_to_list(group_dict, group_id, result_df)
}

print(kl)

# print(key_lst)

# et_full <- iteration_snapshot[["et_norm"]][]
# lz_full <- iteration_snapshot[["lz_norm"]][]

# xmin <- c(-1, -1)
# xmax <- c(1, 0)
# gridsize <- c(1000, 1000)

# it_grp_id <- iteration_snapshot[["group_id"]][]

#  et_grp <- et_full[it_grp_id == group_id]
