# Load necessary libraries
library(hdf5r)
library(jsonlite)
library(ks)
library("RColorBrewer")
library(parallel)

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

  bound_flag <- iteration_snapshot[["bound_flag"]][]

  et_norm <- et_norm[bound_flag == 1]
  lz_norm <- lz_norm[bound_flag == 1]

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
  sum_matrix <- matrix(0, nrow = nrow(h_bandwidth_list[[1]]), ncol = ncol(h_bandwidth_list[[1]]))
  count_matrix <- matrix(0, nrow = nrow(h_bandwidth_list[[1]]), ncol = ncol(h_bandwidth_list[[1]]))

  # Iterate through all matrices in the list
  for (i in 1:num_matrices) {
    matrix_i <- h_bandwidth_list[[i]]

    # Add values to sum_matrix only if they are not NA
    sum_matrix <- sum_matrix + ifelse(!is.na(matrix_i), matrix_i, 0)

    # Count the non-NA values for averaging
    count_matrix <- count_matrix + ifelse(!is.na(matrix_i), 1, 0)
  }

  # Compute the average matrix, dividing sum_matrix by count_matrix to account for NAs
  average_matrix <- sum_matrix / count_matrix

  # Handle the cases where count is zero (i.e., no valid values at that position)
  average_matrix[is.na(average_matrix)] <- 0

  return(average_matrix)
}
##########################################################################

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

  return(list(
    cont_level = cont_level,
    cont_area = cont_area, cont_area_per = cont_area_per,
    den_matrix = den_matrix, density_at_threshold = density_at_threshold
  ))
}

get_kde_from_vals <- function(
    x, y, bandwidth,
    xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
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

  sum(p_masked * log(p_masked / q_masked))
}

append_to_list <- function(lst, key, value) {
  lst[[as.character(key)]] <- append(lst[[as.character(key)]], list(value))
}

######################################################################

get_bandwidth <- function(
    it_lst, sim_dir, snap_lim = 46) {
  proc_path <- file.path(sim_dir, "m12i", "m12i_processed.hdf5")
  snap_path <- file.path(sim_dir, "snapshot_times_public.txt")
  # Open the HDF5 file
  proc_data <- H5File$new(proc_path, mode = "r")

  # Read the file while skipping comment lines
  snap_time_pub <- read.table(snap_path, comment.char = "#", header = FALSE)
  colnames(snap_time_pub) <- c(
    "index", "scale_factor", "redshift",
    "time_Gyr", "lookback_time_Gyr", "time_width_Myr"
  )

  # Filter indices where time is greater than 46
  snap_time_filtered <- snap_time_pub[snap_time_pub$index >= snap_lim, ]
  snap_lst <- snap_time_filtered$index

  cat("\n")
  cat("Retrieving bandwidths:", "\n")

  h_dict <- list()
  start_time <- proc.time()
  for (snap in snap_lst) {
    snap_id <- get_snap_id(snap)

    h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snap)
    average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

    h_dict[[snap_id]] <- average_h_bandwidth
  }
  elapsed_time <- proc.time() - start_time
  cat("Time to retrieve bandwidths:", elapsed_time["elapsed"], "sec\n")
  cat("\n")

  return(h_dict)
}

######################################################################

get_time_dep <- function(
    it, cont_level, h_dict,
    sim, sim_dir, snap_lim = 46) {
  # start_time <- proc.time()

  # proc_path <- file.path(sim_dir, sim, sim, "_processed.hdf5")
  proc_path <- file.path(sim_dir, sim, paste0(sim, "_processed.hdf5"))
  snap_path <- file.path(sim_dir, "snapshot_times_public.txt")

  # Open the HDF5 file
  proc_data <- H5File$new(proc_path, mode = "r")

  # Read the file while skipping comment lines
  snap_time_pub <- read.table(snap_path, comment.char = "#", header = FALSE)
  colnames(snap_time_pub) <- c(
    "index", "scale_factor", "redshift",
    "time_Gyr", "lookback_time_Gyr", "time_width_Myr"
  )

  it_dict <- list() # Create an empty list to store the data
  # Iterate through it_lst first

  it_id <- get_it_id(it)

  cat("Retrieving", it_id, "\n")

  src_dat <- proc_data[[it_id]][["source"]]
  it_grp_id <- src_dat[["group_id"]][]
  group_ids <- unique(abs(it_grp_id[it_grp_id != -2]))

  # Initialize an entry for each it_id in it_dict
  it_dict[[it_id]] <- list()

  # Initialize data structure for each group_id
  group_dict_template <- list(
    snapshot = numeric(),
    density_thresh = numeric(),
    contour_area_per = numeric(),
    kl = numeric(),
    del_lz_norm = numeric(),
    del_et_norm = numeric(),
    med_lz_norm_grp = numeric(),
    med_et_norm_grp = numeric(),
    avg_lz_norm_grp = numeric(),
    avg_et_norm_grp = numeric()
  )

  i <- 1
  # Iterating through all group_ids for the current it_id
  for (group_id in group_ids) {
    # Initialize an empty list for each group_id
    group_dict <- group_dict_template

    src_group_id <- src_dat[["group_id"]][]
    src_an_flag <- src_dat[["analyse_flag"]][]
    src_mask <- ((src_group_id == group_id) & (src_an_flag == 1))

    if (group_id == 0) {
      snap_base <- min(src_dat[["snap_zform"]][src_mask], na.rm = TRUE)
    } else {
      snap_base <- min(src_dat[["snap_acc"]][src_mask], na.rm = TRUE)
    }

    snap_id_lst <- c()
    for (snap_id in names(proc_data[[it_id]][["snapshots"]])) {
      grp_lst <- proc_data[[it_id]][["snapshots"]][[snap_id]][["group_id"]][]
      if (group_id %in% grp_lst) {
        snap <- as.integer(substr(snap_id, 5, nchar(snap_id)))
        # 46 is chosen as from this point AGAMA fitting works
        if (!is.na(snap) && snap >= snap_base && snap >= snap_lim) {
          snap_id_lst <- c(snap_id_lst, snap_id)
        }
      }
    }

    j <- 0
    for (snap_id in snap_id_lst) {
      j <- j + 1
      # cat(
      #   it_id, "- group count (", i, "/", length(group_ids),
      #   ") - snap count (", j, "/", length(snap_id_lst), ")\n"
      # )

      snap <- as.integer(substr(snap_id, 5, nchar(snap_id)))

      bandwidth <- h_dict[[snap_id]]

      snap_dat <- proc_data[[it_id]][["snapshots"]][[snap_id]]

      snap_grp_lst <- abs(snap_dat[["group_id"]][])
      bnd_flag <- snap_dat[["bound_flag"]][]

      lz_norm <- snap_dat[["lz_norm"]][]
      et_norm <- snap_dat[["et_norm"]][]

      lz_norm_grp <- lz_norm[(snap_grp_lst == group_id) & (bnd_flag == 1)]
      et_norm_grp <- et_norm[(snap_grp_lst == group_id) & (bnd_flag == 1)]
      lz_norm_rst <- lz_norm[(snap_grp_lst != group_id) & (bnd_flag == 1)]
      et_norm_rst <- et_norm[(snap_grp_lst != group_id) & (bnd_flag == 1)]

      kde_p <- get_kde_from_vals(
        lz_norm_grp, et_norm_grp, bandwidth
      )

      kde_q <- get_kde_from_vals(
        lz_norm_rst, et_norm_rst, bandwidth
      )

      cont_res <- get_contour_area(kde_p, cont_level)

      den_thresh <- cont_res$density_at_threshold
      cont_area_per <- cont_res$cont_area_per

      time_data <- snap_time_pub[snap_time_pub$index == snap, ]

      # Update the group_dict with results
      group_dict$snapshot <-
        c(group_dict$snapshot, snap)
      group_dict$time <-
        c(group_dict$time, time_data$time_Gyr)
      group_dict$lbt <-
        c(group_dict$lbt, time_data$lookback_time_Gyr)
      group_dict$redshift <-
        c(group_dict$redshift, time_data$redshift)
      group_dict$density_thresh <-
        c(group_dict$density_thresh, den_thresh)
      group_dict$contour_area_per <-
        c(group_dict$contour_area_per, cont_area_per)
      group_dict$kl <-
        c(group_dict$kl, get_kl(kde_p, kde_q))
      group_dict$del_lz_norm <-
        c(group_dict$del_lz_norm, h_dict[[snap_id]][1])
      group_dict$del_et_norm <-
        c(group_dict$del_et_norm, h_dict[[snap_id]][4])
      group_dict$med_lz_norm_grp <-
        c(group_dict$med_lz_norm_grp, median(lz_norm_grp, na.rm = TRUE))
      group_dict$med_et_norm_grp <-
        c(group_dict$med_et_norm_grp, median(et_norm_grp, na.rm = TRUE))
      group_dict$avg_lz_norm_grp <-
        c(group_dict$avg_lz_norm_grp, mean(lz_norm_grp, na.rm = TRUE))
      group_dict$avg_et_norm_grp <-
        c(group_dict$avg_et_norm_grp, mean(et_norm_grp, na.rm = TRUE))
    }

    # After iterating over the snapshots, update the group_dict in it_dict
    it_dict[[it_id]][[as.character(group_id)]] <- group_dict

    i <- i + 1
  }
  # elapsed_time <- proc.time() - start_time
  # cat("Time for", it_id, ":", elapsed_time["elapsed"], "sec\n")
  # cat("\n")

  # h5file$close()
  proc_data$close()
  return(it_dict)
}


##########################################################################

sim <- "m12i"

# Specify the paths
sim_dir <- "/Users/z5114326/Documents/simulations"
# data_dir <- "/Users/z5114326/Documents/GitHub/gc_statistics/data"
# save_dir <- "/Users/z5114326/Documents/simulations/" + sim
save_dir <- file.path(sim_dir, sim)

# Values to iterate over
it_lst <- seq(0, 1, by = 1)
# it_lst <- it_lst[1:8]

cont_level <- 0.75

h_dict <- get_bandwidth(it_lst, sim_dir)

# num_cores <- 8
num_cores <- 2

it_dict_list <- mclapply(it_lst, function(it) {
  get_time_dep(it, cont_level, h_dict, sim, sim_dir)
}, mc.cores = num_cores)

# Flatten the result to remove the extra level
it_dict_list <- unlist(it_dict_list, recursive = FALSE)

# Save the list to a JSON file
output_path <- file.path(save_dir, "kl_data.json")
write_json(it_dict_list, output_path, pretty = TRUE)

##########################################################################
