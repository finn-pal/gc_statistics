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

  # print(sum(p_masked))
  # print(sum(q_masked))

  sum(p_masked * log(p_masked / q_masked))
}

append_to_list <- function(lst, key, value) {
  lst[[as.character(key)]] <- append(lst[[as.character(key)]], list(value))
}

##########################################################################

# Specify the paths
sim_dir <- "/Users/z5114326/Documents/simulations"
data_dir <- "/Users/z5114326/Documents/GitHub/gc_statistics/data"

proc_path <- file.path(sim_dir, "m12i", "m12i_processed.hdf5")
snap_path <- file.path(sim_dir, "snapshot_times_public.txt")

# Open the HDF5 file
proc_data <- H5File$new(proc_path, mode = "r")

h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snapshot)
average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

print(average_h_bandwidth)
print(average_h_bandwidth[1])
print(average_h_bandwidth[4])

# Read the file while skipping comment lines
snap_time_pub <- read.table(snap_path, comment.char = "#", header = FALSE)
colnames(snap_time_pub) <- c(
  "index", "scale_factor", "redshift",
  "time_Gyr", "lookback_time_Gyr", "time_width_Myr"
)

# get_time_dep <- function(it_lst, sim_dir, data_dir) {}

it_lst <- seq(0, 100, by = 1)

# Example list of it_ids
it_lst <- 1:2 # Adjust this to your actual list of it_ids

it_dict <- list() # Create an empty list to store the data

# Iterate through it_lst first
for (it in it_lst) {
  it_id <- get_it_id(it)
  src_dat <- proc_data[[it_id]][["source"]]
  it_grp_id <- src_dat[["group_id"]][]
  group_ids <- unique(it_grp_id[it_grp_id >= 0])

  # Initialize an entry for each it_id in it_dict
  it_dict[[it_id]] <- list()

  i <- 1
  for (group_id in group_ids[1:2]) { # Iterating through all group_ids for the current it_id
    cat(i, "/", length(group_ids), "\n")
    i <- i + 1

    # Initialize an empty list for each group_id in the it_dict under current it_id
    group_dict <- it_dict[[it_id]]
    group_dict[[as.character(group_id)]] <- list() # Initialize for each group_id

    src_mask <- ((src_dat[["group_id"]][] == group_id) & (src_dat[["analyse_flag"]][] == 1)) # nolint

    if (group_id == 0) {
      snap_base <- min(src_dat[["snap_zform"]][src_mask], na.rm = TRUE)
    } else {
      snap_base <- min(src_dat[["snap_acc"]][src_mask], na.rm = TRUE)
    }

    snap_id_lst <- c()
    for (snap_id in names(proc_data[[it_id]][["snapshots"]])) {
      snap_grp_lst <- proc_data[[it_id]][["snapshots"]][[snap_id]][["group_id"]][]
      if (group_id %in% snap_grp_lst) {
        snap <- as.integer(substr(snap_id, 5, nchar(snap_id)))
        # 46 is chosen as from this point AGAMA fitting works
        if (!is.na(snap) && snap >= snap_base && snap >= 46) {
          snap_id_lst <- c(snap_id_lst, snap_id)
        }
      }
    }

    # Create an empty data frame for storing results
    result_df <- data.frame(
      snapshot = integer(),
      time = numeric(),
      lbt = numeric(),
      redshift = numeric(),
      density_thresh = numeric(),
      contour_area_per = numeric(),
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
    for (snap_id in snap_id_lst[1:2]) {
      snap <- as.integer(substr(snap_id, 5, nchar(snap_id)))
      print(snap)

      snap_dat <- proc_data[[it_id]][["snapshots"]][[snap_id]]

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

      den_thresh <- cont_res$density_at_threshold
      cont_area_per <- cont_res$cont_area_per

      temp_row <- data.frame(
        snapshot = snap,
        time = snap_time_pub[snap_time_pub$index == snap, "time_Gyr"],
        lbt = snap_time_pub[snap_time_pub$index == snap, "lookback_time_Gyr"],
        redshift = snap_time_pub[snap_time_pub$index == snap, "redshift"],
        density_thresh = den_thresh,
        contour_area_per = cont_area_per,
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

      result_df <- rbind(result_df, temp_row)
    }

    # After iterating over the keys, append the result_df to the group_dict
    group_dict[[as.character(group_id)]] <- result_df

    # After finishing group_id iteration, update the it_dict with the group_dict for the current it_id
    it_dict[[it_id]] <- group_dict
  }
}

# Save the list to a JSON file
output_path <- file.path(data_dir, "m12i_time_dep.json")
write_json(it_dict, output_path, pretty = TRUE)

h5file$close()
