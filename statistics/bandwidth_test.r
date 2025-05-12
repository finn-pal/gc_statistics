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


# Specify the paths
sim_dir <- "/Users/z5114326/Documents/simulations"
data_dir <- "/Users/z5114326/Documents/GitHub/gc_statistics/data"

# Values to iterate over
it_lst <- seq(0, 100, by = 1)

cont_level <- 0.75

# get_time_dep(it_lst[1:2], cont_level, sim_dir, data_dir)

######################################################################

proc_path <- file.path(sim_dir, "m12i", "m12i_processed.hdf5")
snap_path <- file.path(sim_dir, "snapshot_times_public.txt")

# Open the HDF5 file
proc_data <- H5File$new(proc_path, mode = "r")

src_dat <- proc_data[[it_id]][["source"]]
# snap_zform <- src_dat[["snap_zform"]][]
# print(snap_zform)

# Read the file while skipping comment lines
snap_time_pub <- read.table(snap_path, comment.char = "#", header = FALSE)
colnames(snap_time_pub) <- c(
    "index", "scale_factor", "redshift",
    "time_Gyr", "lookback_time_Gyr", "time_width_Myr"
)

# Filter indices where time is greater than 46
snap_time_filtered <- snap_time_pub[snap_time_pub$index >= 46, ]
snap_lst <- snap_time_filtered$index

h_dict <- list()
for (snap in snap_lst) {
    snap_id <- get_snap_id(snap)

    h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snap)
    average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

    h_dict[[snap_id]] <- average_h_bandwidth
}

print(snap_id)
print(h_dict[[snap_id]])
