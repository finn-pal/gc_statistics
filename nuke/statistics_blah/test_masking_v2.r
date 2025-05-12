# Load necessary libraries
library(hdf5r)
library(ks)
library("RColorBrewer")

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

it <- 0
snapshot <- 600

cont_level <- 0.5

# TESTING ###############################################################

it_id <- get_it_id(it)
snap_id <- get_snap_id(snapshot)
# group_id <- 8580896
group_id <- 0

iteration <- h5file[[it_id]]
iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

et <- iteration_snapshot[["et"]][]
lz <- iteration_snapshot[["lz"]][]

xmin <- c(-bounds$lz_bound, -bounds$et_bound)
xmax <- c(bounds$lz_bound, 0)
gridsize <- c(1000, 1000)

# Set up KDE calculation using average bandwidth
kde_all <- kde(
    x = cbind(lz, et), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

# group_interest <- 8580896
group_interest <- 0

it_grp_id <- iteration_snapshot[["group_id"]][]

et_grp <- it_et[it_grp_id == group_interest]
lz_grp <- it_lz[it_grp_id == group_interest]

# Step 1: Perform KDE for the group (no filtering of points at the start)
kde_grp <- kde(
    x = cbind(lz_grp, et_grp), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

# Step 2: Normalize the density matrix to a proper probability distribution
dx <- diff(kde_grp$eval.points[[1]])[1]
dy <- diff(kde_grp$eval.points[[2]])[1]
cell_area <- dx * dy
density_matrix_grp <- kde_grp$estimate
density_matrix_grp <- density_matrix_grp / sum(density_matrix_grp * cell_area)

# Step 4: Flatten the density matrix and sort values (excluding NA)
density_values <- na.omit(as.vector(density_matrix_grp))
sorted_densities <- sort(density_values, decreasing = TRUE)

# Step 5: Compute cumulative probability
cumulative_probs <- cumsum(sorted_densities * cell_area)

# Step 6: Find the threshold density at 75% cumulative probability
threshold_index <- which(cumulative_probs >= cont_level)[1]
density_at_threshold <- sorted_densities[threshold_index]

# Verify normalization
total_density <- sum(density_matrix_grp * cell_area, na.rm = TRUE)
print(total_density) # Should be approximately 1


# Step 4: Map each point of lz and et to the grid
filtered_lz <- numeric(0)
filtered_et <- numeric(0)

points_within_contour <- logical(length(lz_grp)) # To store the results

for (i in 1:length(lz_grp)) {
    # Find the closest grid index for each lz and et point
    lz_idx <- which.min(abs(kde_grp$eval.points[[1]] - lz_grp[i])) # Find closest lz index
    et_idx <- which.min(abs(kde_grp$eval.points[[2]] - et_grp[i])) # Find closest et index

    # Correct the density matrix indexing
    if (density_matrix_grp[lz_idx, et_idx] >= density_at_threshold) {
        # Add to the filtered list if within the contour
        filtered_lz <- c(filtered_lz, lz_grp[i])
        filtered_et <- c(filtered_et, et_grp[i])
    }
}


###################################

# Step 7: Plot the original points vs the mapped points
# Plot the original points in blue
plot(lz_grp, et_grp,
    col = "red", pch = 16, cex = 0.5,
    xlab = "lz", ylab = "et",
    xlim = c(-1.1 * bounds$lz_bound, 1.1 * bounds$lz_bound),
    ylim = c(-1.1 * bounds$et_bound, 0.1 * bounds$et_bound)
)

# Plot the mapped points in red
points(filtered_lz, filtered_et, col = "blue", pch = 16, cex = 0.5)

contour(kde_grp$eval.points[[1]], kde_grp$eval.points[[2]], density_matrix_grp,
    levels = density_at_threshold, drawlabels = FALSE, lwd = 1, col = "red", add = TRUE
)
