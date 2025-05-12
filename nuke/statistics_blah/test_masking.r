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

# TESTING ###############################################################

it_id <- get_it_id(it)
snap_id <- get_snap_id(snapshot)
group_id <- 8580896

iteration <- h5file[[it_id]]
iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

et <- iteration_snapshot[["et"]][]
lz <- iteration_snapshot[["lz"]][]

# Masking function (already defined)
mask_val <- function(x, y, a, b, c, d) {
    y >= a / ((abs(x) - c)^b) + d
}

# Mask parameters (already defined)
a <- -6997827.890258315
b <- 0.49889733064484165
c <- -363.0791864583225
d <- 24041.97805734597

xmin <- c(-bounds$lz_bound, -bounds$et_bound)
xmax <- c(bounds$lz_bound, 0)
gridsize <- c(1000, 1000)

# Create grid points for KDE
grid_x <- seq(xmin[1], xmax[1], length.out = gridsize[1])
grid_y <- seq(xmin[2], xmax[2], length.out = gridsize[2])

# Apply the mask to the grid
mask <- outer(
    grid_x, grid_y,
    Vectorize(function(x, y) mask_val(x, y, a, b, c, d))
)

# Set up KDE calculation using average bandwidth
kde_all <- kde(
    x = cbind(lz, et), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

# Reshape the KDE estimate into a matrix
density_matrix_all <- matrix(kde_all$estimate, nrow = gridsize[1], ncol = gridsize[2])
density_matrix_all[!mask] <- NA


##### Lets do this now ###########

density_trial <- FALSE

# group_interest <- 8580896
group_interest <- 0

grp_id <- iteration_snapshot[["group_id"]][]

et_grp <- it_et[it_grp_id == group_interest]
lz_grp <- it_lz[it_grp_id == group_interest]

et_rst <- it_et[it_grp_id != group_interest]
lz_rst <- it_lz[it_grp_id != group_interest]

kde_grp <- kde(
    x = cbind(lz_grp, et_grp), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

density_matrix_grp <- matrix(kde_grp$estimate, nrow = gridsize[1], ncol = gridsize[2])
density_matrix_grp[!mask] <- NA

kde_rst <- kde(
    x = cbind(lz_rst, et_rst), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

density_matrix_rst <- matrix(kde_rst$estimate, nrow = gridsize[1], ncol = gridsize[2])
density_matrix_rst[!mask] <- NA

###################################################################################
# Contour plot #############

# Extract grid spacing from kde output
dx <- diff(kde_grp$eval.points[[1]])[1]
dy <- diff(kde_grp$eval.points[[2]])[1]
cell_area <- dx * dy

# Normalize density matrix to a proper probability distribution
density_matrix_grp <- density_matrix_grp / sum(density_matrix_grp * cell_area, na.rm = TRUE)

# Flatten the density matrix and sort values (excluding NA)
sorted_densities <- sort(na.omit(as.vector(density_matrix_grp)), decreasing = TRUE)

# Compute cumulative probability
cumulative_probs <- cumsum(sorted_densities * cell_area)

# Find the threshold density level at 75% cumulative probability
threshold_index <- which(cumulative_probs >= 0.75)[1]
density_at_threshold <- sorted_densities[threshold_index]

# Verify normalization
total_density <- sum(density_matrix_grp * cell_area, na.rm = TRUE)
print(total_density) # Should be approximately 1

# Plot all three options on the same axis
par(mfrow = c(1, 1)) # 1x1 grid for a single plot


# Plot the KDE heatmap with the mask applied
color_palette <- colorRampPalette(brewer.pal(n = 8, name = "Blues"))(100)
image(grid_x, grid_y, density_matrix_all,
    col = color_palette, main = "Masked KDE Density Estimate",
    xlab = "lz", ylab = "et",
    xlim = c(-1.1 * bounds$lz_bound, 1.1 * bounds$lz_bound), # Larger xlim
    ylim = c(-1.1 * bounds$et_bound, 0.1 * bounds$et_bound) # Larger ylim
)

# Add the contour based on the density at 75% cumulative probability (red)
contour(kde_grp$eval.points[[1]], kde_grp$eval.points[[2]], density_matrix_grp,
    levels = density_at_threshold, drawlabels = FALSE, lwd = 1, col = "red", add = TRUE
)

# Plot the KDE contour plot for 75% level (blue)
plot(kde_grp, cont = 75, col = "blue", drawlabels = FALSE, lwd = 2, add = TRUE)
