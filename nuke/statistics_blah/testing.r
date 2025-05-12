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

# TESTING #######

it_id <- get_it_id(it)
snap_id <- get_snap_id(snapshot)
group_id <- 8580896

iteration <- h5file[[it_id]]
iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

et <- iteration_snapshot[["et"]][]
lz <- iteration_snapshot[["lz"]][]

xmin <- c(-bounds$lz_bound, -bounds$et_bound)
xmax <- c(bounds$lz_bound, 0)
gridsize <- c(500, 500)

kde_all <- kde(
    x = cbind(lz, et), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

# Evaluate both KDEs on the common grid
grid_points <- expand.grid(
    seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2])
)

# Define grid points
grid_x <- seq(xmin[1], xmax[1], length.out = gridsize[1])
grid_y <- seq(xmin[2], xmax[2], length.out = gridsize[2])

# Reshape KDE estimate into a matrix
density_matrix <- matrix(kde_all$estimate, nrow = gridsize[1], ncol = gridsize[2])

# Choose a color palette from RColorBrewer
# color_palette <- colorRampPalette(brewer.pal(n = 8, name = "Blues"))(100)
# color_palette <- colorRampPalette(brewer.pal(n = 8, name = "Greens"))(100)
color_palette <- colorRampPalette(brewer.pal(n = 8, name = "Greys"))(100)

# Plot KDE heatmap
image(grid_x, grid_y, density_matrix,
    col = color_palette, main = "KDE Density Estimate",
    xlab = "lz", ylab = "et"
)

# Add contour lines on top
plot(kde_all,
    cont = c(75), drawlabels = FALSE, col = "blue", add = TRUE,
    xlim = c(-bounds$lz_bound, bounds$lz_bound),
    ylim = c(-bounds$et_bound, 0)
)

print(bounds)
