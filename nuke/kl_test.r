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

    # et_norm <- iteration_snapshot[["et"]][]
    # lz_norm <- iteration_snapshot[["lz"]][]

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

snapshot <- 600
snap_id <- get_snap_id(snapshot)

src_dat <- proc_data[[it_id]][["source"]]
it_grp_id <- src_dat[["group_id"]][]

group_id <- 0

####################################################################

xmin <- c(-1, -1)
xmax <- c(1, 0)
gridsize <- c(500, 500)

snap_dat <- proc_data[[it_id]][["snapshots"]][[snap_id]]

et_full <- snap_dat[["et_norm"]][]
lz_full <- snap_dat[["lz_norm"]][]

####################################################################

it_grp_id <- snap_dat[["group_id"]][]

et_grp <- et_full[it_grp_id == group_id]
lz_grp <- lz_full[it_grp_id == group_id]

et_rst <- et_full[it_grp_id != group_id]
lz_rst <- lz_full[it_grp_id != group_id]

kde_p <- kde(
    x = cbind(lz_grp, et_grp), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

kde_q <- kde(
    x = cbind(lz_rst, et_rst), H = average_h_bandwidth,
    xmin = xmin, xmax = xmax, gridsize = gridsize
)

plot(kde_p, cont = 75)
plot(kde_q, cont = 75, add = TRUE)

####################################################################

grid_points <- expand.grid(
    seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2])
)

p <- predict(kde_p, x = grid_points)
q <- predict(kde_q, x = grid_points)

mask <- (p > 0) & (q > 0)

p_masked <- p[mask]
q_masked <- q[mask]

p_masked <- p_masked / sum(p_masked)
q_masked <- q_masked / sum(q_masked)

print(sum(p_masked * log(p_masked / q_masked)))

####################################################################

get_kl <- function(kde_p, kde_q, bounds) {
    # Define the common grid
    # xmin <- c(-bounds$lz_bound, -bounds$et_bound)
    # xmax <- c(bounds$lz_bound, 0)
    # gridsize <- c(100, 100)

    xmin <- c(-1, -1)
    xmax <- c(1, 0)
    gridsize <- c(500, 500)

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

print(get_kl(kde_p, kde_q, bounds))
