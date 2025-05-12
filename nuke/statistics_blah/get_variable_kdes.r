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

groupings_by_density <- function(
    h5file, it_lst, snapshot, group_interest, cont_level,
    h_bandwidth_list, bounds) {
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

    et_grp <- et_full[it_grp_id == group_interest]
    lz_grp <- lz_full[it_grp_id == group_interest]

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

    # Step 4: Flatten the density matrix and sort values
    density_values <- na.omit(as.vector(density_matrix_grp))
    sorted_densities <- sort(density_values, decreasing = TRUE)

    # Step 5: Compute cumulative probability
    cumulative_probs <- cumsum(sorted_densities * cell_area)

    # Step 6: Find the threshold density at 75% cumulative probability
    threshold_index <- which(cumulative_probs >= cont_level)[1]
    density_at_threshold <- sorted_densities[threshold_index]

    # Verify normalization
    total_density <- sum(density_matrix_grp * cell_area)
    print(paste("Total density normalized to", total_density))

    # Step 4: Map each point of lz and et to the grid
    filtered_lz <- numeric(0)
    filtered_et <- numeric(0)

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

    # Step 1: Create an index for the filtered points
    filtered_indices <- which(lz_full %in% filtered_lz & et_full %in% filtered_et)

    # Step 2: Remove the filtered points from lz_full and et_full
    lz_remaining <- lz_full[-filtered_indices]
    et_remaining <- et_full[-filtered_indices]

    return(list(
        cont_level = cont_level,
        lz_filt = filtered_lz, et_filt = filtered_et,
        lz_rest = lz_remaining, et_rest = et_remaining
    ))
}

get_kde <- function(lz, et, bounds, bandwidth) {
    xmin <- c(-bounds$lz_bound, -bounds$et_bound)
    xmax <- c(bounds$lz_bound, 0)
    gridsize <- c(1000, 1000)

    # print(length(lz))
    # print(length(et))

    # Step 1: Perform KDE for the group (no filtering of points at the start)
    kde_return <- kde(
        x = cbind(lz, et), H = bandwidth,
        xmin = xmin, xmax = xmax, gridsize = gridsize
    )

    return(kde_return)
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


# group_interest <- 8580896
group_interest <- 19898495
cont_level <- 0.95

results <- groupings_by_density(
    h5file, it_lst, snapshot, group_interest, cont_level,
    h_bandwidth_list, bounds
)

kde_p <- get_kde(results$lz_filt, results$et_filt, bounds, average_h_bandwidth)
kde_q <- get_kde(results$lz_rest, results$et_rest, bounds, average_h_bandwidth)

kl_pq <- get_kl(kde_p, kde_q, bounds)

# print(kl_pq)


cont_level_lst <- seq(0.1, 0.95, by = 0.05)

kl_pq_list <- list() # Initialize an empty list to store kl_pq values

# Loop through each value in cont_level_lst
for (cont_level in cont_level_lst) {
    # Run the function with the current cont_level
    results <- groupings_by_density(
        h5file, it_lst, snapshot, group_interest, cont_level,
        h_bandwidth_list, bounds
    )

    # Calculate the KDEs
    kde_p <- get_kde(results$lz_filt, results$et_filt, bounds, average_h_bandwidth)
    kde_q <- get_kde(results$lz_rest, results$et_rest, bounds, average_h_bandwidth)

    # Calculate the KL divergence
    kl_pq <- get_kl(kde_p, kde_q, bounds)

    # Save the KL divergence value to the list
    kl_pq_list[[as.character(cont_level)]] <- kl_pq
}

# Convert the kl_pq_list into a numeric vector
kl_pq_values <- unlist(kl_pq_list)

# Create a plot
plot(
    cont_level_lst, # x-axis: contour levels
    kl_pq_values, # y-axis: KL divergence values
    type = "b", # type = "b" for both points and lines
    xlab = "Contour Level", # x-axis label
    ylab = "KL Divergence", # y-axis label
    main = "KL Divergence vs Contour Level", # plot title
    pch = 19, # point character, makes the points filled circles
    col = "blue", # color of the points and line
    lwd = 2 # line width
)
