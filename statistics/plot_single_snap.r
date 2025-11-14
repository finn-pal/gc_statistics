# load necessary libraries
library(argparse)
library(hdf5r)
library(jsonlite)
library(ks)
# library("RColorBrewer")

# useful functions #################################################

get_it_id <- function(it) {
    sprintf("it%03d", it)
}

get_snap_id <- function(snapshot) {
    sprintf("snap%03d", snapshot)
}

# bandwidth functions ##############################################

get_cov <- function(h5file, it, snapshot, var_1 = "lz_norm", var_2 = "et_norm") {
    it_id <- get_it_id(it)
    snap_id <- get_snap_id(snapshot)

    # Access the data
    iteration <- h5file[[it_id]]
    iteration_snapshot <- iteration[["snapshots"]][[snap_id]]

    var_1_vals <- iteration_snapshot[[var_1]][]
    var_2_vals <- iteration_snapshot[[var_2]][]

    bound_flag <- iteration_snapshot[["bound_flag"]][]

    var_1_vals <- var_1_vals[bound_flag == 1]
    var_2_vals <- var_2_vals[bound_flag == 1]

    data_matrix <- cbind(var_1_vals, var_2_vals)

    num_data <- nrow(data_matrix)
    cov_matrix <- cov(data_matrix)

    # Set off-diagonal elements to 0
    cov_matrix[lower.tri(cov_matrix)] <- 0
    cov_matrix[upper.tri(cov_matrix)] <- 0
    cov_matrix / sqrt(num_data)
}

get_all_h_bandwidths <- function(h5file, it_lst, snapshot, var_1 = "lz_norm", var_2 = "et_norm") {
    lapply(it_lst, function(it) get_cov(h5file, it, snapshot, var_1, var_2))
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

    return(average_matrix) # nolint: return_linter.
}

get_bandwidth <- function(it_lst, sim, sim_dir, snap, var_1 = "lz_norm", var_2 = "et_norm") {
    proc_path <- file.path(sim_dir, sim, paste0(sim, "_processed.hdf5"))
    # Open the HDF5 file
    proc_data <- H5File$new(proc_path, mode = "r")

    cat("\n")
    cat("Retrieving bandwidths:", "\n")

    start_time <- proc.time()

    h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snap, var_1, var_2)
    average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

    elapsed_time <- proc.time() - start_time
    cat("Time to retrieve bandwidths:", elapsed_time["elapsed"], "sec\n")
    cat("\n")

    proc_data$close()

    return(average_h_bandwidth)
}

# main functions in determing kl ###################################

get_kde_from_vals <- function(xi, yi, bandwidth, xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
    kde(x = cbind(xi, yi), H = bandwidth, xmin = xmin, xmax = xmax, gridsize = gridsize)
}

process_snapshot_plot <- function(snap, group_ids, group_cols, it_lst, sim, sim_dir, bandwidth, cont = 75,
                                  var_1 = "lz_norm", var_2 = "et_norm",
                                  var_1_axis = expression(xi), var_2_axis = expression(epsilon)) {
    # Open the HDF5 file
    proc_path <- file.path(sim_dir, sim, paste0(sim, "_processed.hdf5"))
    proc_data <- H5File$new(proc_path, mode = "r")

    snap_id <- get_snap_id(snap)

    group_labels <- c("in-situ", letters[1:(length(group_ids) - 1)])

    cat("Starting snapshot", snap, "\n")

    first_group <- TRUE
    for (k in seq_along(group_ids)) {
        group_id <- group_ids[[k]]
        group_col <- group_cols[[k]]

        kde_list <- list()
        x_all_list <- list()

        for (j in seq_along(it_lst)) {
            it <- it_lst[[j]]
            it_id <- get_it_id(it)

            src_dat <- proc_data[[it_id]][["source"]]
            src_grp <- abs(src_dat[["group_id"]][])
            src_ana <- src_dat[["analyse_flag"]][]

            if (group_id == 0) {
                snap_base <- min(src_dat[["snap_zform"]][(src_grp == group_id) & (src_ana == 1)], na.rm = TRUE)
            } else {
                snap_base <- min(src_dat[["snap_acc"]][(src_grp == group_id) & (src_ana == 1)], na.rm = TRUE)
            }

            if (snap < snap_base) next

            snap_data <- proc_data[[it_id]][["snapshots"]][[snap_id]]
            snap_grp <- abs(snap_data[["group_id"]][])
            bnd_flag <- snap_data[["bound_flag"]][]
            acc_flag <- snap_data[["now_accreted"]][]
            var_1_vals <- snap_data[[var_1]][]
            var_2_vals <- snap_data[[var_2]][]
            var_1_grp_vals <- var_1_vals[(snap_grp == group_id) & (bnd_flag == 1) & (acc_flag == 1)]
            var_2_grp_vals <- var_2_vals[(snap_grp == group_id) & (bnd_flag == 1) & (acc_flag == 1)]

            if (length(var_1_grp_vals) < 2 || length(var_2_grp_vals) < 2) next
            if (length(var_1_grp_vals) != length(var_2_grp_vals)) next

            kde_j <- get_kde_from_vals(var_1_grp_vals, var_2_grp_vals, bandwidth)

            kde_list[[j]] <- kde_j$estimate
            x_all_list[[j]] <- kde_j$x
            if (j == 1) eval_points <- kde_j$eval.points
        }

        valid_kdes <- Filter(Negate(is.null), kde_list)
        if (length(valid_kdes) == 0) next

        kde_est_avg <- if (length(valid_kdes) == 1) {
            valid_kdes[[1]]
        } else {
            apply(simplify2array(valid_kdes), c(1, 2), mean)
        }

        avg_kde <- list(
            estimate = kde_est_avg,
            eval.points = eval_points,
            x = do.call(rbind, x_all_list),
            H = bandwidth,
            gridded = TRUE
        )
        class(avg_kde) <- "kde"

        dobs_avg <- predict(avg_kde, x = avg_kde$x)
        hts_avg <- quantile(dobs_avg, prob = (100 - cont) / 100)

        if (first_group) {
            plot.new()
            plot.window(xlim = c(-1, 1), ylim = c(-1, 0), xaxs = "i", yaxs = "i")
            axis(1)
            axis(2)
            box()
            first_group <- FALSE
        }

        contour(avg_kde$eval.points[[1]],
            avg_kde$eval.points[[2]],
            avg_kde$estimate,
            levels = hts_avg,
            drawlabels = FALSE,
            add = TRUE,
            col = group_col,
            lwd = 2
        )
    }

    title(xlab = var_1_axis, ylab = var_2_axis, cex.lab = 1.5)

    offset <- 0.05 * diff(par("usr")[3:4])
    x_range <- diff(par("usr")[1:2])
    for (g in seq_along(group_labels)) {
        if (g == 1) {
            x_pos <- par("usr")[1] - 0.025
        }
        if (g == 2) {
            x_pos <- x_pos + 0.125 * x_range
        } else {
            x_pos <- x_pos + 0.03 * x_range
        }
        text(
            x = x_pos,
            y = par("usr")[3] + offset,
            labels = group_labels[[g]],
            adj = c(0, 1), cex = 1.2, col = group_cols[[g]]
        )
    }

    proc_data$close()
}



# main #############################################################
# set up the parser
parser <- ArgumentParser()

parser$add_argument("-s", "--simulation", required = TRUE, type = "character", help = "simulation name (e.g. m12i)")
parser$add_argument("-a", "--iteration_low_limit", required = TRUE, type = "integer", help = "lower bound iteration")
parser$add_argument("-b", "--iteration_up_limit", required = TRUE, type = "integer", help = "upper bound iteration")
parser$add_argument("-c", "--cores", required = FALSE, type = "integer", help = "number of cores to use", default = 8)
parser$add_argument("-l", "--contour_level", required = FALSE, type = "numeric", help = "contour level", default = 75)
parser$add_argument("-p", "--snapshot", required = FALSE, type = "integer", help = "Snapshot", default = 600)

# parse the command-line arguments
args <- parser$parse_args()

# Access the arguments
sim <- args$simulation
it_min <- args$iteration_low_limit
it_max <- args$iteration_up_limit
cores <- args$cores
contour_level <- args$contour_level
snap <- args$snapshot

sim_dir <- "/Users/z5114326/Documents/simulations"

it_lst <- seq(it_min, it_max, by = 1)
it_df <- read.csv(file.path(sim_dir, "iteration_check.csv"))
it_msk <- it_df[[sim]] != 0 # sim is a column name stored in a variable
it_skip_arr <- it_df$it_id[it_msk]
it_lst <- it_lst[!(get_it_id(it_lst) %in% it_skip_arr)]
cat(it_lst)

# Groups ######################################

gc_file <- file.path(sim_dir, sim, "gc_groups.json")
gc_group_data <- fromJSON(gc_file)
group_map <- gc_group_data[[sim]]

# Extract the group IDs and colors
group_ids <- as.numeric(names(group_map))
group_cols <- as.character(unname(group_map[as.character(group_ids)]))

# Test ######################################

bandwidth <- get_bandwidth(it_lst, sim, sim_dir, snap)

process_snapshot_plot(
    snap = snap,
    group_ids = group_ids,
    group_cols = group_cols,
    it_lst = it_lst,
    sim = sim,
    sim_dir = sim_dir,
    bandwidth = bandwidth
)
