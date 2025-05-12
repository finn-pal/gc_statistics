# load necessary libraries
library(akima)
library(animation)
library(argparse)
library(hdf5r)
library(jsonlite)
library(ks)
library("RColorBrewer")
library(parallel)

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

get_bandwidth <- function(it_lst, sim, sim_dir, snap_lim = 46, var_1 = "lz_norm", var_2 = "et_norm") {
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

    # Filter indices where time is greater than 46
    snap_time_filtered <- snap_time_pub[snap_time_pub$index >= snap_lim, ]
    snap_lst <- snap_time_filtered$index

    cat("\n")
    cat("Retrieving bandwidths:", "\n")

    h_dict <- list()
    start_time <- proc.time()
    for (snap in snap_lst) {
        snap_id <- get_snap_id(snap)

        h_bandwidth_list <- get_all_h_bandwidths(proc_data, it_lst, snap, var_1, var_2)
        average_h_bandwidth <- get_average_matrix(h_bandwidth_list)

        h_dict[[snap_id]] <- average_h_bandwidth
    }
    elapsed_time <- proc.time() - start_time
    cat("Time to retrieve bandwidths:", elapsed_time["elapsed"], "sec\n")
    cat("\n")

    proc_data$close()

    return(h_dict)
}

# main functions in determing kl ###################################

get_kde_from_vals <- function(xi, yi, bandwidth, xmin = c(-1, -1), xmax = c(1, 0), gridsize = c(1000, 1000)) {
    kde(x = cbind(xi, yi), H = bandwidth, xmin = xmin, xmax = xmax, gridsize = gridsize)
}

# main_kde <- function(it_lst, cont_level, h_dict, sim, sim_dir, snap_lim) {}

# plotting functions ###############################################

plot_contours <- function(group_id, it_lst, snap_lst, save_dir, sim, sim_dir, h_dict, cont = 75,
                          var_1 = "lz_norm", var_2 = "et_norm", var_1_axis = expression(epsilon), var_2_axis = "e") {
    # Open the HDF5 file
    proc_path <- file.path(sim_dir, sim, paste0(sim, "_processed.hdf5"))
    proc_data <- H5File$new(proc_path, mode = "r")

    cat(group_id, "\n")

    saveGIF(
        {
            for (i in seq_along(snap_lst)) {
                snap <- snap_lst[[i]]
                snap_id <- get_snap_id(snap)
                # cat(snap_id, "\n")

                time <- snap_pub_data[snap_pub_data$index == snap, "time_Gyr"]

                bandwidth <- h_dict[[snap_id]]

                kde_list <- list()
                x_all_list <- list()
                for (j in seq_along(it_lst)) {
                    it <- it_lst[[j]]
                    it_id <- get_it_id(it)

                    # need to ensure the group in question exists within the main galaxy at this snapshot
                    src_dat <- proc_data[[it_id]][["source"]]
                    src_grp <- abs(src_dat[["group_id"]][])
                    src_ana <- src_dat[["analyse_flag"]][]

                    if (group_id == 0) {
                        snap_base <- min(src_dat[["snap_zform"]][(src_grp == group_id) & (src_ana == 1)], na.rm = TRUE)
                    } else {
                        snap_base <- min(src_dat[["snap_acc"]][(src_grp == group_id) & (src_ana == 1)], na.rm = TRUE)
                    }

                    if (snap < snap_base) next
                    ########################################################################################

                    snap_data <- proc_data[[it_id]][["snapshots"]][[snap_id]]
                    snap_grp <- abs(snap_data[["group_id"]][])
                    bnd_flag <- snap_data[["bound_flag"]][]
                    acc_flag <- snap_data[["now_accreted"]][]
                    var_1_vals <- snap_data[[var_1]][]
                    var_2_vals <- snap_data[[var_2]][]
                    var_1_grp_vals <- var_1_vals[(snap_grp == group_id) & (bnd_flag == 1) & (acc_flag == 1)]
                    var_2_grp_vals <- var_2_vals[(snap_grp == group_id) & (bnd_flag == 1) & (acc_flag == 1)]

                    kde_j <- get_kde_from_vals(var_1_grp_vals, var_2_grp_vals, bandwidth)

                    kde_list[[j]] <- kde_j$estimate
                    x_all_list[[j]] <- kde_j$x

                    if (j == 1) {
                        eval_points <- kde_j$eval.points
                    }
                }

                # if no data then skip to next snapshot
                if (length(kde_list) == 0) next

                kde_array <- simplify2array(kde_list)
                kde_est_avg <- apply(kde_array, c(1, 2), mean)

                x_all <- do.call(rbind, x_all_list)

                avg_kde <- list(
                    estimate = kde_est_avg,
                    eval.points = eval_points,
                    x = x_all,
                    H = bandwidth,
                    gridded = TRUE
                )
                class(avg_kde) <- "kde"

                dobs_avg <- predict(avg_kde, x = avg_kde$x)
                hts_avg <- quantile(dobs_avg, prob = (100 - cont) / 100)

                plot.window(
                    xlim = range(avg_kde$eval.points[[1]]),
                    ylim = range(avg_kde$eval.points[[2]])
                )

                contour(avg_kde$eval.points[[1]],
                    avg_kde$eval.points[[2]],
                    avg_kde$estimate,
                    levels = hts_avg,
                    drawlabels = FALSE,
                    col = "blue",
                    lwd = 2
                )

                title(
                    main = paste("Group ID:", group_id),
                    xlab = var_1_axis,
                    ylab = var_2_axis,
                    cex.lab = 1.5
                )

                mtext(paste("Time:", sprintf("%.3f", time), "Gyr"),
                    side = 3, line = 1.5, adj = 1, cex = 1.2, col = "blue"
                )
            }
        },
        movie.name = file.path(save_dir, paste0("group_", group_id, "_contour.gif")),
        interval = 0.5,
        ani.width = 600,
        ani.height = 600
    )

    # close proc_data
    # proc_data$close()
}

# main #############################################################
# set up the parser
parser <- ArgumentParser()

parser$add_argument("-s", "--simulation", required = TRUE, type = "character", help = "simulation name (e.g. m12i)")
parser$add_argument("-a", "--iteration_low_limit", required = TRUE, type = "integer", help = "lower bound iteration")
parser$add_argument("-b", "--iteration_up_limit", required = TRUE, type = "integer", help = "upper bound iteration")
parser$add_argument("-c", "--cores", required = FALSE, type = "integer", help = "number of cores to use", default = 8)
parser$add_argument("-l", "--contour_level", required = FALSE, type = "numeric", help = "contour level", default = 75)
parser$add_argument("-p", "--snap_limit", required = FALSE, type = "integer", help = "minimum snapshot", default = 46)
parser$add_argument("-o", "--plot", required = FALSE, type = "logical", help = "plot group contours", default = TRUE)
parser$add_argument("-k", "--get_kl", required = FALSE, type = "logical", help = "get kl for groups", default = TRUE)

# parse the command-line arguments
args <- parser$parse_args()

# Access the arguments
sim <- args$simulation
it_min <- args$iteration_low_limit
it_max <- args$iteration_up_limit
cores <- args$cores
contour_level <- args$contour_level
snap_limit <- args$snap_limit
plot <- args$plot
get_kl <- args$get_kl

sim_dir <- "/Users/z5114326/Documents/simulations"

save_dir <- file.path(sim_dir, sim, "contour_plots")
if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
}

snap_pub_dir <- paste0(sim_dir, "/snapshot_times_public.txt")
snap_pub_data <- read.table(snap_pub_dir, comment.char = "#", header = FALSE)
colnames(snap_pub_data) <- c(
    "index",
    "scale_factor",
    "redshift",
    "time_Gyr",
    "lookback_time_Gyr",
    "time_width_Myr"
)
snap_lst <- snap_pub_data[snap_pub_data$index >= snap_limit, "index"]

it_lst <- seq(it_min, it_max, by = 1)

##############################################

if (sim == "m12i") {
    groups <- list(0, 1920378, 4414957, 8580896, 13687078, 16199669, 19898495)
} else {
    cat("need to add relevant groups\n")
}

##############################


# Test ######################################
# snap <- 600
# group_id <- groups[[2]]
# cont <- 75

# var_1 <- "lz_norm"
# var_2 <- "et_norm"

# var_1_axis <- expression(epsilon)
# var_2_axis <- "e"

h_dict <- get_bandwidth(it_lst, sim, sim_dir)

mclapply(groups, function(group_id) {
    plot_contours(group_id, it_lst, snap_lst, save_dir, sim, sim_dir, h_dict, contour_level)
}, mc.cores = cores)

# print(it_lst)
# plot_contours(group_id, it_lst, snap_lst, save_dir, sim, sim_dir, h_dict, contour_level)
