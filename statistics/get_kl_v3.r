# load necessary libraries
library(akima)
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



# main #############################################################
# set up the parser
parser <- ArgumentParser()

parser$add_argument("-s", "--simulation", required = TRUE, type = "character", help = "simulation name (e.g. m12i)")
parser$add_argument("-a", "--iteration_low_limit", required = TRUE, type = "integer", help = "lower bound iteration")
parser$add_argument("-b", "--iteration_up_limit", required = TRUE, type = "integer", help = "upper bound iteration")
parser$add_argument("-c", "--cores", required = FALSE, type = "integer", help = "number of cores to use", default = 8)
parser$add_argument("-l", "--contour_level", required = FALSE, type = "numeric", help = "contour level", default = 0.75)
parser$add_argument("-p", "--snap_limit", required = FALSE, type = "integer", help = "minimum snapshot", default = 46)
parser$add_argument("-g", "--groups", required = FALSE, type = "integer", help = "groups", default = NULL, nargs = "+")
parser$add_argument("-o", "--plot", required = FALSE, type = "logical", help = "plot group contours", default = TRUE)

# parse the command-line arguments
args <- parser$parse_args()

# Access the arguments
sim <- args$simulation
it_min <- args$iteration_low_limit
it_max <- args$iteration_up_limit
cores <- args$cores
contour_level <- args$contour_level
snap_limit <- args$snap_limit
groups <- args$groups

sim_dir <- "/Users/z5114326/Documents/simulations"

save_dir <- file.path(sim_dir, sim, "contour_plots")
if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
}

it_lst <- seq(it_min, it_max, by = 1)

proc_path <- file.path(sim_dir, sim, paste0(sim, "_processed.hdf5"))
snap_path <- file.path(sim_dir, "snapshot_times_public.txt")
# Open the HDF5 file
proc_data <- H5File$new(proc_path, mode = "r")

# # pop <- get_bandwidth(it_lst, sim, sim_dir)
# # print(pop)

if (is.null(groups) || length(groups) == 0) {
    if (sim == "m12i") {
        groups <- list(0, 1920378, 4414957, 8580896, 13687078, 16199669, 19898495)
    } else {
        cat("need to add relevant groups\n")
    }
}

# Test ######################################

h_dict <- get_bandwidth(it_lst, sim, sim_dir)

it0 <- "it000"
it1 <- "it001"

group_id <- 0

snap_dat0 <- proc_data[[it0]][["snapshots"]][["snap600"]]
snap_grp0 <- abs(snap_dat0[["group_id"]][])
bnd_flag0 <- snap_dat0[["bound_flag"]][]
lz_norm0 <- snap_dat0[["lz_norm"]][]
et_norm0 <- snap_dat0[["et_norm"]][]
lz_norm_grp0 <- lz_norm0[(snap_grp0 == group_id) & (bnd_flag0 == 1)]
et_norm_grp0 <- et_norm0[(snap_grp0 == group_id) & (bnd_flag0 == 1)]

snap_dat1 <- proc_data[[it1]][["snapshots"]][["snap600"]]
snap_grp1 <- abs(snap_dat1[["group_id"]][])
bnd_flag1 <- snap_dat1[["bound_flag"]][]
lz_norm1 <- snap_dat1[["lz_norm"]][]
et_norm1 <- snap_dat1[["et_norm"]][]
lz_norm_grp1 <- lz_norm1[(snap_grp1 == group_id) & (bnd_flag1 == 1)]
et_norm_grp1 <- et_norm1[(snap_grp1 == group_id) & (bnd_flag1 == 1)]

bandwidth <- h_dict[["snap600"]]
kde0 <- get_kde_from_vals(lz_norm_grp0, et_norm_grp0, bandwidth)
kde1 <- get_kde_from_vals(lz_norm_grp1, et_norm_grp1, bandwidth)

kde_list <- list()
kde_list[[1]] <- kde0$estimate
kde_list[[2]] <- kde1$estimate

snap_id <- get_snap_id(600)

# print(kde1$cont)
cont <- 75
hts0 <- contourLevels(kde0, prob = (100 - cont) / 100, approx = TRUE)
hts1 <- contourLevels(kde1, prob = (100 - cont) / 100, approx = TRUE)

contour(kde0$eval.points[[1]],
    kde0$eval.points[[2]],
    kde0$estimate,
    levels = hts0,
    drawlabels = TRUE,
    labels = snap_id,
    # add = TRUE,
    col = "blue",
    lwd = 2
)

contour(kde1$eval.points[[1]],
    kde1$eval.points[[2]],
    kde1$estimate,
    levels = hts1,
    drawlabels = FALSE,
    labels = snap_id,
    add = TRUE,
    col = "blue",
    lwd = 2
)

hst_avg <- mean(c(hts0, hts1))

kde_array <- simplify2array(kde_list)
kde_est_avg <- apply(kde_array, c(1, 2), mean)

print(hst_avg)

contour(kde1$eval.points[[1]],
    kde1$eval.points[[2]],
    kde_est_avg,
    levels = hst_avg,
    drawlabels = FALSE,
    labels = snap_id,
    add = TRUE,
    col = "red",
    lwd = 2
)

grid_points <- expand.grid(
    x = kde1$eval.points[[1]],
    y = kde1$eval.points[[2]]
)

# 2. Evaluate kde0 at these grid points
dobs0 <- predict(kde0, x = grid_points)
dobs1 <- predict(kde1, x = grid_points)

# dobs_list <- list()
# dobs_list[[1]] <- dobs0
# dobs_list[[2]] <- dobs1

# dobs_array <- simplify2array(dobs_list)
# dobs_est_avg <- apply(dobs_array, c(1, 2), mean)

dobs_mat <- cbind(dobs0, dobs1)
dobs_est_avg <- rowMeans(dobs_mat)

# Reshape to matrix matching the grid
dobs_est_avg_matrix <- matrix(
    dobs_est_avg,
    nrow = length(kde1$eval.points[[1]]),
    ncol = length(kde1$eval.points[[2]])
)

hts_test <- quantile(kde0$estimate, prob = (100 - cont) / 100)

contour(kde1$eval.points[[1]],
    kde1$eval.points[[2]],
    kde_est_avg,
    levels = hts_test,
    drawlabels = FALSE,
    labels = snap_id,
    add = TRUE,
    col = "black",
    lwd = 1
)

# print(hts0)
# print(hts_test)

dobs <- kde0$estimate
hts <- quantile(dobs, prob = 0.75)

# contour(kde1$eval.points[[1]],
#     kde1$eval.points[[2]],
#     dobs,
#     levels = hts,
#     drawlabels = FALSE,
#     labels = snap_id,
#     add = TRUE,
#     col = "green",
#     lwd = 2
# )


hts_avg_correct <- quantile(kde_est_avg, prob = 0.25)

contour(kde1$eval.points[[1]],
    kde1$eval.points[[2]],
    kde_est_avg,
    levels = hts_avg_correct,
    drawlabels = FALSE,
    labels = snap_id,
    add = TRUE,
    col = "green",
    lwd = 2
)

x_all <- rbind(kde0$x, kde1$x)

avg_kde <- list(
    estimate = kde_est_avg,
    eval.points = kde0$eval.points,
    # x = kde0$x,
    x = x_all,
    H = kde0$H,
    # w = kde0$w,
    gridded = TRUE
)
class(avg_kde) <- "kde"

# Predict at original x positions
dobs_avg <- predict(avg_kde, x = avg_kde$x)

# Now compute the quantile threshold
hts_avg <- quantile(dobs_avg, prob = (100 - cont) / 100)
print(hts_avg)

contour(avg_kde$eval.points[[1]],
    avg_kde$eval.points[[2]],
    avg_kde$estimate,
    levels = hts_avg,
    drawlabels = FALSE,
    labels = snap_id,
    add = TRUE,
    col = "purple",
    lwd = 2
)
