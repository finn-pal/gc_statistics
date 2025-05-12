# Load necessary libraries
library(hdf5r)
library(ks)

# Specify the path to your HDF5 file
file_path <- "/Users/z5114326/Documents/simulations/m12i/m12i_processed.hdf5"

# Open the HDF5 file
h5file <- H5File$new(file_path, mode = "r")

# iteration 0 ###################################################

# Access the data
it0 <- h5file[["it000"]]
it0_snapshots <- it0[["snapshots"]]
it0_600 <- it0_snapshots[["snap600"]]

# Access the 'group_id' dataset (assuming it exists in the group)
it0_group_id <- it0_600[["group_id"]][]
it0_acc_snap <- it0_600[["acc_snap"]][]
it0_et <- it0_600[["et"]][]
it0_lz <- it0_600[["lz"]][]

it0_et_0 <- it0_et[it0_group_id == 0]
it0_lz_0 <- it0_lz[it0_group_id == 0]

# Filter data by group_id == 0
it0_filtered_data <- data.frame(et = it0_et_0, lz = it0_lz_0)
it0_data_matrix <- cbind(it0_filtered_data$lz, it0_filtered_data$et)

it0_data_len <- nrow(it0_filtered_data)
it0_cov_matrix <- cov(it0_data_matrix)
it0_H_bandwidth <- it0_cov_matrix / sqrt(it0_data_len)

# Set off-diagonal elements to 0
# it0_H_bandwidth[upper.tri(it0_H_bandwidth)] <- 0
# it0_H_bandwidth[lower.tri(it0_H_bandwidth)] <- 0

it0_kde_est <- kde(x = cbind(it0_filtered_data$lz, it0_filtered_data$et), H = it0_H_bandwidth)

# iteration 1 ###################################################

it1 <- h5file[["it001"]]
it1_snapshots <- it1[["snapshots"]]
it1_600 <- it1_snapshots[["snap600"]]

# Access the 'group_id' dataset (assuming it exists in the group)
it1_group_id <- it1_600[["group_id"]][]
it1_acc_snap <- it1_600[["acc_snap"]][]
it1_et <- it1_600[["et"]][]
it1_lz <- it1_600[["lz"]][]

it1_et_0 <- it1_et[it1_group_id == 0]
it1_lz_0 <- it1_lz[it1_group_id == 0]

# Filter data by group_id == 0
it1_filtered_data <- data.frame(et = it1_et_0, lz = it1_lz_0)
it1_data_matrix <- cbind(it1_filtered_data$lz, it1_filtered_data$et)

it1_data_len <- nrow(it1_filtered_data)
it1_cov_matrix <- cov(it1_data_matrix)
it1_H_bandwidth <- it1_cov_matrix / sqrt(it1_data_len)

# Set off-diagonal elements to 0
# it1_H_bandwidth[upper.tri(it1_H_bandwidth)] <- 0
# it1_H_bandwidth[lower.tri(it1_H_bandwidth)] <- 0

it1_kde_est <- kde(x = cbind(it1_filtered_data$lz, it1_filtered_data$et), H = it1_H_bandwidth)


# Plotting ################################################

# Plot KDE and scatter plot
par(mfrow = c(1, 1))

# Set up an empty plot with appropriate limits
plot(it0_kde_est,
    cont = c(75), xlab = "lz", ylab = "et",
    drawlabels = FALSE, col = "red", main = "KDE Overlay of lz vs et"
)

# Add the second KDE estimate on top
plot(it1_kde_est, cont = c(75), col = "blue", add = TRUE, drawlabels = FALSE)


# Add scatter points from both iterations
points(it0_filtered_data$lz, it0_filtered_data$et, pch = 19, col = adjustcolor("red", alpha.f = 0.5), cex = 0.5)
points(it1_filtered_data$lz, it1_filtered_data$et, pch = 19, col = adjustcolor("blue", alpha.f = 0.5), cex = 0.5)

# Add a legend to differentiate KDEs
legend("topright",
    legend = c("Iteration 0", "Iteration 1"),
    col = c("red", "blue"), lwd = 2
)

# Close the HDF5 file
h5file$close()
