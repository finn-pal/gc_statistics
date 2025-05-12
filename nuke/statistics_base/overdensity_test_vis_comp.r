library(ks) # Kernel density estimation

# Generate synthetic 2D data
set.seed(42)
n_points <- 500
points <- matrix(rbind(
    MASS::mvrnorm(n_points / 2, c(2, 2), matrix(c(1, 0.5, 0.5, 1), 2, 2)), # Overdense region
    MASS::mvrnorm(n_points / 2, c(-2, -2), matrix(c(1, -0.3, -0.3, 1), 2, 2)) # Background
), ncol = 2)

# Kernel Density Estimation (KDE)
kde_result <- kde(points, H = diag(1, 2))

# Interpolate densities at data points
densities <- predict(kde_result, x = points)

# Function to compute KL divergence
kl_divergence <- function(p, q) {
    p <- p / sum(p) # Normalize distributions
    q <- q / sum(q)
    nonzero <- p > 0 & q > 0 # Avoid log(0) issues
    sum(p[nonzero] * log(p[nonzero] / q[nonzero]))
}

# Function to optimize KL divergence over density threshold
optimize_threshold <- function(densities, points) {
    thresholds <- quantile(densities, seq(0.1, 0.9, by = 0.05)) # Candidate density thresholds
    kl_values <- sapply(thresholds, function(thresh) {
        d1 <- points[densities >= thresh, , drop = FALSE]
        d2 <- points[densities < thresh, , drop = FALSE]

        if (nrow(d1) < 5 || nrow(d2) < 5) {
            return(Inf)
        } # Avoid extreme cases

        kde1 <- kde(d1, H = diag(1, 2))
        kde2 <- kde(d2, H = diag(1, 2))

        grid_points <- points # Use original points as evaluation grid
        p <- predict(kde1, x = grid_points)
        q <- predict(kde2, x = grid_points)

        kl_divergence(p, q) # Compute KL divergence
    })

    best_thresh <- thresholds[which.min(kl_values)]
    return(list(threshold = best_thresh, kl_values = kl_values, thresholds = thresholds))
}

# Run optimization
result <- optimize_threshold(densities, points)
best_threshold <- result$threshold

# Compute the mean density and threshold multiplier
mean_density <- mean(densities)
density_multiplier <- best_threshold / mean_density

cat("Optimal Density Threshold:", best_threshold, "\n")
cat("Mean Density:", mean_density, "\n")
cat("Density Multiplier:", density_multiplier, "x above mean density\n")

# Select some trial thresholds for comparison
trial_thresholds <- c(
    quantile(densities, 0.25), # Lower threshold
    best_threshold, # Optimal threshold
    quantile(densities, 0.75)
) # Higher threshold

# Plot different thresholds
par(mfrow = c(1, length(trial_thresholds))) # One row, three columns

for (thresh in trial_thresholds) {
    plot(points,
        col = ifelse(densities >= thresh, "red", "blue"), pch = 20,
        main = paste("Threshold =", round(thresh, 4))
    )
}

# Reset layout
par(mfrow = c(1, 1))
