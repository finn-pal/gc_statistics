library(ks) # For kernel density estimation

# Generate synthetic 2D data
set.seed(42)
n_points <- 500
points <- matrix(rbind(
    MASS::mvrnorm(n_points / 2, c(2, 2), matrix(c(1, 0.5, 0.5, 1), 2, 2)), # Overdense region
    MASS::mvrnorm(n_points / 2, c(-2, -2), matrix(c(1, -0.3, -0.3, 1), 2, 2)) # Background region
), ncol = 2)

# Kernel Density Estimation (KDE)
kde_result <- kde(points, H = diag(1, 2)) # Bandwidth matrix set to identity for simplicity

# Interpolate densities at data points
densities <- predict(kde_result, x = points)

# Function to compute KL divergence
kl_divergence <- function(p, q) {
    p <- p / sum(p) # Normalize to probability distributions
    q <- q / sum(q)
    nonzero <- p > 0 & q > 0 # Avoid log(0) issues
    sum(p[nonzero] * log(p[nonzero] / q[nonzero]))
}

# Function to minimize KL divergence over density threshold
optimize_threshold <- function(densities, points) {
    thresholds <- quantile(densities, seq(0.1, 0.9, by = 0.05)) # Candidate density thresholds
    kl_values <- sapply(thresholds, function(thresh) {
        d1 <- points[densities >= thresh, , drop = FALSE]
        d2 <- points[densities < thresh, , drop = FALSE]

        if (nrow(d1) < 5 || nrow(d2) < 5) {
            return(Inf)
        } # Avoid extreme cases

        # Compute KDE for both distributions
        kde1 <- kde(d1, H = diag(1, 2))
        kde2 <- kde(d2, H = diag(1, 2))

        # Estimate density values on shared grid
        grid_points <- points # Use original points as evaluation grid
        p <- predict(kde1, x = grid_points)
        q <- predict(kde2, x = grid_points)

        kl_divergence(p, q) # Compute KL divergence
    })

    # Find threshold that minimizes KL divergence
    best_thresh <- thresholds[which.min(kl_values)]
    return(list(threshold = best_thresh, kl_values = kl_values, thresholds = thresholds))
}

# Run optimization
result <- optimize_threshold(densities, points)
best_threshold <- result$threshold

# Plot results
plot(points, col = ifelse(densities >= best_threshold, "red", "blue"), pch = 20, main = "Optimal Overdensity Threshold")
legend("topright", legend = c("D1 (Overdense)", "D2 (Background)"), col = c("red", "blue"), pch = 20)

cat("Optimal Density Threshold:", best_threshold, "\n")
