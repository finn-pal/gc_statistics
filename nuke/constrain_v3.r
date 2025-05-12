# Load necessary libraries
library(ks)
library(entropy)

# Generate 2D sample data from different distributions
set.seed(123)
X <- matrix(rnorm(2000, mean = 0, sd = 1), ncol = 2) # N(0,1) in 2D
Y <- matrix(rnorm(2000, mean = 1, sd = 1), ncol = 2) # N(1,1) in 2D

# Perform KDE for both datasets
kde_X <- kde(X, H = diag(2) * 0.5) # Bandwidth matrix
kde_Y <- kde(Y, H = diag(2) * 0.5)

# Define a circular bounding region
radius <- 2 # Define circular boundary
inside_circle <- function(x, y, r) sqrt(x^2 + y^2) <= r

# Apply the mask to the evaluation grid
grid_points <- kde_X$eval.points # List of x and y evaluation points
x_grid <- grid_points[[1]]
y_grid <- grid_points[[2]]

# Create a mask for points inside the circle
mask <- outer(x_grid, y_grid, Vectorize(function(x, y) inside_circle(x, y, radius)))

# Extract density estimates within the circular region
P <- kde_X$estimate
Q <- kde_Y$estimate

# Normalize densities to sum to 1 within the bounded region
P[!mask] <- NA
Q[!mask] <- NA
P <- P / sum(P, na.rm = TRUE)
Q <- Q / sum(Q, na.rm = TRUE)

print(sum(P, na.rm = TRUE))
print(sum(Q, na.rm = TRUE))

# Avoid division by zero: replace zero values in Q with a small number
Q[Q == 0] <- 1e-10

# Compute KL divergence
kl_div <- sum(P[mask] * log(P[mask] / Q[mask]), na.rm = TRUE)
print(paste("KL Divergence:", kl_div))

# Compute entropy of P
entropy_P <- entropy(P[mask], unit = "log", na.rm = TRUE)
print(paste("Entropy of P:", entropy_P))

# Plot KDEs as heatmaps
par(mfrow = c(1, 2)) # Set up side-by-side plots

image(x_grid, y_grid, P, col = heat.colors(100), main = "Density Estimate P", xlab = "X", ylab = "Y")
contour(x_grid, y_grid, P, add = TRUE) # Add contours

image(x_grid, y_grid, Q, col = heat.colors(100), main = "Density Estimate Q", xlab = "X", ylab = "Y")
contour(x_grid, y_grid, Q, add = TRUE) # Add contours
