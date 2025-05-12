library(ks)

et_bound <- 308100
lz_bound <- 22800

# xmin <- c(0, -et_bound)
xmin <- c(-lz_bound, -et_bound)
xmax <- c(lz_bound, 0)
gridsize <- c(1000, 1000)

# Create grid points
grid_points <- expand.grid(
    seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2])
)

# Create a matrix of 1s corresponding to each grid point
grid_values <- rep(1, nrow(grid_points))

# Create a matrix of values (for visualization purposes)
grid_matrix <- matrix(grid_values, nrow = gridsize[1], ncol = gridsize[2])

# Masking function
mask_val <- function(x, y, a, b, c, d) {
    y >= a / ((abs(x) - c)^b) + d
}

a <- -6997827.890258315
b <- 0.49889733064484165
c <- -363.0791864583225
d <- 24041.97805734597

# Create grid for masking
x_grid <- seq(xmin[1], xmax[1], length.out = gridsize[1])
y_grid <- seq(xmin[2], xmax[2], length.out = gridsize[2])

# Apply the mask to the grid
mask <- outer(
    x_grid, y_grid,
    Vectorize(function(x, y) mask_val(x, y, a, b, c, d))
)

# Apply the mask to grid_matrix by setting values to NA (or any other desired mask value)
grid_matrix[!mask] <- NA

# Plot the heatmap with the mask applied
image(seq(xmin[1], xmax[1], length.out = gridsize[1]),
    seq(xmin[2], xmax[2], length.out = gridsize[2]),
    grid_matrix,
    col = heat.colors(10), main = "Heatmap with Mask Applied",
    xlab = "lz", ylab = "et",
    xlim = c(-1.1 * lz_bound, 1.1 * lz_bound), # Larger xlim
    ylim = c(-1.1 * et_bound, 0.1 * et_bound) # Larger ylim
)
