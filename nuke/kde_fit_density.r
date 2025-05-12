# Load necessary library
library(ks)

# Generate sample data
set.seed(123)
n <- 500
data1 <- cbind(rnorm(n, mean = 0, sd = 1), rnorm(n, mean = 0, sd = 1))
data2 <- cbind(rnorm(n, mean = 2, sd = 1), rnorm(n, mean = -2, sd = 1))

# Compute KDEs using a common grid
H1 <- Hpi(data1)
H2 <- Hpi(data2)

# Define a common grid range
xmin <- c(-4, -6)
xmax <- c(6, 4)
gridsize <- c(100, 100) # Ensures both KDEs use the same resolution

# Compute KDEs
kde1 <- kde(data1, H = H1, xmin = xmin, xmax = xmax, gridsize = gridsize)
kde2 <- kde(data2, H = H2, xmin = xmin, xmax = xmax, gridsize = gridsize)

# Convolve the two KDEs by summing their density estimates
convolved_density <- kde1$estimate + kde2$estimate

# Extract the grid points
x_grid <- kde1$eval.points[[1]]
y_grid <- kde1$eval.points[[2]]

# Find the contour level at 50% of the maximum density
max_density1 <- max(kde1$estimate)
contour_level_50_1 <- max_density1 * 0.75

max_density2 <- max(kde2$estimate)
contour_level_50_2 <- max_density2 * 0.75

max_density_convolved <- max(convolved_density)
contour_level_50_convolved <- max_density_convolved * 0.75

# Plot PDFs
par(mfrow = c(2, 2)) # Arrange plots in a 2x2 grid

# Plot first KDE PDF
image(x_grid, y_grid, kde1$estimate,
    col = terrain.colors(50), main = "KDE PDF 1", xlab = "X", ylab = "Y"
)
contour(x_grid, y_grid, kde1$estimate, levels = contour_level_50_1, add = TRUE, col = "red", drawlabels = FALSE)

# Plot second KDE PDF
image(x_grid, y_grid, kde2$estimate,
    col = terrain.colors(50), main = "KDE PDF 2", xlab = "X", ylab = "Y"
)
contour(x_grid, y_grid, kde2$estimate, levels = contour_level_50_2, add = TRUE, col = "red", drawlabels = FALSE)

# Plot convolved KDE PDF
image(x_grid, y_grid, convolved_density,
    col = terrain.colors(50), main = "Convolved KDE PDF", xlab = "X", ylab = "Y"
)
contour(x_grid, y_grid, convolved_density, levels = contour_level_50_convolved, add = TRUE, col = "red", drawlabels = FALSE)

# Reset plotting layout
par(mfrow = c(1, 1))
