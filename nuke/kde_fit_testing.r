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

# Reshape the density matrix into a vector for further manipulation
convolved_density_vector <- as.vector(convolved_density)

# Sort the density values to find the cumulative density distribution
sorted_density <- sort(convolved_density_vector, decreasing = TRUE)

# Calculate the cumulative sum to find the 75% total probability
total_density <- sum(convolved_density_vector)
cumulative_density <- cumsum(sorted_density) / total_density

# Find the value that corresponds to 75% of the total density
contour_level <- sorted_density[which.min(abs(cumulative_density - 0.25))]

#
contour_75_1 <- contourLevels(kde1, prob = 0.75)
#

# Test

cont_test_vect <- as.vector(kde1$estimate)
sorted_density_test <- sort(cont_test_vect, decreasing = TRUE)
total_density_test <- sum(cont_test_vect)
cumulative_density_test <- cumsum(sorted_density_test) / total_density_test
contour_level_test <- sorted_density_test[which.min(abs(cumulative_density_test - 0.25))]


# Plot PDFs
par(mfrow = c(2, 2)) # Arrange plots in a 2x2 grid

# Plot first KDE PDF
image(x_grid, y_grid, kde1$estimate,
    col = terrain.colors(50), main = "KDE PDF 1", xlab = "X", ylab = "Y"
)
# contour(x_grid, y_grid, kde1$estimate, add = TRUE)

contour(x_grid, y_grid, kde1$estimate, levels = contour_75_1, drawlabels = FALSE, add = TRUE)
contour(x_grid, y_grid, kde1$estimate, levels = contour_level_test, add = TRUE, col = "blue", drawlabels = FALSE)

# Plot second KDE PDF
image(x_grid, y_grid, kde2$estimate,
    col = terrain.colors(50), main = "KDE PDF 2", xlab = "X", ylab = "Y"
)
# contour(x_grid, y_grid, kde2$estimate, add = TRUE)

# Plot convolved KDE PDF
image(x_grid, y_grid, convolved_density,
    col = terrain.colors(50), main = "Convolved KDE PDF", xlab = "X", ylab = "Y"
)
# contour(x_grid, y_grid, convolved_density, add = TRUE)

# Plot the 75% contour for the convolved KDE PDF
contour(x_grid, y_grid, convolved_density, levels = contour_level, add = TRUE, col = "blue", drawlabels = FALSE)

# Plot again to verify consistency


# Reset plotting layout
par(mfrow = c(1, 1))
