# Load necessary library
library(ks)

# Generate sample data
set.seed(123)
n <- 500
data1 <- cbind(rnorm(n, mean = 0, sd = 1), rnorm(n, mean = 0, sd = 1))
data2 <- cbind(rnorm(n, mean = 2, sd = 1), rnorm(n, mean = -2, sd = 1))

# Compute KDEs using a common grid
H <- Hpi(rbind(data1, data2)) # Use a common bandwidth estimate
xmin <- c(-4, -6)
xmax <- c(6, 4)
gridsize <- c(100, 100)

kde1 <- kde(data1, H = H, xmin = xmin, xmax = xmax, gridsize = gridsize)
kde2 <- kde(data2, H = H, xmin = xmin, xmax = xmax, gridsize = gridsize)

# Sum the KDEs (valid because they are on the same grid)
convolved_density <- kde1$estimate + kde2$estimate

# Normalize to ensure it's still a proper probability density function
dx <- diff(kde1$eval.points[[1]])[1]
dy <- diff(kde1$eval.points[[2]])[1]
cell_area <- dx * dy
convolved_density <- convolved_density / sum(convolved_density * cell_area) # Normalize

# Compute 75% probability contour level
sorted_densities <- sort(as.vector(convolved_density), decreasing = TRUE)
cumulative_probs <- cumsum(sorted_densities * cell_area)
threshold <- sorted_densities[which(cumulative_probs >= 0.25)[1]]

# Verify normalization: The sum should be close to 1
total_density <- sum(convolved_density * cell_area)
print(total_density) # Should be approximately 1

# Plot heatmap (image) and overlay the 75% contour
image(kde1$eval.points[[1]], kde1$eval.points[[2]], convolved_density,
    col = terrain.colors(50), main = "75% Contour of Convolved Density",
    xlab = "X", ylab = "Y"
)

# Overlay the 75% contour
contour(kde1$eval.points[[1]], kde1$eval.points[[2]], convolved_density,
    levels = threshold, drawlabels = FALSE, lwd = 2, col = "red", add = TRUE
)
