# Load necessary library
library(ks)

# Generate sample data
set.seed(123)
n <- 500
data1 <- cbind(rnorm(n, mean = 0, sd = 1), rnorm(n, mean = 0, sd = 1))

# Compute KDEs using a common grid
H <- Hpi(data1) # Use a common bandwidth estimate
xmin <- c(-4, -6)
xmax <- c(6, 4)
gridsize <- c(100, 100)

kde1 <- kde(data1, H = H, xmin = xmin, xmax = xmax, gridsize = gridsize)

# Sum the KDEs (valid because they are on the same grid)
kde1_density <- kde1$estimate

# Normalize to ensure it's still a proper probability density function
dx <- diff(kde1$eval.points[[1]])[1]
dy <- diff(kde1$eval.points[[2]])[1]
cell_area <- dx * dy
kde1_density <- kde1_density / sum(kde1_density * cell_area) # Normalize

# Compute 75% probability contour level
sorted_densities <- sort(as.vector(kde1_density), decreasing = TRUE)
cumulative_probs <- cumsum(sorted_densities * cell_area)
threshold <- sorted_densities[which(cumulative_probs >= 0.75)[1]]

# Verify normalization: The sum should be close to 1
total_density <- sum(kde1_density * cell_area)
print(total_density) # Should be approximately 1

# Plot PDFs in a 1x2 grid (two different plots)
par(mfrow = c(1, 1)) # Arrange plots in a 1x2 grid

# First plot: KDE contour plot for 75% level
plot(kde1, cont = 75, col = "blue", drawlabels = FALSE, lwd = 2)

# Second plot: Contour with threshold level
contour(kde1$eval.points[[1]], kde1$eval.points[[2]], kde1_density,
    levels = threshold, drawlabels = FALSE, lwd = 1, col = "red", add = TRUE
)
