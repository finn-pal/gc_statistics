library(ks)
library(sp)

# Generate sample data inside the circle
set.seed(42)
n <- 500
radius <- 1
theta <- runif(n, 0, 2 * pi)
r <- sqrt(runif(n, 0, radius^2)) # Ensures uniform sampling inside the circle
x <- r * cos(theta)
y <- r * sin(theta)
data <- cbind(x, y)

# Define KDE estimation
H <- Hpi(data) # Compute bandwidth matrix
kde_est <- kde(data, H = H, gridsize = c(100, 100), xmin = c(-radius, -radius), xmax = c(radius, radius))

# Extract density estimate
x_eval <- kde_est$eval.points[[1]] # x grid (sorted)
y_eval <- kde_est$eval.points[[2]] # y grid (sorted)
z_eval <- matrix(kde_est$estimate, nrow = length(x_eval), byrow = TRUE) # Reshape density values

par(mfrow = c(1, 1)) # Arrange plots in a 1x2 grid

# Plot the results
plot(data, pch = 16, cex = 0.5, col = "blue", xlab = "X", ylab = "Y", main = "2D KDE within a Circle")
# contour(x_eval, y_eval, z_eval, add = TRUE, levels = pretty(z_eval, 10)) # Adjust contour levels
symbols(0, 0, circles = radius, add = TRUE, lwd = 2, inches = FALSE) # Draw boundary

# Second plot: Contour with threshold level
# contour(x_eval, y_eval, z_eval)

plot(kde_est,
    cont = c(99), xlab = "lz", ylab = "et",
    drawlabels = FALSE, col = "red",
    main = "KDE Overlay of lz vs et", add = TRUE
)
