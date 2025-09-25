col <- "darkgreen" # the color you want

# Save to PNG file
png("point_plot.png", width = 600, height = 600)

# First point
plot(0, 0,
    col = col, # use your chosen color
    pch = 19, # filled circle
    cex = 3,
    xlim = c(-1, 1),
    ylim = c(-1, 1),
    xlab = "X",
    ylab = "Y",
    main = paste("Points in", col)
)

# Add second point
points(0, 1,
    col = "green", # use your chosen color
    pch = 19,
    cex = 3
)

# dev.off()
