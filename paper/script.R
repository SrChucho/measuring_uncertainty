library(vars)

# Load the data
gdp <- read.csv("Indicadores20260410102825.csv", row.names = 1)

# Convert the data to a time series object
gdp_ts <- ts(gdp, start = c(2003, 1), frequency = 4)

# make beautiful time series plot
plot(gdp_ts[,3], 
     main = "Mexico's GDP growth", xlab = "Year", ylab = "% change GDP", col = "blue", lwd = 2)

plot(gdp_ts[, 3],
     main     = "Mexico: Real GDP Growth",
     xlab     = "Year",
     ylab     = "Annual % Change",
     col      = "steelblue",
     lwd      = 2,
     type     = "l",
     bty      = "l",          # cleaner L-shaped axis box
     las      = 1,            # horizontal y-axis tick labels
     cex.main = 1.2,
     cex.lab  = 1.0,
     font.main = 1,           # plain (not bold) title
     yaxt     = "n")          # suppress default y-axis to customize it

# Custom y-axis
axis(2, las = 1)

# Zero reference line
abline(h = 0, col = "red", lty = 2, lwd = 1)


# Add grid
grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)

