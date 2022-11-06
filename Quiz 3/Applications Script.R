# ================================================
# Applications Script
# 
# STAT 460/560 - Statistical Theory I
# 
# Anthony-A. Christidis
# Department of Statistics
# University of British Columbia
# ================================================


# Source cubic spline function
source("Fit_Cubic_Spline.R")
source("Fit_Penalized_Cubic_Spline.R")
source("CV_Penalized_Cubic_Spline.R")

# Simulation of data 
x <- seq (pi, 5*pi, length = 500)
f <- function(x){
  return(sin(2*x)*log(x))
}
set.seed(0)
y <- f(x) + rnorm(length(x), mean = 0, sd = 1)

# ________________________________
# Cubic Splines - No Penalization
# ________________________________

# Spline fit
spline <- Fit_Cubic_Spline(x, y)

# Spline predictions
x_values <- seq(pi, 5*pi, length = 5000)
spline_predictions <- Cubic_Spline_Prediction(spline, x_values)

# Plot fitted cubic spline
plot(x, y, pch = 20)
lines(x_values, spline_predictions, lwd = 2, col = "red")
lines(x_values, f(x_values), lwd = 2, col = "blue")

# __________________________
# Cubic Splines - Penalized
# __________________________

# Spline fit
spline <- Fit_Penalized_Cubic_Spline(x, y, alpha = 0.8)

# Spline predictions
x_values <- seq(pi, 5*pi, length = 5000)
spline_predictions <- Cubic_Spline_Prediction(spline, x_values)

# Plot fitted cubic spline
plot(x, y, pch = 20)
lines(x_values, spline_predictions, lwd = 2, col = "red")
lines(x_values, f(x_values), lwd = 2, col = "blue")

# _____________________________
# Cubic Splines - CV Penalized
# _____________________________

# Spline fit
spline <- CV_Penalized_Cubic_Spline(x, y, 
                                    alpha_grid = seq(0.6, 0.95, by = 0.01), 
                                    n_folds = 5)

# Plot CV Error by alpha value
plot(spline$alpha_grid, spline$CV_prediction_error, pch = 20,
     xlab = "alpha", ylab = "CV Error")
abline(v = spline$alpha_opt)

# Spline predictions
x_values <- seq(pi, 5*pi, length = 5000)
spline_predictions <- Cubic_Spline_Prediction(spline, x_values)

# Plot fitted cubic spline
plot(x, y, pch = 20)
lines(x_values, spline_predictions, lwd = 2, col = "red")
lines(x_values, f(x_values), lwd = 2, col = "blue")

# _____________________________
# An interesting investigation
# _____________________________

# Spline with low alpha value
spline <- CV_Penalized_Cubic_Spline(x, y, 
                                    alpha_grid = seq(0.01, 0.1, by = 0.005), 
                                    n_folds = 5)
spline <- Fit_Penalized_Cubic_Spline(x, y, alpha = spline$alpha_opt)

# Plot the fit
plot(x, y, pch = 20)
spline_predictions <- Cubic_Spline_Prediction(spline, x)
lines(x, spline_predictions, lwd = 2, col = "red")
lines(x, f(x), lwd = 2, col = "blue")

# Plot the residuals from the fit
plot(x, y - Cubic_Spline_Prediction(spline, x), pch = 20,
     xlab = "x", ylab = "Residuals")


