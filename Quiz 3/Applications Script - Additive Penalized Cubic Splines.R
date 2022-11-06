# ==================================================
# Applications Script - Additive Penalized Splinles
# 
# STAT 460/560 - Statistical Theory I
# 
# Anthony-A. Christidis
# Department of Statistics
# University of British Columbia
# ==================================================

# Required libraries
library(plotly)

# Source cubic spline function
source("Fit_Additive_Penalized_Cubic_Spline.R")
source("CV_Additive_Penalized_Cubic_Spline.R")

# __________________________________
# Generate Data From True Function 
# _________________________________

# Setting the seed
set.seed(0)

# Simulation of data 
x1 <- seq(pi, 2*pi, length = 400)[sample(1:400)]
x2 <- seq(3*pi, 4*pi, length = 400)[sample(1:400)]
f <- function(x1, x2){
  
  return(x1^2*cos(3*x1)/10 + x2^2/10)
}
sim_data <- data.frame(x1 = x1, 
                       x2 = x2)

# Simulate noisy data
sim_data$y <- f(sim_data$x1, sim_data$x2) + rnorm(nrow(sim_data), mean = 0, sd = 1/2)

# __________________________________
# Plot True Surface with Noisy Data 
# __________________________________


# Plot true surface with noisy data
true_figure <- plot_ly()
true_figure <- true_figure %>% layout(scene = list(
  xaxis = list(title = "x2"),
  yaxis = list(title = "x1"),
  zaxis = list(title = "y")))
true_figure <- true_figure %>% add_trace(x = sim_data$x1, y = sim_data$x2, 
                                         z = f(sim_data$x1, sim_data$x2),
                                         type = "mesh3d")
true_figure <- true_figure %>% add_markers(x = sim_data$x1, y = sim_data$x2, 
                           z = sim_data$y, 
                           size = 1)
true_figure

# __________________________________
# Additive Penalized Cubic Splines 
# __________________________________

# Spline fit
splines <- Fit_Additive_Penalized_Cubic_Spline(x = cbind(sim_data$x1, sim_data$x2), 
                                               y = sim_data$y, 
                                               alpha = 0.5, 
                                               tolerance = 1e-6)

# Spline predictions
x_values <- expand.grid(seq(pi, 2*pi, length = 50),
                        seq(3*pi, 4*pi, length = 50))
spline_predictions <- Additive_Cubic_Spline_Prediction(splines, x_values)

# Plot fitted surface with noisy data
true_figure <- plot_ly()
true_figure <- true_figure %>% add_trace(x = x_values[, 1], y = x_values[, 2], 
                                         z = spline_predictions,
                                         type = "mesh3d")
true_figure <- true_figure %>% add_markers(x = sim_data$x1, y = sim_data$x2, 
                                           z = sim_data$y, 
                                           size = 1)
true_figure

# _____________________________
# Cubic Splines - CV Additive
# _____________________________

# Spline fit
splines <- CV_Additive_Penalized_Cubic_Spline(x = cbind(sim_data$x1, sim_data$x2), 
                                              y = sim_data$y,
                                              alpha_grid = seq(0.75, 0.99, by = 0.01), 
                                              n_folds = 5,
                                              tolerance = 1e-6)

# Plot CV Error by alpha value
plot(splines$alpha_grid, splines$CV_prediction_error, pch = 20,
     xlab = "alpha", ylab = "CV Error")
abline(v = splines$alpha_opt)

# Spline fit
splines <- Fit_Additive_Penalized_Cubic_Spline(cbind(sim_data$x1, sim_data$x2), sim_data$y, alpha = splines$alpha_opt,
                                               tolerance = 1e-6)

# Spline predictions
x_values <- expand.grid(seq(pi, 2*pi, length = 50),
                        seq(3*pi, 4*pi, length = 50))
spline_predictions <- Additive_Cubic_Spline_Prediction(splines, x_values)

# Plot fitted surface with noisy data
fit_figure <- plot_ly()
fit_figure <- fit_figure %>% layout(scene = list(
                      xaxis = list(title = "x2"),
                      yaxis = list(title = "x1"),
                      zaxis = list(title = "y")))
fit_figure <- fit_figure %>% add_trace(x = x_values[, 1], y = x_values[, 2], 
                                         z = spline_predictions,
                                         type = "mesh3d")
fit_figure <- fit_figure %>% add_markers(x = sim_data$x1, y = sim_data$x2, 
                                           z = sim_data$y, 
                                           size = 1)
fit_figure

