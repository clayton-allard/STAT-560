# ================================================
# Implementation - CV Penalized Cublic Splines
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

# Fitting penalized cubic spline function 
CV_Penalized_Cubic_Spline <- function(x, y, 
                                      alpha_grid = seq(0.5, 0.99, by = 0.01), 
                                      n_folds = 5){
  
  # Check alpha grid
  if(any(alpha_grid <= 0) | any(alpha_grid >= 1))
    stop("The grid for alpha must only contain values greater than 0 and less than 1.")
  
  # Prediction Error
  CV_prediction_error <- numeric(length(alpha_grid))
  
  # Creating folds
  folds <- caret::createFolds(y = 1:length(y), k = n_folds)
  
  # Looping over tuning parameter
  for(alpha in alpha_grid){
    
    # Print iteration
    cat("CV for alpha = ", alpha, "\n", sep = "")
    
    # Looping over folds
    for(fold_id in 1:n_folds){
      
      # Subsetting training and test data
      x_train <- x[-folds[[fold_id]]]; y_train <- y[-folds[[fold_id]]]
      x_test <- x[folds[[fold_id]]]; y_test <- y[folds[[fold_id]]]
      valid_test_samples <- which(x_test >= min(x_train) & x_test <= max(x_train))
      x_test <- x_test[valid_test_samples]; y_test <- y_test[valid_test_samples]
      
      # Penalized cubic spline fit and prediction error
      spline_fit <- Fit_Penalized_Cubic_Spline(x_train, y_train, alpha)
      spline_predictions <- Cubic_Spline_Prediction(spline_fit, x_test)
      CV_prediction_error[which(alpha==alpha_grid)] <- CV_prediction_error[which(alpha==alpha_grid)] + 
        sum((y_test - spline_predictions)^2)
    }
  }
  
  # Optimal alpha parameter
  alpha_opt <- alpha_grid[which.min(CV_prediction_error)]
  
  # Find new interpolation points
  delta <- x[-1] - x[-length(x)]
  R_Matrix <- Get_R_Matrix(x)
  Q_Matrix <- Get_Q_Matrix(x)
  s <- alpha_opt*solve(alpha_opt*diag(x = 1, nrow = length(y)) +
                     (1 - alpha_opt)*t(Q_Matrix) %*% solve(R_Matrix) %*% Q_Matrix) %*% y
  gamma <- solve(R_Matrix) %*% Q_Matrix %*% s
  
  # Return interpolation points and gamma values
  return(list(x = x, s = s, delta = delta, gamma = c(0, gamma, 0), alpha_opt = alpha_opt,
              alpha_grid = alpha_grid, CV_prediction_error = CV_prediction_error))
}





