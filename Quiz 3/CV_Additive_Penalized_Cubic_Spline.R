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
source("Fit_Additive_Penalized_Cubic_Spline.R")

# Fitting penalized cubic spline function 
CV_Additive_Penalized_Cubic_Spline <- function(x, y, 
                                      alpha_grid = seq(0.5, 0.99, by = 0.01), 
                                      n_folds = 5,
                                      tolerance = 1e-6){
  
  # Check alpha grid
  if(any(alpha_grid <= 0) | any(alpha_grid >= 1))
    stop("The grid for alpha must only contain values greater than 0 and less than 1.")
  
  # Prediction Error
  CV_prediction_error <- numeric(length(alpha_grid))
  
  # Separate the data into n_folds groups to do cross validation.
  n=length(y)
  # Putting the samples in a random order then categorizing them.
  samp=sample(n)
  size=n/n_folds
  folds=sort(samp[1:size])
  for (i in 1:n_folds-1) {
    folds=cbind(folds,sort(samp[(size*i+1):(size*(i+1))]))
  }
  #for some reason it appends an extra column which is annoying.
  folds=folds[,-1]
  folds
  
  # Looping over tuning parameter
  for(alpha in alpha_grid){
    
    # Print iteration
    cat("CV for alpha = ", alpha, "\n", sep = "")
    
    # Looping over folds
    for(i in 1:n_folds){
      
      # Subsetting training and test data
      x_train <- x[-folds[,i],]; y_train <- y[-folds[,i]]
      x_test <- x[folds[,i],]; y_test <- y[folds[,i]]
      # want to make sure the test data is within range.
      valid_test_samples <- which(x_test[,1] >= min(x_train[,1]) & x_test[,1] <= max(x_train[,1])&
                                    x_test[,2] >= min(x_train[,2]) & x_test[,2] <= max(x_train[,2]))
      x_test <- x_test[valid_test_samples,]; y_test <- y_test[valid_test_samples]
      
      # Penalized cubic spline fit and prediction error
      spline_fit <- Fit_Additive_Penalized_Cubic_Spline(x_train, y_train, alpha, tolerance)
      spline_predictions <- Additive_Cubic_Spline_Prediction(spline_fit, x_test)
      #Keep track of the amount of error between the points.
      CV_prediction_error[which(alpha==alpha_grid)] <- CV_prediction_error[which(alpha==alpha_grid)] + 
        sum((y_test - spline_predictions)^2)
    }
  }
  
  # Optimal alpha parameter
  alpha_opt <- alpha_grid[which.min(CV_prediction_error)]
  
  return(list(alpha_opt = alpha_opt,
         alpha_grid = alpha_grid, 
         CV_prediction_error = CV_prediction_error))
}





