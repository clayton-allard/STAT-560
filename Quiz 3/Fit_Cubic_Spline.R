# ================================================
# Implementation - Cublic Splines
# 
# STAT 460/560 - Statistical Theory I
# 
# Anthony-A. Christidis
# Department of Statistics
# University of British Columbia
# ================================================

# R Matrix
Get_R_Matrix <- function(x){
  
  delta <- x[-1] - x[-length(x)]
  R_Matrix <- matrix(0, nrow = length(delta) - 1, ncol = length(delta) - 1)
  for(row_id in 2:(length(delta) - 2))
    R_Matrix[row_id, c(row_id - 1, row_id, row_id + 1)] <- c(delta[row_id], 
                                                             2*(delta[row_id] + delta[row_id + 1]),
                                                             delta[row_id + 1])
  R_Matrix[1, c(1:2)] <- c(2*(delta[1] + delta[2]), 
                           delta[2])
  R_Matrix[nrow(R_Matrix), c(nrow(R_Matrix) - 1, nrow(R_Matrix))] <- c(delta[length(delta) - 1], 
                           2*(delta[length(delta) - 1] + delta[length(delta)]))
  
  R_Matrix <- R_Matrix/6
    
  return(R_Matrix)
}

# Q Matrix
Get_Q_Matrix <- function(x){
  
  delta <- x[-1] - x[-length(x)]
  Q_Matrix <- matrix(0, nrow = length(delta) - 1, ncol = length(delta) + 1)
  for(row_id in 1:(length(delta) - 1))
    Q_Matrix[row_id, c(row_id, row_id + 1, row_id + 2)] <- c(1/delta[row_id], 
                                                             -1/delta[row_id] - 1/delta[row_id + 1],
                                                             1/delta[row_id + 1])
  return(Q_Matrix)
}

# Fitting cubic spline function (no penalization)
Fit_Cubic_Spline <- function(x, y){
  
  # Sort the data
  data_order <- order(x)
  x <- x[data_order]
  y <- y[data_order]
  
  delta <- x[-1] - x[-length(x)]
  R_Matrix <- Get_R_Matrix(x)
  Q_Matrix <- Get_Q_Matrix(x)
  gamma <- solve(R_Matrix) %*% Q_Matrix %*% y
  
  # Return interpolation points and gamma values
  return(list(x = x, s = y, delta = delta, gamma = c(0, gamma, 0)))
}

# Predict value from cubic spline
Cubic_Spline_Prediction <- function(spline, x_values){
  
  if(any(x_values < min(spline$x)) || any(x_values > max(spline$x)))
    stop("Outside of range of spline.")
  
  # # Sort the data
  # data_order <- order(x_values)
  # x_values <- x_values[data_order]

  predictions <- numeric(length(x_values))
  
  for(x_val in x_values){
    
    x_ind <- ifelse(x_val==max(spline$x), which(x_values==max(spline$x))[1], min(which(x_val < spline$x)) - 1)
    predictions[which(x_val==x_values)] <- 
      spline$gamma[x_ind]/(6*spline$delta[x_ind])*(spline$x[x_ind + 1]- x_val)^3 + 
      spline$gamma[x_ind + 1]/(6*spline$delta[x_ind])*(x_val - spline$x[x_ind])^3 + 
      (spline$s[x_ind + 1]/spline$delta[x_ind] - spline$gamma[x_ind + 1]/6*spline$delta[x_ind])*(x_val - spline$x[x_ind]) + 
      (spline$s[x_ind]/spline$delta[x_ind] - spline$gamma[x_ind]/6*spline$delta[x_ind])*(spline$x[x_ind + 1] - x_val)
  }
  
  max_ind <- which(x_values==max(spline$x))
  predictions[max_ind] <- spline$s[max_ind]
  
  # # Sort back to original data order
  # predictions <- predictions[data_order]
  
  return(predictions)
}


