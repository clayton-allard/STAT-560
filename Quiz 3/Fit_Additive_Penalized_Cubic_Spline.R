# ================================================
# Implementation - Fit_Additive_Penalized_Cubic_Spline
# 
# STAT 460/560 - Statistical Theory I
# 
# Clayton Allard
# Credit: Anthony Christidis
# Department of Statistics
# University of British Columbia
# ================================================

# Fit the penalized additive spline by creating two splines and adding them.
Fit_Additive_Penalized_Cubic_Spline <-
  function(x, y, alpha, tolerance) {
    #separating x1 and x2.
    x1 = x[, 1]
    x2 = x[, 2]
    
    #Must keep consistent ordering for everything so 
    #that no data points get misplaced.
    
    #need to order x2 so that we can get the R and Q matrix.
    data_order2 <- order(x2)
    x1 <- x1[data_order2]
    x2 <- x2[data_order2]
    y <- y[data_order2]
    # Find new interpolation points
    R_Matrix2 <- Get_R_Matrix(x2)
    Q_Matrix2 <- Get_Q_Matrix(x2)
    
    #need to order x2 so that we can get the R and Q matrix.
    data_order1 <- order(x1)
    x1 <- x1[data_order1]
    x2 <- x2[data_order1]
    y <- y[data_order1]
    # Find new interpolation points
    R_Matrix1 <- Get_R_Matrix(x1)
    Q_Matrix1 <- Get_Q_Matrix(x1)
    
    # Initializing splines.
    s1 <- rep(0, length(y))
    s2 <- rep(1, length(y))
    
    #these are the spline values for each iteration.
    m <- cbind(rep(0, length(y)), s1 + s2)
    iter = 1
    
    while (max(abs(m[, iter + 1] - m[, iter])) >= tolerance) {
      #sorting by x2 so that we can do the spline step for s2.
      data_order2 <- order(x2)
      x1 <- x1[data_order2]
      x2 <- x2[data_order2]
      s1 <- s1[data_order2]
      s2 <- s2[data_order2]
      y <- y[data_order2]
      delta2 <- x2[-1] - x2[-length(y)]
      #same formula as with the single spline, but we subtract s1.
      s2 <- alpha * solve(
        alpha * diag(x = 1, nrow = length(y)) +
          (1 - alpha) * t(Q_Matrix2) %*% solve(R_Matrix2) %*% Q_Matrix2
      ) %*% (y - s1)
      gamma2 <- c(0, solve(R_Matrix2) %*% Q_Matrix2 %*% s2, 0)
      
      #sorting by x1 so that we can do the spline step for s1.
      data_order1 <- order(x1)
      x1 <- x1[data_order1]
      x2 <- x2[data_order1]
      s1 <- s1[data_order1]
      s2 <- s2[data_order1]
      y <- y[data_order1]
      delta1 <- x1[-1] - x1[-length(y)]
      #same formula as with the single spline, but we subtract s2.
      s1 <- alpha * solve(
        alpha * diag(x = 1, nrow = length(y)) +
          (1 - alpha) * t(Q_Matrix1) %*% solve(R_Matrix1) %*% Q_Matrix1
      ) %*% (y - s2)
      gamma1 <- c(0, solve(R_Matrix1) %*% Q_Matrix1 %*% s1, 0)
      
      m <- cbind(m, s1 + s2)
      iter = iter + 1
    }
    cat("Complete after", iter - 1, "iterations\n")
    # Return interpolation points and gamma values
    return(cbind(
      list(
        x = x1,
        s = s1[, 1],
        delta = delta1,
        gamma = gamma1
      ),
      list(
        x = x2,
        s = s2,
        delta = delta2,
        gamma = gamma2
      )
    ))
  }


#Compute the function value of the additive spline.
Additive_Cubic_Spline_Prediction <- function(splines, x_values) {
  #Just need to calculate each separately then add them.
  predictions1 = Cubic_Spline_Prediction(splines[, 1], x_values[, 1])
  predictions2 = Cubic_Spline_Prediction(splines[, 2], x_values[, 2])
  
  return(predictions1 + predictions2)
}


# Predict value from a one dimensional cubic spline. This is a helper function.
Cubic_Spline_Prediction <- function(spline, x_values) {
  if (any(x_values < min(spline$x)) || any(x_values > max(spline$x)))
    stop("Outside of range of spline.")
  
  # Sort the data
  spline_order <- order(spline$x)
  x <- spline$x[spline_order]
  s <- spline$s[spline_order]
  
  # Keep track of the predictions.
  predictions <- numeric(length(x_values))
  
  for (x_val in x_values) {
    n = length(x_values)
    #remove the maximum since it does not have its own corresponding polynomial.
    k = x[-which(x == max(x))]
    #find the corresponding polynomial.
    x_ind = which(k == max(k[x_val >= k]))
    #Store the prediction for all the applicable values.
    predictions[which(x_val == x_values)] <-
      spline$gamma[x_ind] / (6 * spline$delta[x_ind]) * (x[x_ind + 1] - x_val) ^
      3 +
      spline$gamma[x_ind + 1] / (6 * spline$delta[x_ind]) * (x_val - x[x_ind]) ^
      3 +
      (s[x_ind + 1] / spline$delta[x_ind] - spline$gamma[x_ind + 1] / 6 *
         spline$delta[x_ind]) * (x_val - x[x_ind]) +
      (s[x_ind] / spline$delta[x_ind] - spline$gamma[x_ind] / 6 * spline$delta[x_ind]) *
      (x[x_ind + 1] - x_val)
  }
  
  return(predictions)
}


# R Matrix. Copied from Anthony's code. This is a helper function.
Get_R_Matrix <- function(x) {
  delta <- x[-1] - x[-length(x)]
  R_Matrix <-
    matrix(0, nrow = length(delta) - 1, ncol = length(delta) - 1)
  for (row_id in 2:(length(delta) - 2))
    R_Matrix[row_id, c(row_id - 1, row_id, row_id + 1)] <-
    c(delta[row_id],
      2 * (delta[row_id] + delta[row_id + 1]),
      delta[row_id + 1])
  R_Matrix[1, c(1:2)] <- c(2 * (delta[1] + delta[2]),
                           delta[2])
  R_Matrix[nrow(R_Matrix), c(nrow(R_Matrix) - 1, nrow(R_Matrix))] <-
    c(delta[length(delta) - 1],
      2 *
        (delta[length(delta) - 1] + delta[length(delta)]))
  
  R_Matrix <- R_Matrix / 6
  
  return(R_Matrix)
}


# Q Matrix. Copied from Anthony's code. This is a helper function.
Get_Q_Matrix <- function(x) {
  delta <- x[-1] - x[-length(x)]
  Q_Matrix <-
    matrix(0, nrow = length(delta) - 1, ncol = length(delta) + 1)
  for (row_id in 1:(length(delta) - 1))
    Q_Matrix[row_id, c(row_id, row_id + 1, row_id + 2)] <-
    c(1 / delta[row_id],-1 / delta[row_id] - 1 / delta[row_id + 1],
      1 / delta[row_id + 1])
  return(Q_Matrix)
}